#!/usr/bin/env python3
"""FPSim2 similarity index server for the NovoMCP Open Corpus.

Serves the ``molecule-index`` contract the NovoMCP engine already proxies to
when ``NOVOMCP_MOLECULE_INDEX_URL`` is set (engine config maps the service name
``molecule-index`` -> that URL). The engine forwards these results through
unchanged, so the response shapes here ARE the tool output the sideloads
render — top-level ``results`` with per-item ``cid/smiles/similarity/
molecular_weight/xlogp/qed`` (+ raw ``has_pains``). Do not rename keys.

Endpoints match the engine executors verbatim:
  * POST /api/search/similar   (search_similar)   -> exact-Tanimoto over the corpus
  * POST /api/search/filter    (filter_molecules) -> property scan over the sidecar
  * GET  /health

Honesty rules (open-corpus-index-spec.md §4c):
  * No compliance verdict. No compliance provider without the closed
    compliance layer, so that column stays blank — never a fake verdict.
    ``has_pains`` is a raw structural-alert field, not "compliance".
  * ``exclude_controlled`` / ``exclude_flagged`` are no-ops — the open corpus
    has no controlled/flagged columns. Accepted, not applied, said so in
    ``notes``. Never silently pretend to filter.
  * A SMILES that does not parse -> 400 with a clear message, not empty 200.

Config via env:
  INDEX_DIR   dir holding index.h5 / meta.parquet / manifest.json (default /data/index)
  INDEX_MODE  "in-memory" (fast, ~31 GB RAM at 122M) or "on-disk" (low RAM, ~30 s/query)
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import pandas as pd
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel, ConfigDict, Field
from rdkit import Chem

from FPSim2 import FPSim2Engine

INDEX_DIR = Path(os.getenv("INDEX_DIR", "/data/index"))
INDEX_MODE = os.getenv("INDEX_MODE", "in-memory").lower()
_IN_MEMORY = INDEX_MODE != "on-disk"

app = FastAPI(title="NovoMCP similarity index", version="1.0.0")

_engine: FPSim2Engine | None = None
_meta: pd.DataFrame | None = None
_manifest: dict = {}


class SimilarRequest(BaseModel):
    smiles: str
    limit: int = Field(10, ge=1, le=100)
    threshold: float = Field(0.7, ge=0.0, le=1.0)
    exclude_controlled: bool = False  # accepted, never applied (no-op)
    exclude_flagged: bool = False      # accepted, never applied (no-op)


class FilterRequest(BaseModel):
    # The engine flattens its `filters` dict into the top level, so accept the
    # known property bounds explicitly and tolerate extras.
    model_config = ConfigDict(extra="allow")
    mw_min: float | None = None
    mw_max: float | None = None
    qed_min: float | None = None
    qed_max: float | None = None
    logp_min: float | None = None
    logp_max: float | None = None
    tpsa_min: float | None = None
    tpsa_max: float | None = None
    limit: int = Field(10, ge=1, le=100)
    offset: int = Field(0, ge=0)
    exclude_controlled: bool = False
    exclude_flagged: bool = False


@app.on_event("startup")
def _load() -> None:
    global _engine, _meta, _manifest
    h5 = INDEX_DIR / "index.h5"
    if not h5.exists():
        raise RuntimeError(f"index not found at {h5} — run build_index.py first")
    _engine = FPSim2Engine(str(h5), in_memory_fps=_IN_MEMORY)
    _meta = pd.read_parquet(INDEX_DIR / "meta.parquet").set_index("cid")
    mf = INDEX_DIR / "manifest.json"
    if mf.exists():
        _manifest = json.loads(mf.read_text())


@app.get("/health")
def health() -> dict:
    return {
        "status": "ok" if _engine is not None else "loading",
        "n_molecules": _manifest.get("n_molecules", 0 if _meta is None else len(_meta)),
        "mode": "in-memory" if _IN_MEMORY else "on-disk",
        "fp": _manifest.get("fp", "morgan-r2-2048"),
    }


def _row(cid: int, m, similarity: float | None) -> dict:
    """Build one result item in the shape the sideloads render."""
    item = {
        "cid": cid,
        "smiles": m.get("smiles"),
        "molecular_weight": _f(m.get("molecular_weight")),
        # Two provenance-distinct logP columns, never coalesced: xlogp is
        # PubChem XLogP3 (null-honest where PubChem has no value); logp is
        # RDKit Crippen (populated). Both carried; the renderer aliases
        # xlogp -> logp, so a null xlogp falls through to the populated logp.
        "xlogp": _f(m.get("xlogp")),
        "logp": _f(m.get("logp")),
        "qed": _f(m.get("qed")),
        "has_pains": bool(m.get("has_pains")) if m.get("has_pains") is not None else None,
        # No compliance field by design — see honesty rules.
    }
    if similarity is not None:
        item["similarity"] = round(float(similarity), 4)
    return item


@app.post("/api/search/similar")
def search_similar(req: SimilarRequest) -> dict:
    if _engine is None or _meta is None:
        raise HTTPException(status_code=503, detail="index still loading")
    if Chem.MolFromSmiles(req.smiles) is None:
        raise HTTPException(status_code=400, detail=f"could not parse SMILES: {req.smiles!r}")

    # Exact Tanimoto over the whole FP set. The call differs by mode: an on-disk
    # engine raises "FPs not loaded into memory" if you call .similarity().
    if _IN_MEMORY:
        hits = _engine.similarity(req.smiles, threshold=req.threshold)
    else:
        hits = _engine.on_disk_similarity(req.smiles, threshold=req.threshold)

    notes: list[str] = []
    if req.exclude_controlled or req.exclude_flagged:
        notes.append("exclude_controlled/exclude_flagged ignored: no compliance provider in this index")

    results: list[dict] = []
    for row in hits:  # structured array: mol_id (=CID), coeff (=Tanimoto), sorted desc
        cid = int(row["mol_id"])
        try:
            results.append(_row(cid, _meta.loc[cid], float(row["coeff"])))
        except KeyError:
            continue
        if len(results) >= req.limit:
            break

    return {"query_smiles": req.smiles, "count": len(results), "results": results, "notes": notes}


@app.post("/api/search/filter")
def filter_molecules(req: FilterRequest) -> dict:
    if _meta is None:
        raise HTTPException(status_code=503, detail="index still loading")

    df = _meta
    mask = pd.Series(True, index=df.index)
    if req.mw_min is not None:
        mask &= df["molecular_weight"] >= req.mw_min
    if req.mw_max is not None:
        mask &= df["molecular_weight"] <= req.mw_max
    if req.qed_min is not None:
        mask &= df["qed"] >= req.qed_min
    if req.qed_max is not None:
        mask &= df["qed"] <= req.qed_max
    # Filter on the COMPLETE RDKit `logp`, not the sparse PubChem `xlogp` —
    # filtering on xlogp would silently drop every molecule PubChem left null,
    # returning a biased subset with no error.
    if req.logp_min is not None:
        mask &= df["logp"] >= req.logp_min
    if req.logp_max is not None:
        mask &= df["logp"] <= req.logp_max
    if req.tpsa_min is not None:
        mask &= df["tpsa"] >= req.tpsa_min
    if req.tpsa_max is not None:
        mask &= df["tpsa"] <= req.tpsa_max

    notes: list[str] = []
    if req.exclude_controlled or req.exclude_flagged:
        notes.append("exclude_controlled/exclude_flagged ignored: no compliance provider in this index")

    hits = df[mask]
    total = int(len(hits))
    page = hits.iloc[req.offset : req.offset + req.limit]
    results = [_row(int(cid), page.loc[cid], None) for cid in page.index]
    return {"count": len(results), "total_matches": total, "results": results, "notes": notes}


def _f(v) -> float | None:
    """Coerce a possibly-NaN/None meta value to a JSON-safe float or None."""
    if v is None:
        return None
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return None if f != f else round(f, 4)  # drop NaN
