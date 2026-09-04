#!/usr/bin/env python3
"""Build an FPSim2 similarity index from the NovoMCP Open Corpus.

Reads the published open-corpus parquet (local path, glob, or an
``s3://novomcp-open-corpus/...`` anonymous glob), streams ``(smiles, cid)``
into an FPSim2 Morgan-r2-2048 fingerprint database, and writes a small
``meta.parquet`` display/filter sidecar keyed by CID plus a ``manifest.json``.

Design (per open-corpus-index-spec.md):
  * FPSim2 does the heavy exact-Tanimoto indexing; we do not hand-roll search.
  * The FPSim2 db holds fingerprints + CIDs. The sidecar holds the handful of
    display/filter fields the sideloads render.
  * meta.parquet is a DEFAULT-DENY allowlist: it emits ONLY the six fields
    below and drops every other input column, so even a richer parquet cannot
    smuggle controlled/compliance columns through. Build only from the published
    open corpus (already stripped of closed columns); never point this at an
    internal corpus.

Verified against FPSim2 0.7.4 + RDKit 2026.03: the Morgan fp param is
``fpSize`` (not the legacy ``nBits``); ``mol_format='smiles'`` with a
``[smiles, id]`` iterable passes the corpus CID straight through as the
FPSim2 mol-id, so search results map back to CID with no extra join.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Iterator

import pyarrow.dataset as ds
from FPSim2.io import create_db_file

# The ONLY columns meta.parquet may carry. Fail closed: anything not here is
# dropped. Keep in sync with the /search response contract in server.py.
#
# Two logP columns, kept DISTINCT by provenance — never coalesced into one:
#   * xlogp = PubChem XLogP3. Genuinely absent for a real slice of PubChem
#     (some elements, many charged/zwitterionic species, salts/mixtures), so
#     we carry its null faithfully rather than fabricate a value.
#   * logp  = RDKit Crippen, computed in-house, populated for ~every valid
#     molecule. This is the column property FILTERS use (see server.py), so a
#     logp filter doesn't silently drop the null-xlogp molecules.
# Filling `xlogp` from `logp` would be a provenance smear (different methods,
# same field name) — the exact thing the open corpus is built to avoid.
META_ALLOWLIST = ["cid", "smiles", "molecular_weight", "xlogp", "logp", "tpsa", "qed", "has_pains"]

# Input corpus column aliases -> our canonical names. The published corpus has
# settled column names, but tolerate the common variants so a build never
# silently drops a field to null. xlogp and logp are resolved SEPARATELY.
COLUMN_ALIASES = {
    "cid": ["cid", "CID", "pubchem_cid"],
    "smiles": ["smiles", "canonical_smiles", "SMILES"],
    "molecular_weight": ["molecular_weight", "mw", "MolWt"],
    "xlogp": ["xlogp", "XLogP", "xlogp3"],
    "logp": ["logp", "crippen_logp", "MolLogP", "rdkit_logp"],
    "tpsa": ["tpsa", "TPSA"],
    "qed": ["qed", "QED"],
    "has_pains": ["has_pains", "pains", "PAINS"],
}

DRUGLIKE_PRESET = {"mw_max": 600.0, "qed_min": 0.3}


def _resolve_columns(schema_names: list[str]) -> dict[str, str]:
    """Map our canonical field -> the actual column name present in the parquet."""
    lower = {n.lower(): n for n in schema_names}
    resolved: dict[str, str] = {}
    for canonical, candidates in COLUMN_ALIASES.items():
        for c in candidates:
            if c.lower() in lower:
                resolved[canonical] = lower[c.lower()]
                break
    for required in ("cid", "smiles"):
        if required not in resolved:
            sys.exit(f"FATAL: corpus is missing a '{required}' column (looked for {COLUMN_ALIASES[required]})")
    return resolved


def _open_dataset(input_spec: str):
    """Open a parquet dataset from a local path/glob or an s3:// glob (anon)."""
    if input_spec.startswith("s3://"):
        import pyarrow.fs as pafs

        fs = pafs.S3FileSystem(anonymous=True)
        path = input_spec[len("s3://") :]
        return ds.dataset(path, filesystem=fs, format="parquet")
    return ds.dataset(input_spec, format="parquet")


def _iter_rows(dataset, cols: dict[str, str], preset: dict | None, limit: int | None) -> Iterator[list]:
    """Yield rows as dicts of canonical field -> value, applying an optional preset."""
    read_cols = list(cols.values())
    n = 0
    for batch in dataset.to_batches(columns=read_cols):
        table = batch.to_pydict()
        # map actual column name -> list, then walk row-wise
        inv = {v: k for k, v in cols.items()}
        canonical_cols = {inv[actual]: table[actual] for actual in read_cols}
        for i in range(len(table[read_cols[0]])):
            row = {c: canonical_cols[c][i] for c in canonical_cols}
            if preset:
                mw = row.get("molecular_weight")
                qed = row.get("qed")
                if mw is not None and mw > preset["mw_max"]:
                    continue
                if qed is not None and qed < preset["qed_min"]:
                    continue
            yield row
            n += 1
            if limit and n >= limit:
                return


def main() -> None:
    p = argparse.ArgumentParser(description="Build an FPSim2 index from the NovoMCP Open Corpus.")
    p.add_argument("--input", required=True, help="parquet path, glob, or s3:// glob (anonymous read)")
    p.add_argument("--out", default="/data/index", help="output directory")
    p.add_argument("--limit", type=int, default=None, help="quick-test slice (dev/CI only) — first N molecules")
    p.add_argument("--preset", choices=["druglike"], default=None, help="druglike = MW<600 & QED>0.3")
    p.add_argument("--fp-radius", type=int, default=2)
    p.add_argument("--fp-size", type=int, default=2048)
    args = p.parse_args()

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    h5_path = out_dir / "index.h5"
    meta_path = out_dir / "meta.parquet"
    manifest_path = out_dir / "manifest.json"

    preset = DRUGLIKE_PRESET if args.preset == "druglike" else None

    dataset = _open_dataset(args.input)
    cols = _resolve_columns(dataset.schema.names)
    print(f"Resolved corpus columns: {cols}", file=sys.stderr)
    if args.limit:
        print(f"WARNING: --limit {args.limit} builds a QUICK-TEST SLICE, not the full corpus.", file=sys.stderr)

    # Stream once. We feed FPSim2 a GENERATOR of [smiles, cid] (never a
    # materialized 122M-row list) and write the meta sidecar incrementally to
    # parquet in batches, so peak RAM is FPSim2's own fingerprint set, not our
    # input. The generator side-writes meta as FPSim2 consumes it — single pass,
    # single read of the corpus.
    import pyarrow as pa
    import pyarrow.parquet as pq

    META_BATCH = 1_000_000
    meta_buf: dict[str, list] = {k: [] for k in META_ALLOWLIST}
    writer: pq.ParquetWriter | None = None
    counter = {"kept": 0}

    def _flush_meta() -> None:
        nonlocal writer
        if not meta_buf[META_ALLOWLIST[0]]:
            return
        table = pa.table(meta_buf)
        if writer is None:
            writer = pq.ParquetWriter(str(meta_path), table.schema)
        writer.write_table(table)
        for k in meta_buf:
            meta_buf[k].clear()

    def mol_gen():
        for row in _iter_rows(dataset, cols, preset, args.limit):
            smiles, cid = row.get("smiles"), row.get("cid")
            if not smiles or cid is None:
                continue
            for field in META_ALLOWLIST:
                meta_buf[field].append(row.get(field))
            counter["kept"] += 1
            if counter["kept"] % META_BATCH == 0:
                _flush_meta()
                print(f"  ...{counter['kept']:,} molecules ingested", file=sys.stderr)
            yield [smiles, int(cid)]

    # FPSim2 computes + stores fingerprints sorted by popcount for fast search,
    # keyed by the CID we pass as mol-id. fpSize (not nBits) per 0.7.4.
    print("Building FPSim2 fingerprint database (this is the slow step)...", file=sys.stderr)
    create_db_file(
        mols_source=mol_gen(),
        filename=str(h5_path),
        mol_format="smiles",
        fp_type="Morgan",
        fp_params={"radius": args.fp_radius, "fpSize": args.fp_size},
    )
    _flush_meta()
    if writer is not None:
        writer.close()

    kept = counter["kept"]
    if kept == 0:
        sys.exit("FATAL: no molecules ingested — check --input and column names.")
    print(f"Ingested {kept:,} molecules.", file=sys.stderr)

    manifest = {
        "n_molecules": kept,
        "fp": f"morgan-r{args.fp_radius}-{args.fp_size}",
        "source": args.input,
        "preset": args.preset,
        "is_slice": bool(args.limit),
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Wrote {h5_path}, {meta_path}, {manifest_path}", file=sys.stderr)
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
