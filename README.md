# NovoMCP similarity index

Exact-Tanimoto similarity search over the **NovoMCP Open Corpus** (122M PubChem
compounds), self-hosted, free. This is the reference implementation of the
`NOVOMCP_MOLECULE_INDEX_URL` service — point the NovoMCP engine at it to light
up `search_similar`, `filter_molecules`, `vector_search`, and the tree tools.

We host nothing. You build the index once from the public corpus (already on
Kaggle / Zenodo / AWS Open Data) and serve it on your own hardware.

## What it is

- **Method:** Morgan fingerprint (radius 2, 2048 bits) + **exact Tanimoto**, via
  [FPSim2](https://github.com/chembl/FPSim2) (MIT, EMBL-EBI) — popcount-bounds
  pruning over a memory-mapped fingerprint database. Deterministic: same SMILES
  in, same result out, no model drift.
- **In-memory mode** — full FP db in RAM. The flagship "full corpus, fast,
  exact" path. **Latency is threshold-dependent** (tighter thresholds prune far
  more): measured on a ~40M index on one r6i.4xlarge, median query time was
  **65 ms at threshold 0.9, 283 ms at 0.7, 791 ms at 0.5**. Sub-second at the
  thresholds people actually search (≥0.7); query time grows roughly linearly
  with corpus size, so at the full 122M a loose 0.5-threshold search is a few
  seconds, not sub-second. RAM scales with the FP set (~31 GB at 122M).
- **On-disk mode** — for machines that can't hold it in RAM; much slower per
  query, so it ships as a batch/slice path, not behind an interactive card.

## Build

```bash
# From the public corpus on S3 (anonymous read), or a local parquet path/glob:
docker run -v "$PWD/corpus:/data" novomcp/similarity-index \
  build --input 's3://novomcp-open-corpus/novomcp-open-corpus-lite/*.parquet' --out /data/index

# Laptop-friendly slice to try it out in a couple of minutes:
docker run -v "$PWD/corpus:/data" novomcp/similarity-index \
  build --input '/data/corpus/*.parquet' --out /data/index --preset druglike --limit 200000
```

`--preset druglike` keeps `MW < 600 & QED > 0.3`. `--limit N` builds a
quick-test slice (dev/CI only — not the deliverable).

Outputs `index.h5` (FPSim2 db), `meta.parquet` (display/filter sidecar), and
`manifest.json`.

## Serve

```bash
docker run -p 8080:8080 -v "$PWD/corpus:/data" novomcp/similarity-index \
  serve --index /data/index                 # add --mode on-disk for low RAM

# Point the engine at it:
export NOVOMCP_MOLECULE_INDEX_URL=http://localhost:8080
```

Endpoints match what the NovoMCP engine's `molecule-index` proxy already calls:

- `POST /api/search/similar` ← `{ "smiles": "...", "limit": 10, "threshold": 0.7 }`
  → `{ query_smiles, count, results: [{ cid, smiles, similarity, molecular_weight, xlogp, logp, qed, has_pains }], notes }`
- `POST /api/search/filter` ← flat property bounds `{ "mw_max": 500, "qed_min": 0.3, "logp_max": 5, "limit": 10, "offset": 0 }`
  → `{ count, total_matches, results: [{ cid, smiles, molecular_weight, xlogp, logp, qed, has_pains }], notes }`

**Provenance — every column records its source, blanks stay blank.** `logp`
traces to **RDKit Crippen** and is complete; `xlogp` carries **PubChem XLogP3**
and is left null where XLogP3 is absent (~18% of the corpus). Filters run on the
complete `logp`; both columns are emitted, neither is coalesced. (Same rule the
open corpus is built on: the compliance column is likewise left blank where no
open source exists.)
- `GET /health` → mode + molecule count.

Point the engine at this service with `NOVOMCP_MOLECULE_INDEX_URL`; it's already
wired to route the `search_similar` / `filter_molecules` tools here.

## Honesty rules

- **No compliance verdict.** This stack has no open compliance provider, so the
  compliance column is left blank — never a fake certified verdict. `has_pains`
  is a raw structural-alert field, explicitly not "compliance".
- **`exclude_controlled` is a no-op.** The open corpus has no controlled-
  substance flags (those live in the closed compliance layer). The param is
  accepted, not applied, and the response says so in `notes`.
- **Default-deny fields.** The index is built only from the published,
  already-stripped open corpus, and `meta.parquet` emits only a fixed allowlist
  of display/filter fields — every other input column is dropped.

## License

The service code is Apache-2.0 (see `LICENSE`). The corpus it indexes is the
NovoMCP Open Corpus, published under **CC-BY-4.0** — you provide it. FPSim2 is
MIT (EMBL-EBI).
