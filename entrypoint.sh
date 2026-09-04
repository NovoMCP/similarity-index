#!/usr/bin/env bash
# Dispatch: `build ...` runs the indexer; `serve --index DIR` runs the server.
set -euo pipefail

cmd="${1:-serve}"
shift || true

case "$cmd" in
  build)
    exec python build_index.py "$@"
    ;;
  serve)
    # Accept `--index DIR` as sugar for INDEX_DIR; pass INDEX_MODE via env.
    while [ $# -gt 0 ]; do
      case "$1" in
        --index) export INDEX_DIR="$2"; shift 2 ;;
        --mode)  export INDEX_MODE="$2"; shift 2 ;;
        *) shift ;;
      esac
    done
    exec uvicorn server:app --host 0.0.0.0 --port "${PORT:-8080}"
    ;;
  *)
    echo "usage: [build ...] | [serve --index DIR [--mode in-memory|on-disk]]" >&2
    exit 2
    ;;
esac
