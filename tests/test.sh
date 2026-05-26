#!/usr/bin/env bash
set -euo pipefail

PROJECT_PATH="$(pwd)"
CONTAINER_PATH="$PROJECT_PATH/bin/container"
SIF="${CONTAINER_PATH}/bioc_cosia_1.10.0.sif"
SIF_URI="${COSIA_SIF_URI:-docker://lizzyr/bioc_cosia:1.10.0}"

ensure_sif() {
    mkdir -p "$CONTAINER_PATH"
    if [[ -f "$SIF" ]]; then
        return
    fi
    echo "SIF not found: $SIF"
    echo "Pulling from $SIF_URI ..."
    singularity pull "$SIF" "$SIF_URI"
}

CACHE_BASE="/data/user/$USER/bioconductor/bioc-cache"
TMP_BASE="/data/user/$USER/bioconductor/bioc-tmp"

ensure_sif

mkdir -p \
  "$CACHE_BASE/home" \
  "$CACHE_BASE/annotationhub" \
  "$CACHE_BASE/experimenthub" \
  "$CACHE_BASE/biocfilecache" \
  "$CACHE_BASE/xdg" \
  "$TMP_BASE"

singularity exec --cleanenv \
  --containall \
  --home "${CACHE_BASE}/home:/home/${USER}" \
  -B "${PROJECT_PATH}" \
  -B "${CACHE_BASE}:${CACHE_BASE}" \
  -B "${TMP_BASE}:/tmp" \
  -B "${TMP_BASE}:/var/tmp" \
  --env TMPDIR="/tmp" \
  --env TMP="/tmp" \
  --env TEMP="/tmp" \
  --env XDG_CACHE_HOME="${CACHE_BASE}/xdg" \
  --env ANNOTATION_HUB_CACHE="${CACHE_BASE}/annotationhub" \
  --env EXPERIMENT_HUB_CACHE="${CACHE_BASE}/experimenthub" \
  --env BFC_CACHE="${CACHE_BASE}/biocfilecache" \
  --pwd "${PROJECT_PATH}" \
  "$SIF" \
  Rscript -e "devtools::test()"
