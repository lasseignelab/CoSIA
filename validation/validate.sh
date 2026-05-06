#!/usr/bin/env bash
set -euo pipefail

PROJECT_PATH="/data/user/$USER/CoSIA"
DEV_PATH="/data/user/$USER/CoSIA"
CONTAINER_PATH="$PROJECT_PATH/bin/container"

CACHE_BASE="/data/user/$USER/bioconductor/bioc-cache"
TMP_BASE="/data/user/$USER/bioconductor/bioc-tmp"

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
  -B "${DEV_PATH}" \
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
  "${CONTAINER_PATH}/bioc_cosia_1.10.0.sif" \
  Rscript validation/01_run_bioc.R

singularity exec --cleanenv \
  --containall \
  --home "${CACHE_BASE}/home:/home/${USER}" \
  -B "${DEV_PATH}" \
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
  "${CONTAINER_PATH}/bioc_cosia_1.10.0.sif" \
  Rscript validation/02_run_dev.R

singularity exec --cleanenv \
  --containall \
  --home "${CACHE_BASE}/home:/home/${USER}" \
  -B "${DEV_PATH}" \
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
  "${CONTAINER_PATH}/bioc_cosia_1.10.0.sif" \
  Rscript validation/03_compare.R
