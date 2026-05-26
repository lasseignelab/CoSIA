# validation/01_run_bioc.R
#
# Captures outputs from the Bioconductor production version (1.11.1) which
# contains the CV_Tissue bug. Metrics other than CV_Tissue are expected to
# succeed and provide reference values for comparison with 1.11.2.
#
# NOTE ON VERSION LOADING (Singularity)
# --------------------------------------
# If your container bind-mounts a cache home directory and you have previously
# run devtools::install(), a newer dev version may shadow the container's
# installed version. Check what is loaded:
#
#   Rscript -e "cat(find.package('CoSIA'),'\n'); cat(as.character(packageVersion('CoSIA')),'\n')"
#
# To remove a shadowing dev install from your cache home:
#   find ${CACHE_BASE}/home/R -type d -name "CoSIA" -exec rm -rf {} +
#
# Usage (from repo root, with production 1.11.1 loaded):
#   Rscript validation/01_run_bioc.R
#
# Output: validation/output/bioc_results.rds

library(CoSIA)

loaded_version <- as.character(packageVersion("CoSIA"))
loaded_path    <- find.package("CoSIA")

cat("CoSIA version  :", loaded_version, "\n")
cat("CoSIA path     :", loaded_path, "\n\n")

if (loaded_version != "1.11.1") {
    warning("Expected production version 1.11.1 but found ", loaded_version,
            ". See NOTE ON VERSION LOADING above.\n",
            "Continuing — results will document whatever version is loaded.")
}

dir.create("validation/output", showWarnings = FALSE, recursive = TRUE)

GENE_SET   <- c("ENSG00000008710", "ENSG00000118762", "ENSG00000152217")
I_SPECIES  <- "h_sapiens"
INPUT_ID   <- "Ensembl_id"
OUTPUT_IDS <- c("Ensembl_id", "Symbol")
TOOL       <- "annotationDBI"
ORTHO_DB   <- "HomoloGene"
TISSUE     <- "heart"
MULTI_SPECIES <- c("h_sapiens", "m_musculus", "r_norvegicus")

make_obj <- function(map_species, map_tissues, metric_type,
                     o_species = map_species) {
    obj <- CoSIAn(
        gene_set          = GENE_SET,
        i_species         = I_SPECIES,
        input_id          = INPUT_ID,
        o_species         = o_species,
        output_ids        = OUTPUT_IDS,
        mapping_tool      = TOOL,
        ortholog_database = ORTHO_DB,
        map_tissues       = map_tissues,
        map_species       = map_species,
        metric_type       = metric_type
    )
    getConversion(obj)
}

run_metric <- function(label, ...) {
    cat("---", label, "---\n")
    tryCatch({
        obj    <- make_obj(...)
        obj    <- getGExMetrics(obj)
        metric <- obj@metric
        cat("  SUCCESS — Rows:", nrow(metric), " Cols:", ncol(metric), "\n")
        cat("  Columns:", paste(colnames(metric), collapse = ", "), "\n\n")
        metric
    }, error = function(e) {
        cat("  ERROR:", conditionMessage(e), "\n\n")
        NULL
    })
}

results <- list(
    version = loaded_version,
    path    = loaded_path,

    # CV_Tissue — expected to crash in 1.11.1 (the known bug)
    cv_tissue_single = run_metric(
        "CV_Tissue single species",
        map_species  = I_SPECIES,
        map_tissues  = TISSUE,
        metric_type  = "CV_Tissue"
    ),
    cv_tissue_multi = run_metric(
        "CV_Tissue multi-species",
        map_species  = MULTI_SPECIES,
        map_tissues  = TISSUE,
        metric_type  = "CV_Tissue",
        o_species    = MULTI_SPECIES
    ),

    # CV_Species — expected to succeed; values used as reference
    cv_species = run_metric(
        "CV_Species",
        map_species  = MULTI_SPECIES,
        map_tissues  = TISSUE,
        metric_type  = "CV_Species",
        o_species    = MULTI_SPECIES
    ),

    # DS_Gene — expected to succeed; values used as reference
    ds_gene = run_metric(
        "DS_Gene",
        map_species  = MULTI_SPECIES,
        map_tissues  = TISSUE,
        metric_type  = "DS_Gene",
        o_species    = MULTI_SPECIES
    ),

    # DS_Gene_all — expected to succeed; values used as reference
    ds_gene_all = run_metric(
        "DS_Gene_all",
        map_species  = MULTI_SPECIES,
        map_tissues  = TISSUE,
        metric_type  = "DS_Gene_all",
        o_species    = MULTI_SPECIES
    ),

    # DS_Tissue — expected to succeed; values used as reference
    ds_tissue = run_metric(
        "DS_Tissue",
        map_species  = MULTI_SPECIES,
        map_tissues  = TISSUE,
        metric_type  = "DS_Tissue",
        o_species    = MULTI_SPECIES
    ),

    # DS_Tissue_all — expected to succeed; values used as reference
    ds_tissue_all = run_metric(
        "DS_Tissue_all",
        map_species  = MULTI_SPECIES,
        map_tissues  = TISSUE,
        metric_type  = "DS_Tissue_all",
        o_species    = MULTI_SPECIES
    )
)

saveRDS(results, "validation/output/bioc_results.rds")
cat("Saved to validation/output/bioc_results.rds\n")
