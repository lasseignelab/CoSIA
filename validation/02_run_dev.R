# validation/02_run_dev.R
#
# Captures outputs from the fixed development version (1.11.2).
# All metric types including CV_Tissue are expected to succeed.
#
# Usage (from repo root):
#   Rscript validation/02_run_dev.R
#
# Output: validation/output/dev_results.rds

devtools::load_all(".")

loaded_version <- as.character(packageVersion("CoSIA"))
loaded_path    <- find.package("CoSIA")

cat("CoSIA version  :", loaded_version, "\n")
cat("CoSIA path     :", loaded_path, "\n\n")

if (loaded_version != "1.11.2") {
    warning("Expected dev version 1.11.2 but found ", loaded_version, ".")
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

    cv_species = run_metric(
        "CV_Species",
        map_species  = MULTI_SPECIES,
        map_tissues  = TISSUE,
        metric_type  = "CV_Species",
        o_species    = MULTI_SPECIES
    ),

    ds_gene = run_metric(
        "DS_Gene",
        map_species  = MULTI_SPECIES,
        map_tissues  = TISSUE,
        metric_type  = "DS_Gene",
        o_species    = MULTI_SPECIES
    ),

    ds_gene_all = run_metric(
        "DS_Gene_all",
        map_species  = MULTI_SPECIES,
        map_tissues  = TISSUE,
        metric_type  = "DS_Gene_all",
        o_species    = MULTI_SPECIES
    ),

    ds_tissue = run_metric(
        "DS_Tissue",
        map_species  = MULTI_SPECIES,
        map_tissues  = TISSUE,
        metric_type  = "DS_Tissue",
        o_species    = MULTI_SPECIES
    ),

    ds_tissue_all = run_metric(
        "DS_Tissue_all",
        map_species  = MULTI_SPECIES,
        map_tissues  = TISSUE,
        metric_type  = "DS_Tissue_all",
        o_species    = MULTI_SPECIES
    )
)

saveRDS(results, "validation/output/dev_results.rds")
cat("Saved to validation/output/dev_results.rds\n")
