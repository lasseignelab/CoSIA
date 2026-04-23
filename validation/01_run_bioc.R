# validation/01_run_bioc.R
#
# Run this script FIRST using the unmodified Bioconductor release version of
# CoSIA to capture reference outputs.
#
# Usage (from repo root, with Bioconductor CoSIA installed):
#   Rscript validation/01_run_bioc.R
#
# Output: validation/output/bioc_results.rds

stopifnot(
    "CoSIA must be installed from Bioconductor, not loaded from source.
     Run: BiocManager::install('CoSIA') in a clean R session first." =
        "CoSIA" %in% installed.packages()[, "Package"]
)

library(CoSIA)

cat("CoSIA version:", as.character(packageVersion("CoSIA")), "\n")
cat("Running with Bioconductor release version.\n\n")

dir.create("validation/output", showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Shared inputs
# Three kidney genes used in the CoSIA vignette (human Ensembl IDs)
# ---------------------------------------------------------------------------
GENE_SET   <- c("ENSG00000008710", "ENSG00000118762", "ENSG00000152217")
I_SPECIES  <- "h_sapiens"
INPUT_ID   <- "Ensembl_id"
OUTPUT_IDS <- c("Ensembl_id", "Symbol")
TOOL       <- "annotationDBI"
ORTHO_DB   <- "HomoloGene"
TISSUE     <- "heart"

run_cv_tissue <- function(o_species, map_species, label) {
    cat("--- Running CV_Tissue:", label, "---\n")
    obj <- CoSIAn(
        gene_set         = GENE_SET,
        i_species        = I_SPECIES,
        input_id         = INPUT_ID,
        o_species        = o_species,
        output_ids       = OUTPUT_IDS,
        mapping_tool     = TOOL,
        ortholog_database = ORTHO_DB,
        map_tissues      = TISSUE,
        map_species      = map_species,
        metric_type      = "CV_Tissue"
    )
    obj <- getConversion(obj)
    obj <- getGExMetrics(obj)
    metric <- obj@metric
    cat("  Rows:", nrow(metric), " Cols:", ncol(metric), "\n")
    cat("  Columns:", paste(colnames(metric), collapse = ", "), "\n\n")
    metric
}

results <- list(

    # Case 1: single species — the code path that was unaffected by the
    # merge bug; results here MUST match the dev version exactly.
    cv_tissue_single = run_cv_tissue(
        o_species  = "h_sapiens",
        map_species = "h_sapiens",
        label      = "single species (h_sapiens)"
    ),

    # Case 2: two species — simplest multi-species case affected by the bug.
    cv_tissue_two = tryCatch(
        run_cv_tissue(
            o_species  = c("h_sapiens", "m_musculus"),
            map_species = c("h_sapiens", "m_musculus"),
            label      = "two species (h_sapiens + m_musculus)"
        ),
        error = function(e) {
            cat("  ERROR:", conditionMessage(e), "\n\n")
            NULL
        }
    ),

    # Case 3: three species — the case shown in the bug report.
    cv_tissue_three = tryCatch(
        run_cv_tissue(
            o_species  = c("h_sapiens", "m_musculus", "r_norvegicus"),
            map_species = c("h_sapiens", "m_musculus", "r_norvegicus"),
            label      = "three species (h_sapiens + m_musculus + r_norvegicus)"
        ),
        error = function(e) {
            cat("  ERROR:", conditionMessage(e), "\n\n")
            NULL
        }
    )
)

saveRDS(results, "validation/output/bioc_results.rds")
cat("Saved to validation/output/bioc_results.rds\n")
