# validation/03_compare.R
#
# Run this script THIRD to compare Bioconductor vs dev outputs.
#
# Usage (from repo root):
#   Rscript validation/03_compare.R
#
# Requires both validation/output/bioc_results.rds and
# validation/output/dev_results.rds to exist.

stopifnot(
    "Run 01_run_bioc.R first" = file.exists("validation/output/bioc_results.rds"),
    "Run 02_run_dev.R first"  = file.exists("validation/output/dev_results.rds")
)

bioc <- readRDS("validation/output/bioc_results.rds")
dev  <- readRDS("validation/output/dev_results.rds")

pass <- 0L
fail <- 0L

check <- function(label, expr) {
    result <- tryCatch(expr, error = function(e) FALSE)
    if (isTRUE(result)) {
        cat("  PASS:", label, "\n")
        pass <<- pass + 1L
    } else {
        cat("  FAIL:", label, "\n")
        fail <<- fail + 1L
    }
}

structural_checks <- function(df, label) {
    cat("\n  Structural checks for", label, "\n")
    check("is a data.frame",         is.data.frame(df))
    check("has rows",                nrow(df) > 0)
    check("has Anatomical_entity_name column",
          "Anatomical_entity_name" %in% colnames(df))
    check("has at least one CV_tissue column",
          any(grepl("CV_[Tt]issue", colnames(df))))
    check("no duplicate rows",       nrow(df) == nrow(unique(df)))
    check("no duplicate col names",  !anyDuplicated(colnames(df)))

    cv_cols <- grep("CV_[Tt]issue", colnames(df), value = TRUE)
    for (col in cv_cols) {
        vals <- df[[col]][!is.na(df[[col]])]
        check(paste("CV values non-negative in", col), all(vals >= 0))
    }
}

# ---------------------------------------------------------------------------
cat("========================================\n")
cat("CASE 1: Single species (h_sapiens)\n")
cat("Expected: dev result MUST match Bioconductor exactly.\n")
cat("========================================\n")

b <- bioc$cv_tissue_single
d <- dev$cv_tissue_single

if (is.null(b)) {
    cat("  SKIP: Bioconductor version errored — cannot compare.\n")
} else if (is.null(d)) {
    cat("  FAIL: Dev version errored on single-species case.\n")
    fail <- fail + 1L
} else {
    structural_checks(d, "dev single-species")

    # Sort both the same way before comparing values
    sort_df <- function(df) {
        df[do.call(order, as.list(df[, seq_len(min(2, ncol(df)))])), ]
    }
    b_sorted <- sort_df(b)
    d_sorted <- sort_df(d)

    cat("\n  Value comparison (single-species must match exactly):\n")
    eq <- all.equal(b_sorted, d_sorted, check.attributes = FALSE)
    if (isTRUE(eq)) {
        cat("  PASS: dev result is identical to Bioconductor result.\n")
        pass <- pass + 1L
    } else {
        cat("  FAIL: differences found:\n")
        cat("   ", paste(eq, collapse = "\n    "), "\n")
        fail <- fail + 1L
    }
}

# ---------------------------------------------------------------------------
cat("\n========================================\n")
cat("CASE 2: Two species (h_sapiens + m_musculus)\n")
cat("Expected: Bioconductor may crash; dev should produce valid output.\n")
cat("========================================\n")

b2 <- bioc$cv_tissue_two
d2 <- dev$cv_tissue_two

if (is.null(b2)) {
    cat("  NOTE: Bioconductor version crashed (expected — known bug).\n")
} else {
    cat("  NOTE: Bioconductor version did not crash.\n")
    structural_checks(b2, "bioc two-species")
}

if (is.null(d2)) {
    cat("  FAIL: Dev version errored on two-species case.\n")
    fail <- fail + 1L
} else {
    structural_checks(d2, "dev two-species")

    if (!is.null(b2)) {
        cat("\n  Value comparison (two-species):\n")
        eq <- all.equal(b2, d2, check.attributes = FALSE)
        if (isTRUE(eq)) {
            cat("  PASS: results are identical.\n")
            pass <- pass + 1L
        } else {
            cat("  INFO: results differ (may be expected if old code was wrong):\n")
            cat("   ", paste(eq, collapse = "\n    "), "\n")
        }
    }
}

# ---------------------------------------------------------------------------
cat("\n========================================\n")
cat("CASE 3: Three species (h_sapiens + m_musculus + r_norvegicus)\n")
cat("Expected: Bioconductor crashes; dev produces valid output.\n")
cat("========================================\n")

b3 <- bioc$cv_tissue_three
d3 <- dev$cv_tissue_three

if (is.null(b3)) {
    cat("  NOTE: Bioconductor version crashed (expected — known bug).\n")
} else {
    cat("  NOTE: Bioconductor version did not crash.\n")
    structural_checks(b3, "bioc three-species")
}

if (is.null(d3)) {
    cat("  FAIL: Dev version errored on three-species case.\n")
    fail <- fail + 1L
} else {
    structural_checks(d3, "dev three-species")
}

# ---------------------------------------------------------------------------
cat("\n========================================\n")
cat("Manual spot-check instructions\n")
cat("========================================\n")
cat("
To verify a CV value by hand:

  1. Load the ExperimentHub data for the species of interest:
       eh <- ExperimentHub::ExperimentHub()
       hs_data <- eh[['EH7858']][[1]]   # h_sapiens

  2. Filter to one gene + one tissue:
       gene_data <- dplyr::filter(hs_data,
           Ensembl_ID == '<ensembl_id>',
           Anatomical_entity_name == 'heart')

  3. Compute CV manually (sd / median):
       vst <- as.numeric(unlist(gene_data$VST))
       manual_cv <- sd(vst) / median(vst)

  4. Compare to the value in dev$cv_tissue_single or dev$cv_tissue_three.

ExperimentHub IDs:
  h_sapiens       EH7858
  m_musculus      EH7859
  r_norvegicus    EH7860
  d_rerio         EH7861
  d_melanogaster  EH7862
  c_elegans       EH7863
")

# ---------------------------------------------------------------------------
cat("========================================\n")
cat(sprintf("Summary: %d passed, %d failed\n", pass, fail))
cat("========================================\n")

if (fail > 0L) quit(status = 1L)
