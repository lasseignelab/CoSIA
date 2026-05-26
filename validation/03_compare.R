# validation/03_compare.R
#
# Compares production (1.11.1) vs fixed dev (1.11.2) outputs.
#
# Strategy:
#   - Metrics that work in 1.11.1 (CV_Species, DS_*): values must match exactly
#   - CV_Tissue: 1.11.1 crashes (confirms bug), 1.11.2 succeeds (confirms fix)
#     and is validated structurally + by manual spot-check
#
# Usage (from repo root, after running 01 and 02):
#   Rscript validation/03_compare.R

stopifnot(
    "Run 01_run_bioc.R first" = file.exists("validation/output/bioc_results.rds"),
    "Run 02_run_dev.R first"  = file.exists("validation/output/dev_results.rds")
)

bioc <- readRDS("validation/output/bioc_results.rds")
dev  <- readRDS("validation/output/dev_results.rds")

cat("Production version :", bioc$version, "\n")
cat("Dev version        :", dev$version, "\n\n")

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

compare_exact <- function(b, d, label) {
    if (is.null(b)) {
        cat("  SKIP:", label, "— production version errored, no reference\n")
        return()
    }
    if (is.null(d)) {
        cat("  FAIL:", label, "— dev version errored\n")
        fail <<- fail + 1L
        return()
    }
    sort_df <- function(df) df[do.call(order, lapply(df, identity)), ]
    eq <- all.equal(sort_df(b), sort_df(d), check.attributes = FALSE)
    if (isTRUE(eq)) {
        check(paste(label, "— values identical between versions"), TRUE)
    } else {
        cat("  FAIL:", label, "— values differ:\n")
        cat("   ", paste(utils::head(eq, 5), collapse = "\n    "), "\n")
        fail <<- fail + 1L
    }
}

structural_checks <- function(df, label) {
    check(paste(label, "is a data.frame"),
          is.data.frame(df))
    check(paste(label, "has rows"),
          nrow(df) > 0)
    check(paste(label, "has Anatomical_entity_name"),
          "Anatomical_entity_name" %in% colnames(df))
    check(paste(label, "has CV_tissue column"),
          any(grepl("CV_[Tt]issue", colnames(df), ignore.case = TRUE)))
    check(paste(label, "no duplicate rows"),
          nrow(df) == nrow(unique(df)))
    check(paste(label, "no duplicate column names"),
          !anyDuplicated(colnames(df)))
    check(paste(label, "no all-NA columns"),
          !any(vapply(df, function(x) all(is.na(x)), logical(1))))
    cv_cols <- grep("CV_[Tt]issue", colnames(df), value = TRUE, ignore.case = TRUE)
    for (col in cv_cols) {
        vals <- df[[col]][!is.na(df[[col]])]
        check(paste(label, "—", col, "values non-negative"),
              length(vals) > 0 && all(vals >= 0))
    }
    n_genes   <- length(unique(df[[1]]))
    n_tissues <- length(unique(df$Anatomical_entity_name))
    cat("  INFO:", n_genes, "genes x", n_tissues, "tissues =",
        nrow(unique(df[, 1:2])), "gene-tissue rows (", nrow(df), "total)\n")
}

# ---------------------------------------------------------------------------
cat("========================================\n")
cat("PART 1: CV_Tissue — confirm bug and fix\n")
cat("========================================\n\n")

for (case in c("cv_tissue_single", "cv_tissue_multi")) {
    label <- gsub("cv_tissue_", "", case)
    b <- bioc[[case]]
    d <- dev[[case]]

    cat("--- CV_Tissue", label, "---\n")
    if (is.null(b)) {
        cat("  PASS: production version crashed (bug confirmed)\n")
        pass <- pass + 1L
    } else {
        cat("  NOTE: production version did not crash —",
            "bug may already be fixed upstream\n")
    }
    if (is.null(d)) {
        cat("  FAIL: dev version crashed — fix not working\n\n")
        fail <- fail + 1L
    } else {
        structural_checks(d, paste("dev CV_Tissue", label))
        cat("\n")
    }
}

# Cross-case consistency: h_sapiens CV values must be the same whether run
# alone or alongside other species (expression data doesn't change)
d_single <- dev$cv_tissue_single
d_multi  <- dev$cv_tissue_multi
if (!is.null(d_single) && !is.null(d_multi)) {
    hs_single <- grep("h_sapiens_CV", colnames(d_single), value = TRUE)
    hs_multi  <- grep("h_sapiens_CV", colnames(d_multi),  value = TRUE)
    id_col    <- grep("ensembl_id",   colnames(d_single), value = TRUE)[1]
    if (length(hs_single) == 1 && length(hs_multi) == 1) {
        merged <- merge(
            d_single[, c(id_col, "Anatomical_entity_name", hs_single)],
            d_multi[,  c(id_col, "Anatomical_entity_name", hs_multi)],
            by = c(id_col, "Anatomical_entity_name")
        )
        check(
            "h_sapiens CV_tissue values identical in single vs multi-species run",
            nrow(merged) > 0 && isTRUE(all.equal(merged[[3]], merged[[4]],
                                                   tolerance = 1e-10))
        )
    }
}
cat("\n")

# ---------------------------------------------------------------------------
cat("========================================\n")
cat("PART 2: Unchanged metrics must match exactly\n")
cat("========================================\n\n")

compare_exact(bioc$cv_species,   dev$cv_species,   "CV_Species")
compare_exact(bioc$ds_gene,      dev$ds_gene,      "DS_Gene")
compare_exact(bioc$ds_gene_all,  dev$ds_gene_all,  "DS_Gene_all")
compare_exact(bioc$ds_tissue,    dev$ds_tissue,    "DS_Tissue")
compare_exact(bioc$ds_tissue_all, dev$ds_tissue_all, "DS_Tissue_all")

# ---------------------------------------------------------------------------
cat("\n========================================\n")
cat("PART 3: Manual spot-check for CV_Tissue values\n")
cat("========================================\n")
cat("
To verify a CV value against raw ExperimentHub data:

  devtools::load_all('.')
  eh      <- ExperimentHub::ExperimentHub()
  hs_data <- eh[['EH7858']][[1]]

  gene_tissue <- dplyr::filter(hs_data,
      Ensembl_ID == 'ENSG00000008710',
      Anatomical_entity_name == 'heart')

  vst       <- as.numeric(unlist(gene_tissue$VST))
  manual_cv <- sd(vst) / median(vst)
  cat('Manual CV:', manual_cv, '\\n')

  # Compare to dev result:
  d <- readRDS('validation/output/dev_results.rds')
  d$cv_tissue_single[d$cv_tissue_single[[1]] == 'ENSG00000008710', ]

ExperimentHub IDs:
  h_sapiens EH7858 | m_musculus EH7859 | r_norvegicus EH7860
  d_rerio   EH7861 | d_melanogaster EH7862 | c_elegans EH7863
")

# ---------------------------------------------------------------------------
cat("========================================\n")
cat(sprintf("Summary: %d passed, %d failed\n", pass, fail))
cat("========================================\n")

if (fail > 0L) quit(status = 1L)
