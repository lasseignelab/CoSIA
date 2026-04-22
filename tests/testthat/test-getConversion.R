# getConversion uses AnnotationDbi (local org.*.eg.db packages) and HomoloGene
# (bundled data), so these tests do not require network access.

test_that("getConversion returns a CoSIAn object", {
    obj <- make_cosian(o_species = "h_sapiens", map_species = "h_sapiens")
    result <- getConversion(obj)
    expect_s4_class(result, "CoSIAn")
})

test_that("getConversion populates the converted_id slot", {
    obj <- make_cosian(o_species = "h_sapiens", map_species = "h_sapiens")
    result <- getConversion(obj)
    expect_s3_class(result@converted_id, "data.frame")
    expect_gt(nrow(result@converted_id), 0)
})

test_that("same-species conversion produces expected column names", {
    obj <- make_cosian(
        o_species = "h_sapiens",
        output_ids = c("Ensembl_id", "Symbol"),
        map_species = "h_sapiens"
    )
    result <- getConversion(obj)
    cols <- colnames(result@converted_id)
    expect_true(any(grepl("h_sapiens_ensembl_id", cols)))
    expect_true(any(grepl("h_sapiens_symbol", cols)))
})

test_that("cross-species conversion with HomoloGene produces output columns for both species", {
    obj <- make_cosian(
        o_species = c("h_sapiens", "m_musculus"),
        output_ids = "Ensembl_id",
        map_species = "m_musculus"
    )
    result <- getConversion(obj)
    cols <- colnames(result@converted_id)
    expect_true(any(grepl("h_sapiens_ensembl_id", cols)))
    expect_true(any(grepl("m_musculus_ensembl_id", cols)))
})

test_that("getConversion deduplicates gene_set input", {
    obj <- make_cosian(
        gene_set = c("ENSG00000008710", "ENSG00000008710"),
        o_species = "h_sapiens",
        map_species = "h_sapiens"
    )
    result <- getConversion(obj)
    hs_col <- result@converted_id[["h_sapiens_ensembl_id"]]
    expect_equal(sum(hs_col == "ENSG00000008710", na.rm = TRUE), 1L)
})

test_that("getConversion with NCBIOrtho requires network", {
    skip_if_offline()
    obj <- make_cosian(
        o_species = c("h_sapiens", "m_musculus"),
        output_ids = "Ensembl_id",
        ortholog_database = "NCBIOrtho",
        map_species = "m_musculus"
    )
    result <- getConversion(obj)
    expect_s4_class(result, "CoSIAn")
    expect_gt(nrow(result@converted_id), 0)
})
