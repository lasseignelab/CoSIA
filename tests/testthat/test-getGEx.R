# getGEx fetches data from ExperimentHub and requires network on first run.
# After the initial download, ExperimentHub caches locally.

test_that("getGEx returns a CoSIAn object", {
    skip_if_no_experimenthub()
    obj <- make_cosian(
        o_species = c("h_sapiens", "m_musculus"),
        output_ids = "Ensembl_id",
        map_tissues = "heart",
        map_species = "m_musculus"
    )
    obj <- getConversion(obj)
    result <- getGEx(obj)
    expect_s4_class(result, "CoSIAn")
})

test_that("getGEx populates the gex slot", {
    skip_if_no_experimenthub()
    obj <- make_cosian(
        o_species = c("h_sapiens", "m_musculus"),
        output_ids = "Ensembl_id",
        map_tissues = "heart",
        map_species = "m_musculus"
    )
    obj <- getConversion(obj)
    result <- getGEx(obj)
    expect_s3_class(result@gex, "data.frame")
    expect_gt(nrow(result@gex), 0)
})

test_that("getGEx gex slot has expected columns", {
    skip_if_no_experimenthub()
    obj <- make_cosian(
        o_species = c("h_sapiens", "m_musculus"),
        output_ids = "Ensembl_id",
        map_tissues = "heart",
        map_species = "m_musculus"
    )
    obj <- getConversion(obj)
    result <- getGEx(obj)
    expected_cols <- c(
        "Anatomical_entity_name", "Ensembl_ID", "Sample_size",
        "VST", "Experiment_ID", "Anatomical_entity_ID", "Species"
    )
    expect_true(all(expected_cols %in% colnames(result@gex)))
})

test_that("getGEx filters results to requested tissues", {
    skip_if_no_experimenthub()
    obj <- make_cosian(
        o_species = c("h_sapiens", "m_musculus"),
        output_ids = "Ensembl_id",
        map_tissues = "heart",
        map_species = "m_musculus"
    )
    obj <- getConversion(obj)
    result <- getGEx(obj)
    expect_true(all(result@gex$Anatomical_entity_name == "heart"))
})

test_that("getGEx filters results to requested species", {
    skip_if_no_experimenthub()
    obj <- make_cosian(
        o_species = c("h_sapiens", "m_musculus"),
        output_ids = "Ensembl_id",
        map_tissues = "heart",
        map_species = "m_musculus"
    )
    obj <- getConversion(obj)
    result <- getGEx(obj)
    expect_true(all(result@gex$Species == "Mus_musculus"))
})
