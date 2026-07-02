# getGExMetrics requires ExperimentHub data; tests skip when offline.
# Helper to build a converted CoSIAn object for a given metric type.
make_converted <- function(metric_type, map_tissues = "heart",
                           map_species = "m_musculus") {
    obj <- make_cosian(
        o_species = c("h_sapiens", "m_musculus"),
        output_ids = "Ensembl_id",
        map_tissues = map_tissues,
        map_species = map_species,
        metric_type = metric_type
    )
    getConversion(obj)
}

test_that("getGExMetrics returns a CoSIAn object", {
    skip_if_no_experimenthub()
    obj <- make_converted("DS_Gene")
    result <- getGExMetrics(obj)
    expect_s4_class(result, "CoSIAn")
})

test_that("DS_Gene metric slot is a non-empty data frame", {
    skip_if_no_experimenthub()
    obj <- make_converted("DS_Gene")
    result <- getGExMetrics(obj)
    expect_s3_class(result@metric, "data.frame")
    expect_gt(nrow(result@metric), 0)
})

test_that("DS_Gene metric has Ensembl_ID, Specificity, Diversity, Species columns", {
    skip_if_no_experimenthub()
    obj <- make_converted("DS_Gene")
    result <- getGExMetrics(obj)
    expect_true(all(c("Ensembl_ID", "Specificity", "Diversity", "Species") %in%
        colnames(result@metric)))
})

test_that("DS_Tissue metric has Anatomical_entity_name, Specificity, Diversity, Species columns", {
    skip_if_no_experimenthub()
    obj <- make_converted(
        "DS_Tissue",
        map_tissues = c("heart", "liver", "kidney")
    )
    result <- getGExMetrics(obj)
    expect_s3_class(result@metric, "data.frame")
    expect_true(all(c("Anatomical_entity_name", "Specificity", "Diversity", "Species") %in%
        colnames(result@metric)))
})

test_that("DS_Tissue errors with fewer than two genes", {
    skip_if_no_experimenthub()
    obj <- CoSIAn(
        gene_set = "ENSG00000008710",
        i_species = "h_sapiens",
        input_id = "Ensembl_id",
        o_species = c("h_sapiens", "m_musculus"),
        output_ids = "Ensembl_id",
        mapping_tool = "annotationDBI",
        ortholog_database = "HomoloGene",
        map_tissues = c("heart", "liver"),
        map_species = "m_musculus",
        metric_type = "DS_Tissue"
    )
    obj <- getConversion(obj)
    expect_error(getGExMetrics(obj), regexp = "more than one gene")
})

test_that("DS_Gene_all metric slot is a non-empty data frame", {
    skip_if_no_experimenthub()
    obj <- make_converted("DS_Gene_all")
    result <- getGExMetrics(obj)
    expect_s3_class(result@metric, "data.frame")
    expect_gt(nrow(result@metric), 0)
})

test_that("DS_Tissue_all metric slot is a non-empty data frame", {
    skip_if_no_experimenthub()
    obj <- make_converted(
        "DS_Tissue_all",
        map_tissues = c("heart", "liver", "kidney")
    )
    result <- getGExMetrics(obj)
    expect_s3_class(result@metric, "data.frame")
    expect_gt(nrow(result@metric), 0)
})

test_that("CV_Tissue metric slot is a non-empty data frame", {
    skip_if_no_experimenthub()
    obj <- make_converted("CV_Tissue")
    result <- getGExMetrics(obj)
    expect_s3_class(result@metric, "data.frame")
    expect_gt(nrow(result@metric), 0)
})

test_that("CV_Species metric slot is a non-empty data frame", {
    skip_if_no_experimenthub()
    obj <- make_converted("CV_Species")
    result <- getGExMetrics(obj)
    expect_s3_class(result@metric, "data.frame")
    expect_gt(nrow(result@metric), 0)
})

test_that("getGExMetrics errors on invalid metric_type", {
    skip_if_no_experimenthub()
    obj <- make_converted("DS_Gene")
    obj@metric_type <- "invalid_metric"
    expect_error(getGExMetrics(obj))
})
