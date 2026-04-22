test_that("CoSIAn constructor returns a CoSIAn object", {
    obj <- make_cosian()
    expect_s4_class(obj, "CoSIAn")
})

test_that("constructor initializes input slots correctly", {
    obj <- make_cosian()
    expect_identical(obj@gene_set, .minimal_args$gene_set)
    expect_identical(obj@i_species, "h_sapiens")
    expect_identical(obj@input_id, "Ensembl_id")
    expect_identical(obj@o_species, c("h_sapiens", "m_musculus"))
    expect_identical(obj@output_ids, c("Ensembl_id", "Symbol"))
    expect_identical(obj@mapping_tool, "annotationDBI")
    expect_identical(obj@ortholog_database, "HomoloGene")
    expect_identical(obj@map_tissues, "heart")
    expect_identical(obj@map_species, "m_musculus")
    expect_identical(obj@metric_type, "DS_Gene")
})

test_that("constructor initializes output slots as empty data frames", {
    obj <- make_cosian()
    expect_s3_class(obj@converted_id, "data.frame")
    expect_s3_class(obj@gex, "data.frame")
    expect_s3_class(obj@metric, "data.frame")
})

test_that("mapping_tool defaults to annotationDBI", {
    obj <- CoSIAn(
        gene_set = .minimal_args$gene_set,
        i_species = "h_sapiens",
        input_id = "Ensembl_id",
        o_species = "h_sapiens",
        output_ids = "Ensembl_id",
        map_tissues = "heart",
        map_species = "h_sapiens",
        metric_type = "DS_Gene"
    )
    expect_identical(obj@mapping_tool, "annotationDBI")
})

test_that("ortholog_database defaults to HomoloGene", {
    obj <- CoSIAn(
        gene_set = .minimal_args$gene_set,
        i_species = "h_sapiens",
        input_id = "Ensembl_id",
        o_species = "h_sapiens",
        output_ids = "Ensembl_id",
        map_tissues = "heart",
        map_species = "h_sapiens",
        metric_type = "DS_Gene"
    )
    expect_identical(obj@ortholog_database, "HomoloGene")
})

test_that("constructor errors on missing i_species", {
    expect_error(
        make_cosian(i_species = character(0)),
        regexp = "i_species"
    )
})

test_that("constructor errors when i_species has length > 1", {
    expect_error(
        make_cosian(i_species = c("h_sapiens", "m_musculus")),
        regexp = "i_species"
    )
})

test_that("constructor errors on missing input_id", {
    expect_error(
        make_cosian(input_id = character(0)),
        regexp = "input_id"
    )
})

test_that("constructor errors when input_id has length > 1", {
    expect_error(
        make_cosian(input_id = c("Ensembl_id", "Symbol")),
        regexp = "input_id"
    )
})

test_that("constructor errors when mapping_tool has length > 1", {
    expect_error(
        make_cosian(mapping_tool = c("annotationDBI", "biomaRt")),
        regexp = "mapping_tool"
    )
})

test_that("constructor errors when ortholog_database has length > 1", {
    expect_error(
        make_cosian(ortholog_database = c("HomoloGene", "NCBIOrtho")),
        regexp = "ortholog_database"
    )
})

test_that("constructor errors on empty gene_set", {
    expect_error(
        make_cosian(gene_set = character(0)),
        regexp = "gene_set"
    )
})
