# Shared fixtures for CoSIA tests

.minimal_args <- list(
    gene_set = c("ENSG00000008710", "ENSG00000118762", "ENSG00000152217"),
    i_species = "h_sapiens",
    input_id = "Ensembl_id",
    o_species = c("h_sapiens", "m_musculus"),
    output_ids = c("Ensembl_id", "Symbol"),
    mapping_tool = "annotationDBI",
    ortholog_database = "HomoloGene",
    map_tissues = "heart",
    map_species = "m_musculus",
    metric_type = "DS_Gene"
)

# Build a CoSIAn with any slot overridden
make_cosian <- function(...) {
    args <- utils::modifyList(.minimal_args, list(...))
    do.call(CoSIAn, args)
}

# Skip if ExperimentHub is not reachable
skip_if_no_experimenthub <- function() {
    skip_if_offline()
    tryCatch(
        ExperimentHub::ExperimentHub(),
        error = function(e) skip(paste("ExperimentHub unavailable:", conditionMessage(e)))
    )
}
