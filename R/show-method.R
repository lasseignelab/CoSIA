setMethod("showCosiaAnnotate", "CosiaAnnotate", function(object) {
    cat(is(object)[[1]], "\n", "  input: ", object@input, "\n", "  input_species:  ", object@input_species, "\n", "  input_id: ", object@input_id,
        "\n", "  output_ids:  ", object@output_ids, "\n", "  output_species: ", object@output_species, "\n", "  tool:  ", object@tool, "\n",
        "  ortholog_database: ", object@ortholog_database, "\n", sep = "")
})

setMethod("showCosiaExpressSpecies", "CosiaExpressSpecies", function(object) {
    cat(is(object)[[1]], "\n", "  list_of_ensembl_ids: ", object@list_of_ensembl_ids, "\n", "  list_of_respective_species:  ", object@list_of_respective_species,
        "\n", "  single_tissue: ", object@single_tissue, "\n", "  pathToData:  ", object@pathToData, "\n", "  dataframe: ", object@dataframe,
        "\n", sep = "")
})

setMethod("showCosiaExpressTissue", "CosiaExpressTissue", function(object) {
    cat(is(object)[[1]], "\n", "  single_gene: ", object@single_gene, "\n", "  gene_species:  ", object@gene_species, "\n", "  tissues: ", object@tissues,
        "\n", "  pathToData:  ", object@pathToData, "\n", "  dataframe: ", object@dataframe, "\n", sep = "")
})
