setMethod("showCosiaAnnotate", "CosiaAnnotate", function(object) {
    cat(is(object)[[1]], "\n",
        "  input: ", object@input, "\n",
        "  input_species:  ", object@input_species, "\n",
        "  input_id: ", object@input_id, "\n",
        "  output_ids:  ", object@output_ids, "\n",
        "  output_species: ", object@output_species, "\n",
        "  tool:  ", object@tool, "\n",
        "  ortholog_database: ", object@ortholog_database, "\n",
        sep = ""
    )
  })