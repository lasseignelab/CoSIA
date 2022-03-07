setMethod("showCosiaExpressSpecies", "CosiaExpressSpecies", function(object) {
  cat(is(object)[[1]], "\n",
      "  list_of_ensembl_ids: ", object@list_of_ensembl_ids, "\n",
      "  list_of_respective_species:  ", object@list_of_respective_species, "\n",
      "  single_tissue: ", object@single_tissue, "\n",
      "  pathToData:  ", object@pathToData, "\n",
      "  dataframe: ", object@dataframe, "\n",
      sep = ""
  )
})