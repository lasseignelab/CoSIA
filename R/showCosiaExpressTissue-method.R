setMethod("showCosiaExpressTissue", "CosiaExpressTissue", function(object) {
  cat(is(object)[[1]], "\n",
      "  single_gene: ", object@single_gene, "\n",
      "  gene_species:  ", object@gene_species, "\n",
      "  tissues: ", object@tissues, "\n",
      "  pathToData:  ", object@pathToData, "\n",
      "  dataframe: ", object@dataframe, "\n",
      sep = ""
  )
})