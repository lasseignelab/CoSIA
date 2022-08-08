# All the package Generics

# Show Generics
setGeneric("showCosiaAnnotate", function(object) standardGeneric("showCosiaAnnotate"))
setGeneric("showCosiaExpressSpecies", function(object) standardGeneric("showCosiaExpressSpecies"))
setGeneric("showCosiaExpressTissue", function(object) standardGeneric("showCosiaExpressTissue"))

# Accessor Generics for CosiaAnnotate Class

setGeneric("input", function(x) standardGeneric("input"))
setGeneric("input<-", function(x, value) standardGeneric("input<-"))

setGeneric("input_species", function(x) standardGeneric("input_species"))
setGeneric("input_species<-", function(x, value) standardGeneric("input_species<-"))

setGeneric("input_id", function(x) standardGeneric("input_id"))
setGeneric("input_id<-", function(x, value) standardGeneric("input_id<-"))

setGeneric("output_species", function(x) standardGeneric("output_species"))
setGeneric("output_species<-", function(x, value) standardGeneric("output_species<-"))

setGeneric("output_ids", function(x) standardGeneric("output_ids"))
setGeneric("output_ids<-", function(x, value) standardGeneric("output_ids<-"))

setGeneric("tool", function(x) standardGeneric("tool"))
setGeneric("tool<-", function(x, value) standardGeneric("tool<-"))

setGeneric("ortholog_database", function(x) standardGeneric("ortholog_database"))
setGeneric("ortholog_database<-", function(x, value) standardGeneric("ortholog_database<-"))


# Accessor Generics for CosiaExpressSpecies Class
setGeneric("list_of_ensembl_ids", function(x) standardGeneric("list_of_ensembl_ids"))
setGeneric("list_of_ensembl_ids<-", function(x, value) standardGeneric("list_of_ensembl_ids<-"))

setMethod("showCosiaExpressSpecies", "CosiaExpressSpecies", function(object) {
    cat(is(object)[[1]], "\n", "  list_of_ensembl_ids: ", object@list_of_ensembl_ids, "\n", "  list_of_respective_species:  ", object@list_of_respective_species,
        "\n", "  single_tissue: ", object@single_tissue, "\n", "  pathToData:  ", object@pathToData, "\n", "  dataframe: ", object@dataframe,
        "\n", sep = "")
})

setMethod("showCosiaExpressTissue", "CosiaExpressTissue", function(object) {
    cat(is(object)[[1]], "\n", "  single_gene: ", object@single_gene, "\n", "  gene_species:  ", object@gene_species, "\n", "  tissues: ", object@tissues,
        "\n", "  pathToData:  ", object@pathToData, "\n", "  dataframe: ", object@dataframe, "\n", sep = "")
})







# Method Generics
setGeneric("getConversion", function(object) standardGeneric("getConversion"))

setGeneric("getTissueExpression", function(object) standardGeneric("getTissueExpression"))

setGeneric("getSpeciesExpression", function(object) standardGeneric("getSpeciesExpression"))
