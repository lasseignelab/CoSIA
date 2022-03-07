#All the package Generics

#Show Generics for CosiaAnnotate Class
setGeneric("showCosiaAnnotate", function(x) standardGeneric("showCosiaAnnotate"))

# Accessor Generics for CosiaAnnotate Class

setGeneric("get_input", function(x) standardGeneric("get_input"))
setGeneric("set_input<-", function(x, value) standardGeneric("set_input<-"))

setGeneric("get_input_species", function(x) standardGeneric("get_input_species"))
setGeneric("set_input_species<-", function(x, value) standardGeneric("set_input_species<-"))

setGeneric("get_input_id", function(x) standardGeneric("get_input_id"))
setGeneric("set_input_id<-", function(x, value) standardGeneric("set_input_id<-"))

setGeneric("get_output_ids", function(x) standardGeneric("get_output_ids"))
setGeneric("set_output_ids<-", function(x, value) standardGeneric("set_output_ids<-"))

setGeneric("get_output_species", function(x) standardGeneric("get_output_species"))
setGeneric("set_output_species<-", function(x, value) standardGeneric("set_output_species<-"))

setGeneric("get_tool", function(x) standardGeneric("get_tool"))
setGeneric("set_tool<-", function(x, value) standardGeneric("set_tool<-"))

setGeneric("get_ortholog_database", function(x) standardGeneric("get_ortholog_database"))
setGeneric("set_ortholog_database", function(x, value) standardGeneric("set_ortholog_database"))

# Method Generics
setGeneric("getConversion", function(object) standardGeneric("getConversion"))

setGeneric("getTissueExpression", function(object) standardGeneric("getTissueExpression"))

setGeneric("getSpeciesExpression", function(object) standardGeneric("getSpeciesExpression"))