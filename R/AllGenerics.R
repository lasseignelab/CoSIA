#All the package Generics

#Show Generics for CosiaAnnotate Class
setGeneric("showCosiaAnnotate", function(x) standardGeneric("showCosiaAnnotate"))
setGeneric("showCosiaExpressTissue", function(x) standardGeneric("showCosiaExpressSpecies"))
setGeneric("showCosiaExpressTissue", function(x) standardGeneric("showCosiaExpressTissue"))

# Accessor Generics for CosiaAnnotate Class
setGeneric("get", function(x,y) standardGeneric("get_input"))
setGeneric("set<-", function(x,y, value) standardGeneric("set_input<-"))


# Method Generics
setGeneric("getConversion", function(object) standardGeneric("getConversion"))

setGeneric("getTissueExpression", function(object) standardGeneric("getTissueExpression"))

setGeneric("getSpeciesExpression", function(object) standardGeneric("getSpeciesExpression"))