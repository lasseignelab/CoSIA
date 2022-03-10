#All the package Generics

#Show Generics for CosiaAnnotate Class
setGeneric("showCosiaAnnotate", function(x) standardGeneric("showCosiaAnnotate"))
setGeneric("showCosiaExpressSpecies", function(x) standardGeneric("showCosiaExpressSpecies"))
setGeneric("showCosiaExpressTissue", function(x) standardGeneric("showCosiaExpressTissue"))

# Accessor Generics for CosiaAnnotate Class
setGeneric("getSlot", function(object) standardGeneric("getSlot"))

setGeneric("setSlot<-", function(object, value) standardGeneric("setSlot<-"))


# Method Generics
setGeneric("getConversion", function(object) standardGeneric("getConversion"))

setGeneric("getTissueExpression", function(object) standardGeneric("getTissueExpression"))

setGeneric("getSpeciesExpression", function(object) standardGeneric("getSpeciesExpression"))