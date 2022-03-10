#All the package Generics

#Show Generics for CosiaAnnotate Class
setGeneric("showCosiaAnnotate", function(object) standardGeneric("showCosiaAnnotate"))
setGeneric("showCosiaExpressSpecies", function(object) standardGeneric("showCosiaExpressSpecies"))
setGeneric("showCosiaExpressTissue", function(object) standardGeneric("showCosiaExpressTissue"))

# Accessor Generics for CosiaAnnotate Class
setGeneric("getAnnotateSlot", function(object, slot) standardGeneric("getAnnotateSlot"))
setGeneric("setAnnotateSlot<-", function(object, slot, value) standardGeneric("setAnnotateSlot<-"))

# Accessor Generics for CosiaExpress Class
setGeneric("getExpressSlot", function(object, slot) standardGeneric("getExpressSlot"))
setGeneric("setExpressSlot<-", function(object, slot, value) standardGeneric("setExpressSlot<-"))


# Method Generics
setGeneric("getConversion", function(object) standardGeneric("getConversion"))

setGeneric("getTissueExpression", function(object) standardGeneric("getTissueExpression"))

setGeneric("getSpeciesExpression", function(object) standardGeneric("getSpeciesExpression"))