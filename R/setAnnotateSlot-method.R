setMethod("setAnnotateSlot<-", "CosiaAnnotate", function(object, slot, value) {
  object@slot <- value
  object
})