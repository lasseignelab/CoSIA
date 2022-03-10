setMethod("setSlot<-", "CosiaAnnotate", function(object, slot, value) { # user's input of the function
  object@slot <- value
  object@slot
}
)