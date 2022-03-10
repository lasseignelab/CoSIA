setMethod("setExpressSlot<-", "CosiaExpress", function(object, slot, value) {
  object@slot <- value
  object
})