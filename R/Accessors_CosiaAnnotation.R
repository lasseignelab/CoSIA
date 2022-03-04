setMethod("get_input", "CosiaAnnotate", function(x) x@input)
setMethod("set_input<-", "CosiaAnnotate", function(x, value) {
  x@input <- value
  x
})
