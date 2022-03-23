#Getter and Setter for CosiaAnnotate
setMethod("input", "CosiaAnnotate", function(x) x@input)
setMethod("input<-", "CosiaAnnotate", function(x, value) {
  x@input <- value
  x
})

setMethod("input_species", "CosiaAnnotate", function(x) x@input_species)
setMethod("input_species<-", "CosiaAnnotate", function(x, value) {
  x@input_species <- value
  x
})

setMethod("input_id", "CosiaAnnotate", function(x) x@input_id)
setMethod("input_id<-", "CosiaAnnotate", function(x, value) {
  x@input_id <- value
  x
})

setMethod("output_species", "CosiaAnnotate", function(x) x@output_species)
setMethod("output_species<-", "CosiaAnnotate", function(x, value) {
  x@output_species <- value
  x
})

setMethod("output_ids", "CosiaAnnotate", function(x) x@output_ids)
setMethod("output_ids<-", "CosiaAnnotate", function(x, value) {
  x@output_ids <- value
  x
})

setMethod("tool", "CosiaAnnotate", function(x) x@tool)
setMethod("tool<-", "CosiaAnnotate", function(x, value) {
  x@tool <- value
  x
})

setMethod("ortholog_database", "CosiaAnnotate", function(x) x@ortholog_database)
setMethod("ortholog_database<-", "CosiaAnnotate", function(x, value) {
  x@ortholog_database <- value
  x
})


#Getter and Setter for CosiaAnnotate
setMethod("input", "CosiaAnnotate", function(x) x@input)
setMethod("input<-", "CosiaAnnotate", function(x, value) {
  x@input <- value
  x
})

setMethod("input_species", "CosiaAnnotate", function(x) x@input_species)
setMethod("input_species<-", "CosiaAnnotate", function(x, value) {
  x@input_species <- value
  x
})

setMethod("input_id", "CosiaAnnotate", function(x) x@input_id)
setMethod("input_id<-", "CosiaAnnotate", function(x, value) {
  x@input_id <- value
  x
})

setMethod("output_species", "CosiaAnnotate", function(x) x@output_species)
setMethod("output_species<-", "CosiaAnnotate", function(x, value) {
  x@output_species <- value
  x
})

setMethod("output_ids", "CosiaAnnotate", function(x) x@output_ids)
setMethod("output_ids<-", "CosiaAnnotate", function(x, value) {
  x@output_ids <- value
  x
})

setMethod("tool", "CosiaAnnotate", function(x) x@tool)
setMethod("tool<-", "CosiaAnnotate", function(x, value) {
  x@tool <- value
  x
})

setMethod("ortholog_database", "CosiaAnnotate", function(x) x@ortholog_database)
setMethod("ortholog_database<-", "CosiaAnnotate", function(x, value) {
  x@ortholog_database <- value
  x
})



