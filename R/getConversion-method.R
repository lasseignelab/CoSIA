
setMethod("getConversion", signature(object = "CosiaAnnotate"), function(object) { # user's input of the function
  if (object@input_species=="homo_sapiens"){ #code follows this path if the user chooses homo_sapiens as their input species
    human_data<-homo_sapiens(object@input_id,object@input,object@output_ids,object@output_species, object@tool, object@ortholog_database)
    return(human_data)
  }
  if (object@input_species=="d_melanogaster"){ #code follows this path if the user chooses d_melanogaster as their input species
    fly_data<-d_melanogaster(object@input_id,object@input,object@output_ids,object@output_species, object@tool, object@ortholog_database)
    return(fly_data)
  }
  if(object@input_species=="mus_musculus"){ #code follows this path if the user chooses mus_musculus as their input species
    mouse_data<-mus_musculus(object@input_id,object@input,object@output_ids,object@output_species, object@tool, object@ortholog_database)
    return(mouse_data)
  }
  if(object@input_species=="danio_rerio"){  #code follows this path if the user chooses danio_rerio as their input species
    zebrafish_data<-danio_rerio(object@input_id,object@input,object@output_ids,object@output_species, object@tool, object@ortholog_database)
    return(zebrafish_data)
  }
  if(object@input_species=="c_elegans"){  #code follows this path if the user chooses c_elegans as their input species
    celegan_data<-c_elegans(object@input_id,object@input,object@output_ids,object@output_species, object@tool, object@ortholog_database)
    return(celegan_data)
  }
  if(object@input_species=="r_norvegicus"){  #code follows this path if the user chooses r_norvegicus as their input species
    celegan_data<-r_norvegicus(object@input_id,object@input,object@output_ids,object@output_species, object@tool, object@ortholog_database)
    return(celegan_data)
  }
  else{ #code follows this path if the inout species is not correctly written or is not one of the available species.
    stop("Error. Invalid input species. Make sure it matches the proper format")
  }
})
