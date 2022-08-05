#' Converting Gene Identifiers between the same species and different species
#'
#' @param object CosiaAnnotate. 
#'
#' @return A dataframe with input identifiers and matching gene identifier conversion of the output species 
#' @export
#'
#' @examples
#' getConversion(mouse_data)
#' getConversion(rat_data)

setMethod("getConversion", signature(object = "CosiaAnnotate"), function(object) { # user's input of the function
  #set the parts of an object into variables
  input_species<-object@input_species
  input_id<-object@input_id
  input<-object@input
  output_ids<-object@output_ids
  output_species<-object@output_species
  tool<-object@tool
  ortholog_database<-object@ortholog_database
  if (input_species=="homo_sapiens"){ #code follows this path if the user chooses homo_sapiens as their input species
    human_data<-homo_sapiens(input_id,input,output_ids,output_species, tool, ortholog_database)
    #take this out if you are unable to fix it and make sure it works before you publish
    return(human_data)
  }
  if (input_species=="d_melanogaster"){ #code follows this path if the user chooses d_melanogaster as their input species
    fly_data<-d_melanogaster(input_id,input,output_ids,output_species, tool, ortholog_database)
    return(fly_data)
  }
  if(input_species=="mus_musculus"){ #code follows this path if the user chooses mus_musculus as their input species
    mouse_data<-mus_musculus(input_id,input,output_ids,output_species, tool, ortholog_database)
    return(mouse_data)
  }
  if(input_species=="danio_rerio"){  #code follows this path if the user chooses danio_rerio as their input species
    zebrafish_data<-danio_rerio(input_id,input,output_ids,output_species, tool, ortholog_database)
    return(zebrafish_data)
  }
  if(input_species=="c_elegans"){  #code follows this path if the user chooses c_elegans as their input species
    celegan_data<-c_elegans(input_id,input,output_ids,output_species, tool, ortholog_database)
    return(celegan_data)
  }
  if(input_species=="r_norvegicus"){  #code follows this path if the user chooses r_norvegicus as their input species
    rat_data<-r_norvegicus(input_id,input,output_ids,output_species, tool, ortholog_database)
    return(rat_data)
  }
  else{ #code follows this path if the inout species is not correctly written or is not one of the available species.
    stop("Error. Invalid input species. Make sure it matches the proper format")
  }
})
