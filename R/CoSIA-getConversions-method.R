#' getConversions Generic
#'
#' @param object
#'
#' @return
#' @export
#'
#' @examples

setGeneric("getConversions", function(object) standardGeneric("getConversions"))

#' getConversions Method
#'
#' @param object CoSIAn. 
#'
#' @return 
#' @export
#'
#' @examples

setMethod("getConversions", signature(object = "CoSIAn"), function(object) { # user's input of the function
  #Set each part of the object that this method uses into their own variable that will be used inside the code
    input_species<-object@i_species
    input_id<-object@input_id
    input<-object@gene_set
    output_ids<-object@o_ids
    output_species<-object@o_species
    tool<-object@mapping_tool
    ortholog_database<-object@ortholog_database
    Filter_I_Species <- switch(input_species,
                               h_sapiens = {
                                 species_data<-h_sapiens(input_id,input,output_ids,output_species, tool, ortholog_database)
                               },
                               m_musculus = {
                                 species_data<-m_musculus(input_id,input,output_ids,output_species, tool, ortholog_database)
                               },
                               r_norvegicus = {
                                 species_data<-r_norvegicus(input_id,input,output_ids,output_species, tool, ortholog_database)
                               },
                               d_rerio = {
                                 species_data<-d_rerio(input_id,input,output_ids,output_species, tool, ortholog_database)
                               },
                               c_elegans = {
                                 species_data<-c_elegans(input_id,input,output_ids,output_species, tool, ortholog_database)
                               },
                               d_melanogaster = {
                                 species_data<-d_melanogaster(input_id,input,output_ids,output_species, tool, ortholog_database)
                               },
                               stop("Error: Invalid i_species in CoSIAn Object. Make sure the species in the i_species slot is an avalible model 
                                    organism and is in the correct format.")
                               
    )
    
    slot(object, converted_id, check = TRUE) <- data.frame(species_data)
    converted_id<-slot(object, converted_id)
    return(converted_id)
})

### getConversions Functions

## Species Functions

# Caenorhabditis elegans - Roundworm

c_elegans<-function(input_id,input_dataset,output_ids,output_species, tool, ortholog_database) { #input data funneling in from the cross species conversion function
  if (tool=="biomaRt"){ # if the user has choose the tool biomart this it the path the codes follows
    Filter_BO_CE <- switch(output_species,
                           c_elegans = {
                             output_data<-CoSIA::biomaRt(input_id,input_dataset, output_ids, input_species="c_elegans", output_species, 6239, "celegans_gene_ensembl", "celegans_gene_ensembl", ortholog_database)
                             },
                           d_melanogaster = {
                             output_data<-CoSIA::biomaRt(input_id,input_dataset, output_ids, input_species="c_elegans", output_species, 7227, "celegans_gene_ensembl", "dmelanogaster_gene_ensembl", ortholog_database)
                             },
                           m_musculus = {
                             output_data<-CoSIA::biomaRt(input_id,input_dataset, output_ids,input_species="c_elegans", output_species, 10090, "celegans_gene_ensembl", "mmusculus_gene_ensembl", ortholog_database)
                             },
                           d_rerio = {
                             output_data<-CoSIA::biomaRt(input_id,input_dataset, output_ids, input_species="c_elegans",output_species, 7955, "celegans_gene_ensembl", "drerio_gene_ensembl", ortholog_database)
                             },
                           h_sapien = {
                             output_data<-CoSIA::biomaRt(input_id,input_dataset, output_ids,input_species="c_elegans", output_species, 9606, "celegans_gene_ensembl", "hsapiens_gene_ensembl", ortholog_database)
                             },
                           r_norvegicus ={
                             output_data<-CoSIA::biomaRt(input_id,input_dataset, output_ids,input_species="c_elegans", output_species, 10116, "celegans_gene_ensembl", "rnorvegicus_gene_ensembl", ortholog_database)
                             },
                             stop("Error: Invalid o_species in CoSIAn Object. Make sure the species in the o_species slot is an avalible model 
                                    organism and is in the correct format.")
                           )
  }
  if(tool=="annotationDBI"){# code follows this path if the user chooses annotationDBI as their tool of choice
    Filter_AO_CE <- switch(output_species,
                           h_sapiens = {
                             output_data<-CoSIA::annotateDBI(input_id,input_dataset,output_ids,input_species="c_elegans",output_species,9606, org.Ce.eg.db::org.Ce.eg.db , org.Hs.eg.db::org.Hs.eg.db, ortholog_database)
                           },
                           m_musculus = {
                             output_data<-CoSIA::annotateDBI(input_id,input_dataset,output_ids,input_species="c_elegans",output_species,10090, org.Ce.eg.db::org.Ce.eg.db , org.Mm.eg.db::org.Mm.eg.db, ortholog_database)
                           },
                           d_melanogaster = {
                             output_data<-CoSIA::annotateDBI(input_id,input_dataset,output_ids,input_species="c_elegans",output_species,7227, org.Ce.eg.db::org.Ce.eg.db , org.Dm.eg.db::org.Dm.eg.db, ortholog_database)
                           },
                           danio_rerio = {
                             output_data<-CoSIA::annotateDBI(input_id,input_dataset,output_ids,input_species="c_elegans",output_species,7955, org.Ce.eg.db::org.Ce.eg.db , org.Dr.eg.db::org.Dr.eg.db, ortholog_database)
                           },
                           c_elegans = {
                             output_data<-CoSIA::annotateDBI(input_id,input_dataset,output_ids,input_species="c_elegans",output_species,6239, org.Ce.eg.db::org.Ce.eg.db , org.Ce.eg.db::org.Ce.eg.db, ortholog_database)
                           },
                           r_norvegicus ={
                             output_data<-CoSIA::annotateDBI(input_id,input_dataset,output_ids,input_species="c_elegans",output_species,10116, org.Ce.eg.db::org.Ce.eg.db , org.Rn.eg.db::org.Rn.eg.db,ortholog_database)
                           },
                           stop("Error: Invalid o_species in CoSIAn Object. Make sure the species in the o_species slot is an avalible model 
                                    organism and is in the correct format.")
                           )
  }
  else{ # code follows this path if the tool that was inputted into the system does not match.
    stop("Error: Invalid tool in CoSIAn Object. Make sure the tool is either annotationDBI or biomaRt.")
  }
}

#Drosophila Melanogaster - Fruit Fly

d_melanogaster<-function(input_id,input_dataset,output_ids,output_species, tool, ortholog_database) { #input data funneling in from the cross species conversion function
  if (tool=="biomaRt"){ # if the user has choose the tool biomart this it the path the codes follows
    Filter_BO_DM <- switch(output_species,
                           c_elegans = {
                             output_data<-CoSIA::biomaRt(input_id,input_dataset, output_ids, input_species="d_melanogaster",output_species, 6239, "dmelanogaster_gene_ensembl", "celegans_gene_ensembl", ortholog_database)
                           },
                           d_melanogaster = {
                             output_data<-CoSIA::biomaRt(input_id,input_dataset, output_ids,input_species="d_melanogaster", output_species, 7227, "dmelanogaster_gene_ensembl", "dmelanogaster_gene_ensembl", ortholog_database)
                           },
                           m_musculus = {
                             output_data<-CoSIA::biomaRt(input_id,input_dataset, output_ids,input_species="d_melanogaster", output_species, 10090, "dmelanogaster_gene_ensembl", "mmusculus_gene_ensembl", ortholog_database)
                           },
                           d_rerio = {
                             output_data<-CoSIA::biomaRt(input_id,input_dataset, output_ids,input_species="d_melanogaster", output_species, 7955, "dmelanogaster_gene_ensembl", "drerio_gene_ensembl", ortholog_database)
                           },
                           h_sapien = {
                             output_data<-CoSIA::biomaRt(input_id,input_dataset, output_ids,input_species="d_melanogaster", output_species, 9606, "dmelanogaster_gene_ensembl", "hsapiens_gene_ensembl", ortholog_database)
                           },
                           r_norvegicus ={
                             output_data<-CoSIA::biomaRt(input_id,input_dataset, output_ids,input_species="d_melanogaster", output_species, 10116, "dmelanogaster_gene_ensembl", "rnorvegicus_gene_ensembl", ortholog_database)
                           },
                           stop("Error: Invalid o_species in CoSIAn Object. Make sure the species in the o_species slot is an avalible model 
                                    organism and is in the correct format.")
    )
  }
  if(tool=="annotationDBI"){# code follows this path if the user chooses annotationDBI as their tool of choice
    Filter_AO_DM <- switch(output_species,
                           h_sapiens = {
                             output_data<-CoSIA::annotateDBI(input_id,input_dataset,output_ids,input_species="d_melanogaster",output_species,9606, org.Dm.eg.db::org.Dm.eg.db , org.Hs.eg.db::org.Hs.eg.db, ortholog_database)
                           },
                           m_musculus = {
                             output_data<-CoSIA::annotateDBI(input_id,input_dataset,output_ids,input_species="d_melanogaster",output_species,10090, org.Dm.eg.db::org.Dm.eg.db , org.Mm.eg.db::org.Mm.eg.db, ortholog_database)
                           },
                           d_melanogaster = {
                             output_data<-CoSIA::annotateDBI(input_id,input_dataset,output_ids,input_species="d_melanogaster",output_species,7227, org.Dm.eg.db::org.Dm.eg.db , org.Dm.eg.db::org.Dm.eg.db, ortholog_database)
                           },
                           danio_rerio = {
                             output_data<-CoSIA::annotateDBI(input_id,input_dataset,output_ids,input_species="d_melanogaster",output_species,7955, org.Dm.eg.db::org.Dm.eg.db , org.Dr.eg.db::org.Dr.eg.db, ortholog_database)
                           },
                           c_elegans = {
                             output_data<-CoSIA::annotateDBI(input_id,input_dataset,output_ids,input_species="d_melanogaster",output_species,6239, org.Dm.eg.db::org.Dm.eg.db , org.Ce.eg.db::org.Ce.eg.db, ortholog_database)
                           },
                           r_norvegicus ={
                             output_data<-CoSIA::annotateDBI(input_id,input_dataset,output_ids,input_species="d_melanogaster",output_species,10116, org.Dm.eg.db::org.Dm.eg.db , org.Rn.eg.db::org.Rn.eg.db, ortholog_database)
                           },
                           stop("Error: Invalid o_species in CoSIAn Object. Make sure the species in the o_species slot is an avalible model 
                                    organism and is in the correct format.")
    )
  }
  else{ # code follows this path if the tool that was inputted into the system does not match.
    stop("Error: Invalid tool in CoSIAn Object. Make sure the tool is either annotationDBI or biomaRt.")
  }
}

#Danio_Rerio - Zebrafish

d_rerio<-function(input_id,input_dataset,output_ids,output_species, tool, ortholog_database) { #input data funneling in from the cross species conversion function
  if (tool=="biomaRt"){ # if the user has choose the tool biomart this it the path the codes follows
    Filter_BO_DR <- switch(output_species,
                           c_elegans = {
                             output_data <- CoSIA::biomaRt(input_id, input_dataset, output_ids, input_species = "d_rerio", output_species, 6239, "drerio_gene_ensembl",
                                                        "celegans_gene_ensembl", ortholog_database)                           },
                           d_melanogaster = {
                             output_data <- CoSIA::biomaRt(input_id, input_dataset, output_ids, input_species = "d_rerio", output_species, 7227, "drerio_gene_ensembl",
                                                        "dmelanogaster_gene_ensembl", ortholog_database)                           },
                           m_musculus = {
                             output_data <- CoSIA::biomaRt(input_id, input_dataset, output_ids, input_species = "d_rerio", output_species, 10090, "drerio_gene_ensembl",
                                                        "mmusculus_gene_ensembl", ortholog_database)                           },
                           d_rerio = {
                             output_data <- CoSIA::biomaRt(input_id, input_dataset, output_ids, input_species = "d_rerio", output_species, 7955, "drerio_gene_ensembl","drerio_gene_ensembl", ortholog_database)                           },
                           h_sapien = {
                             output_data <- CoSIA::biomaRt(input_id, input_dataset, output_ids, input_species = "d_rerio", output_species, 9606, "drerio_gene_ensembl",
                                                        "hsapiens_gene_ensembl", ortholog_database)                           },
                           r_norvegicus ={
                             output_data <- CoSIA::biomaRt(input_id, input_dataset, output_ids, input_species = "d_rerio", output_species, 10116, "drerio_gene_ensembl",
                                                        "rnorvegicus_gene_ensembl", ortholog_database)                           },
                           stop("Error: Invalid o_species in CoSIAn Object. Make sure the species in the o_species slot is an avalible model 
                                    organism and is in the correct format.")
    )
  }
  if(tool=="annotationDBI"){# code follows this path if the user chooses annotationDBI as their tool of choice
    Filter_AO_DR <- switch(output_species,
                           h_sapiens = {
                             output_data <- CoSIA::annotateDBI(input_id, input_dataset, output_ids, input_species = "d_rerio", output_species, 9606, org.Dr.eg.db::org.Dr.eg.db,
                                                               org.Hs.eg.db::org.Hs.eg.db, ortholog_database)                           },
                           m_musculus = {
                             output_data <- CoSIA::annotateDBI(input_id, input_dataset, output_ids, input_species = "d_rerio", output_species, 10090, org.Dr.eg.db::org.Dr.eg.db,
                                                               org.Mm.eg.db::org.Mm.eg.db, ortholog_database)                           },
                           d_melanogaster = {
                             output_data <- CoSIA::annotateDBI(input_id, input_dataset, output_ids, input_species = "d_rerio", output_species, 7227, org.Dr.eg.db::org.Dr.eg.db,
                                                               org.Dm.eg.db::org.Dm.eg.db, ortholog_database)                           },
                           danio_rerio = {
                             output_data <- CoSIA::annotateDBI(input_id, input_dataset, output_ids, input_species = "d_rerio", output_species, 7955, org.Dr.eg.db::org.Dr.eg.db,
                                                               org.Dr.eg.db::org.Dr.eg.db, ortholog_database)                           },
                           c_elegans = {
                             output_data <- CoSIA::annotateDBI(input_id, input_dataset, output_ids, input_species = "d_rerio", output_species, 6239, org.Dr.eg.db::org.Dr.eg.db,
                                                               org.Ce.eg.db::org.Ce.eg.db, ortholog_database)                           },
                           r_norvegicus ={
                             output_data <- CoSIA::annotateDBI(input_id, input_dataset, output_ids, input_species = "d_rerio", output_species, 6239, org.Dr.eg.db::org.Dr.eg.db,
                                                               org.Ce.eg.db::org.Ce.eg.db, ortholog_database)                           },
                           stop("Error: Invalid o_species in CoSIAn Object. Make sure the species in the o_species slot is an avalible model 
                                    organism and is in the correct format.")
    )
  }
  else{ # code follows this path if the tool that was inputted into the system does not match.
    stop("Error: Invalid tool in CoSIAn Object. Make sure the tool is either annotationDBI or biomaRt.")
  }
}

#Homo sapiens - Humans

h_sapiens<-function(input_id,input_dataset,output_ids,output_species, tool, ortholog_database) { #input data funneling in from the cross species conversion function
  if (tool=="biomaRt"){ # if the user has choose the tool biomart this it the path the codes follows
    Filter_BO_HS <- switch(output_species,
                           c_elegans = {
                             output_data <- CoSIA::biomaRt(input_id, input_dataset, output_ids, input_species = "homo_sapiens", output_species, 6239, "hsapiens_gene_ensembl",
                                                        "drerio_gene_ensembl", ortholog_database)                           },
                           d_melanogaster = {
                             output_data <- CoSIA::biomaRt(input_id, input_dataset, output_ids, input_species = "homo_sapiens", output_species, 7227, "hsapiens_gene_ensembl",
                                                        "dmelanogaster_gene_ensembl", ortholog_database)                          },
                           m_musculus = {
                             output_data <- CoSIA::biomaRt(input_id, input_dataset, output_ids, input_species = "homo_sapiens", output_species, 10090, "hsapiens_gene_ensembl",
                                                        "mmusculus_gene_ensembl", ortholog_database)                           },
                           d_rerio = {
                             output_data <- CoSIA::biomaRt(input_id, input_dataset, output_ids, input_species = "homo_sapiens", output_species, 7955, "hsapiens_gene_ensembl",
                                                        "drerio_gene_ensembl", ortholog_database)                               },
                           h_sapien = {
                             output_data <- CoSIA::biomaRt(input_id, input_dataset, output_ids, input_species = "homo_sapiens", output_species, 9606, "hsapiens_gene_ensembl",
                                                        "hsapiens_gene_ensembl", ortholog_database)                           },
                           r_norvegicus ={
                             output_data <- CoSIA::biomaRt(input_id, input_dataset, output_ids, input_species = "homo_sapiens", output_species, 10116, "hsapiens_gene_ensembl",
                                                        "rnorvegicus_gene_ensembl", ortholog_database)                           },
                           stop("Error: Invalid o_species in CoSIAn Object. Make sure the species in the o_species slot is an avalible model 
                                    organism and is in the correct format.")
    )
  }
  if(tool=="annotationDBI"){# code follows this path if the user chooses annotationDBI as their tool of choice
    Filter_AO_HS <- switch(output_species,
                           h_sapiens = {
                             output_data <- CoSIA::annotateDBI(input_id, input_dataset, output_ids, input_species = "homo_sapiens", output_species, 9606, org.Hs.eg.db::org.Hs.eg.db,
                                                               org.Hs.eg.db::org.Hs.eg.db, ortholog_database)                          },
                           m_musculus = {
                             output_data <- CoSIA::annotateDBI(input_id, input_dataset, output_ids, input_species = "homo_sapiens", output_species, 10090, org.Hs.eg.db::org.Hs.eg.db,
                                                               org.Mm.eg.db::org.Mm.eg.db, ortholog_database)                           },
                           d_melanogaster = {
                             output_data <- CoSIA::annotateDBI(input_id, input_dataset, output_ids, input_species = "homo_sapiens", output_species, 7227, org.Hs.eg.db::org.Hs.eg.db,
                                                               org.Dm.eg.db::org.Dm.eg.db, ortholog_database)                           },
                           danio_rerio = {
                             output_data <- CoSIA::annotateDBI(input_id, input_dataset, output_ids, input_species = "homo_sapiens", output_species, 7955, org.Hs.eg.db::org.Hs.eg.db,
                                                               org.Dr.eg.db::org.Dr.eg.db, ortholog_database)                           },
                           c_elegans = {
                             output_data <- CoSIA::annotateDBI(input_id, input_dataset, output_ids, input_species = "homo_sapiens", output_species, 6239, org.Hs.eg.db::org.Hs.eg.db,
                                                               org.Ce.eg.db::org.Ce.eg.db, ortholog_database)                           },
                           r_norvegicus ={
                             output_data <- CoSIA::annotateDBI(input_id, input_dataset, output_ids, input_species = "homo_sapiens", output_species, 10116, org.Hs.eg.db::org.Hs.eg.db,
                                                               org.Rn.eg.db::org.Rn.eg.db, ortholog_database)                           },
                           stop("Error: Invalid o_species in CoSIAn Object. Make sure the species in the o_species slot is an avalible model 
                                    organism and is in the correct format.")
    )
  }
  else{ # code follows this path if the tool that was inputted into the system does not match.
    stop("Error: Invalid tool in CoSIAn Object. Make sure the tool is either annotationDBI or biomaRt.")
  }
}

#Mus musculus - Mouse

m_musculus<-function(input_id,input_dataset,output_ids,output_species, tool, ortholog_database) { #input data funneling in from the cross species conversion function
  if (tool=="biomaRt"){ # if the user has choose the tool biomart this it the path the codes follows
    Filter_BO_MM <- switch(output_species,
                           c_elegans = {
                             output_data <- CoSIA::biomaRt(input_id, input_dataset, output_ids, input_species = "mus_musculus", output_species, 6239, "mmusculus_gene_ensembl",
                                                        "celegans_gene_ensembl", ortholog_database)                           },
                           d_melanogaster = {
                             output_data <- CoSIA::biomaRt(input_id, input_dataset, output_ids, input_species = "mus_musculus", output_species, 7227, "mmusculus_gene_ensembl",
                                                        "dmelanogaster_gene_ensembl", ortholog_database)                         },
                           m_musculus = {
                             output_data <- CoSIA::biomaRt(input_id, input_dataset, output_ids, input_species = "mus_musculus", output_species, 10090, "mmusculus_gene_ensembl",
                                                        "mmusculus_gene_ensembl", ortholog_database)                           },
                           d_rerio = {
                             output_data <- CoSIA::biomaRt(input_id, input_dataset, output_ids, input_species = "mus_musculus", output_species, 7955, "mmusculus_gene_ensembl",
                                                        "drerio_gene_ensembl", ortholog_database)                               },
                           h_sapien = {
                             output_data <- CoSIA::biomaRt(input_id, input_dataset, output_ids, input_species = "mus_musculus", output_species, 9606, "mmusculus_gene_ensembl",
                                                        "hsapiens_gene_ensembl", ortholog_database)                           },
                           r_norvegicus ={
                             output_data <- CoSIA::biomaRt(input_id, input_dataset, output_ids, input_species = "mus_musculus", output_species, 10116, "mmusculus_gene_ensembl",
                                                        "rnorvegicus_gene_ensembl", ortholog_database)                           },
                           stop("Error: Invalid o_species in CoSIAn Object. Make sure the species in the o_species slot is an avalible model 
                                    organism and is in the correct format.")
    )
  }
  if(tool=="annotationDBI"){# code follows this path if the user chooses annotationDBI as their tool of choice
    Filter_AO_MM <- switch(output_species,
                           h_sapiens = {
                             output_data <- CoSIA::annotateDBI(input_id, input_dataset, output_ids, input_species = "mus_musculus", output_species, 9606, org.Mm.eg.db::org.Mm.eg.db,
                                                               org.Hs.eg.db::org.Hs.eg.db, ortholog_database)                          },
                           m_musculus = {
                             output_data <- CoSIA::annotateDBI(input_id, input_dataset, output_ids, input_species = "mus_musculus", output_species, 10090, org.Mm.eg.db::org.Mm.eg.db,
                                                               org.Mm.eg.db::org.Mm.eg.db, ortholog_database)                          },
                           d_melanogaster = {
                             output_data <- CoSIA::annotateDBI(input_id, input_dataset, output_ids, input_species = "mus_musculus", output_species, 7227, org.Mm.eg.db::org.Mm.eg.db,
                                                               org.Dm.eg.db::org.Dm.eg.db, ortholog_database)                           },
                           danio_rerio = {
                             output_data <- CoSIA::annotateDBI(input_id, input_dataset, output_ids, input_species = "mus_musculus", output_species, 7955, org.Mm.eg.db::org.Mm.eg.db,
                                                               org.Dr.eg.db::org.Dr.eg.db, ortholog_database)                           },
                           c_elegans = {
                             output_data <- CoSIA::annotateDBI(input_id, input_dataset, output_ids, input_species = "mus_musculus", output_species, 6239, org.Mm.eg.db::org.Mm.eg.db,
                                                               org.Ce.eg.db::org.Ce.eg.db, ortholog_database)                           },
                           r_norvegicus ={
                             output_data <- CoSIA::annotateDBI(input_id, input_dataset, output_ids, input_species = "mus_musculus", output_species, 10116, org.Mm.eg.db::org.Mm.eg.db,
                                                               org.Rn.eg.db::org.Rn.eg.db, ortholog_database)                           },
                           stop("Error: Invalid o_species in CoSIAn Object. Make sure the species in the o_species slot is an avalible model 
                                    organism and is in the correct format.")
    )
  }
  else{ # code follows this path if the tool that was inputted into the system does not match.
    stop("Error: Invalid tool in CoSIAn Object. Make sure the tool is either annotationDBI or biomaRt.")
  }
}

#Rattus norvegicus- Rat

r_norvegicus<-function(input_id,input_dataset,output_ids,output_species, tool, ortholog_database) { #input data funneling in from the cross species conversion function
  if (tool=="biomaRt"){ # if the user has choose the tool biomart this it the path the codes follows
    Filter_BO_RN <- switch(output_species,
                           c_elegans = {
                             output_data <- CoSIA::biomaRt(input_id, input_dataset, output_ids, input_species = "r_norvegicus", output_species, 6239, "rnorvegicus_gene_ensembl",
                                                        "celegans_gene_ensembl", ortholog_database)                         },
                           d_melanogaster = {
                             output_data <- CoSIA::biomaRt(input_id, input_dataset, output_ids, input_species = "r_norvegicus", output_species, 7227, "rnorvegicus_gene_ensembl",
                                                        "dmelanogaster_gene_ensembl", ortholog_database)                         },
                           m_musculus = {
                             output_data <- CoSIA::biomaRt(input_id, input_dataset, output_ids, input_species = "r_norvegicus", output_species, 10090, "rnorvegicus_gene_ensembl",
                                                        "mmusculus_gene_ensembl", ortholog_database)                           },
                           d_rerio = {
                             output_data <- CoSIA::biomaRt(input_id, input_dataset, output_ids, input_species = "r_norvegicus", output_species, 7955, "rnorvegicus_gene_ensembl",
                                                        "drerio_gene_ensembl", ortholog_database)                               },
                           h_sapien = {
                             output_data <- CoSIA::biomaRt(input_id, input_dataset, output_ids, input_species = "r_norvegicus", output_species, 9606, "rnorvegicus_gene_ensembl",
                                                        "hsapiens_gene_ensembl", ortholog_database)                          },
                           r_norvegicus ={
                             output_data <- CoSIA::biomaRt(input_id, input_dataset, output_ids, input_species = "r_norvegicus", output_species, 10116, "rnorvegicus_gene_ensembl",
                                                        "rnorvegicus_gene_ensembl", ortholog_database)                           },
                           stop("Error: Invalid o_species in CoSIAn Object. Make sure the species in the o_species slot is an avalible model 
                                    organism and is in the correct format.")
    )
  }
  if(tool=="annotationDBI"){# code follows this path if the user chooses annotationDBI as their tool of choice
    Filter_AO_RN <- switch(output_species,
                           h_sapiens = {
                             output_data <- CoSIA::annotateDBI(input_id, input_dataset, output_ids, input_species = "r_norvegicus", output_species, 9606, org.Rn.eg.db::org.Rn.eg.db,
                                                               org.Hs.eg.db::org.Hs.eg.db, ortholog_database)                          },
                           m_musculus = {
                             output_data <- CoSIA::annotateDBI(input_id, input_dataset, output_ids, input_species = "r_norvegicus", output_species, 10090, org.Rn.eg.db::org.Rn.eg.db,
                                                               org.Mm.eg.db::org.Mm.eg.db, ortholog_database)                          },
                           d_melanogaster = {
                             output_data <- CoSIA::annotateDBI(input_id, input_dataset, output_ids, input_species = "r_norvegicus", output_species, 7227, org.Rn.eg.db::org.Rn.eg.db,
                                                               org.Dm.eg.db::org.Dm.eg.db, ortholog_database)                           },
                           danio_rerio = {
                             output_data <- CoSIA::annotateDBI(input_id, input_dataset, output_ids, input_species = "r_norvegicus", output_species, 7955, org.Rn.eg.db::org.Rn.eg.db,
                                                               org.Dr.eg.db::org.Dr.eg.db, ortholog_database)                           },
                           c_elegans = {
                             output_data <- CoSIA::annotateDBI(input_id, input_dataset, output_ids, input_species = "r_norvegicus", output_species, 6239, org.Rn.eg.db::org.Rn.eg.db,
                                                               org.Ce.eg.db::org.Ce.eg.db, ortholog_database)                           },
                           r_norvegicus ={
                             output_data <- CoSIA::annotateDBI(input_id, input_dataset, output_ids, input_species = "r_norvegicus", output_species, 10116, org.Rn.eg.db::org.Rn.eg.db,
                                                               org.Rn.eg.db::org.Rn.eg.db, ortholog_database)                           },
                           stop("Error: Invalid o_species in CoSIAn Object. Make sure the species in the o_species slot is an avalible model 
                                    organism and is in the correct format.")
    )
  }
  else{ # code follows this path if the tool that was inputted into the system does not match.
    stop("Error: Invalid tool in CoSIAn Object. Make sure the tool is either annotationDBI or biomaRt.")
  }
}

#Tool functions 

#annotationDBI

annotateDBI <- function(input_id, input_dataset, output_ids, input_species, output_species, species_number, input_org, output_org, ortholog_database) {
  # This function uses annotationDBI to convert between gene identifier
  # Switch function below is used to transition the id names into the proper format for using annotationdbi
  input_id <- as.character(input_id)
  output_ids <- as.character(output_ids)
  ID_SWITCH <- Vectorize(vectorize.args = "ids", FUN = function(ids) {
    switch(as.character(ids), 
           Entrez_id = "ENTREZID", 
           Ensembl_id = "ENSEMBL", 
           Ensembl_id_version = "ENSEMBLIDVERSION", 
           Gene_name = "GENENAME",
           Symbol = "SYMBOL")
  })
  
  input_id <- ID_SWITCH(ids = input_id)
  output_ids <- ID_SWITCH(ids = output_ids)
  if (output_species == input_species) {
    # code follows this path if the user chooses the same input species as their output species code follows this path if the user
    # chooses ENSEMBLIDVERSION as their input id
    if (input_id == "ENSEMBLIDVERSION") {
      input_dataset_new <- CoSIA::remove_version_number(input_dataset, input_species)  #puts the genes through the remove version number function to remove the version number (everything after the dot)
      # return(input_dataset_new)
      annotationDbi_data <- AnnotationDbi::select(input_org, keys = as.character(as.matrix(input_dataset_new$ENSEMBL)), columns = output_ids,
                                                  keytype = "ENSEMBL")  # after the version number has been removed run the gene through the annotationDBI package using ensembl id as the gene identifer
      merge <- merge.data.frame(data.frame(annotationDbi_data), data.frame(input_dataset_new), by = "ENSEMBL")  # merge the input and output data together
      return(merge)  #works
    } else {
      # code follows this path if the user chooses an ID other than ENSEMBLIDVERSION as their input id
      annotationDbi_data <- AnnotationDbi::select(input_org, keys = input_dataset, columns = output_ids, keytype = input_id)  #run the code through annotationDBI
      return(annotationDbi_data)  #works
    }
  } else {
    # this code follows this path if the input gene is different than the output gene code follows this path if the user chooses
    # ENSEMBLIDVERSION as their input id
    if (input_id == "ENSEMBLIDVERSION") {
      input_dataset_new <- CoSIA::remove_version_number(input_dataset, input_species)  #puts the genes through the remove version number function to remove the version number (everything after the dot)
      input_dataset_new <- data.frame(input_dataset_new)  #sets the new list of genes without versionIDs as a dataframe
      output_data <- AnnotationDbi::select(input_species, keys = as.character(as.matrix(input_dataset_new)), columns = "ENTREZID", "ENSEMBL",
                                           keytype = "ENSEMBL")  #runs the gene through annotationdbi to get EntrezID conversion in order to conduct ortholog mapping
      output_data <- na.omit(output_data)  #omits the missing values that were not mapped
      merged_data <- merge.data.frame(data.frame(output_data), data.frame(input_dataset_new), by = "ENSEMBL")  #merges the original df with the entrezid mapped data set
      id <- CoSIA::homolog(output_data$ENTREZID, species_number, ortholog_database)  # run the EntrezID found in the new dataset through the homolog function
      id <- data.frame(id)  #makes the homolog output into a dataframe
      names(merged_data)[names(merged_data) == "ENTREZID"] <- "species_one"  #changes the name of the entrezid column to mathc the id dataframe as species one
      merged_data <- merge.data.frame(data.frame(merged_data), data.frame(id), by = "species_one")  # merges the merged data and the homolog values that were found
      colnames(merged_data)[which(names(merged_data) == "ENTREZID")] <- paste(input_species, "Entrez_ID", sep = "_")  #change the names to more formal names for the columns
      colnames(merged_data)[which(names(merged_data) == "ENSEMBL")] <- paste(input_species, "Ensembl_ID", sep = "_")  #change the names to more formal names for the columns
      colnames(merged_data)[which(names(merged_data) == "SYMBOL")] <- paste(input_species, "Symbol", sep = "_")  #change the names to more formal names for the columns
      colnames(merged_data)[which(names(merged_data) == "ENSEMBLIDVERSION")] <- paste(input_species, "Ensembl_ID_with_Version_ID", sep = "_")  #change the names to more formal names for the columns
      # names(output_data)[names(output_data) == 'ENTREZID'] <- 'species_one' #change the name of the output dataframe to species_one
      # for merging later
      ortho <- AnnotationDbi::select(output_org, keys = as.character(as.matrix(merged_data$species_two)), columns = c(output_ids, "ENTREZID"),
                                     keytype = "ENTREZID")  # next we are going to take the entrezids of the orthologs and return the output values of the gene
      non_na <- which(is.na(ortho$ENTREZID) == FALSE)  # determine the indices for the non-NA genes
      ortho <- ortho[non_na, ]  #return only the genes with annotations using indices
      non_dups <- which(duplicated(ortho$ENTREZID) == FALSE)  #determine the indices for the non-duplicated genes
      ortho <- ortho[non_dups, ]  #return only the non-duplicated genes using indices
      ortho <- data.frame(ortho)  #make this into a dataframe
      names(ortho)[names(ortho) == "ENTREZID"] <- "species_two"  # set the entrezID into a dataframe
      merged_data <- merge.data.frame(data.frame(ortho), data.frame(merged_data), by = "species_two")  #merge the ortholog conversion dataframe with the merge dataframe
      # names(merged_data)[names(merged_data) == 'species_two'] <- paste(output_species,'Entrez_ID',sep='_') #changes name to more
      # formal names names(merged_data)[names(merged_mouse_data) == 'species_one'] <- paste(input_species,'Entrez_ID',sep='_')
      # #changes name to more formal names
      colnames(merged_data)[which(names(merged_data) == "ENTREZID")] <- paste(output_species, "Entrez_ID", sep = "_")  #changes name to more formal names
      colnames(merged_data)[which(names(merged_data) == "ENSEMBL")] <- paste(output_species, "Ensembl_ID", sep = "_")  #changes name to more formal names
      colnames(merged_data)[which(names(merged_data) == "SYMBOL")] <- paste(output_species, "Symbol", sep = "_")  #changes name to more formal names
      colnames(merged_data)[which(names(merged_data) == "ENSEMBLIDVERSION")] <- paste(output_species, "Ensembl_ID_with_Version_ID", sep = "_")  #changes name to more formal names
      merged_data <- merged_data[!duplicated(as.list(merged_data))]  # removes duplicates from the data
      merged_data = subset(merged_data, select = -c(species_two, species_one))  # removes the species two and species one parts of the data
      return(merged_data)  #works
    } else {
      # code follows this path if the user chooses an ID other than ENSEMBLIDVERSION as their input id
      output_data <- AnnotationDbi::select(input_org, keys = as.character(as.matrix(input_dataset)), columns = "ENTREZID", input_id, keytype = input_id)  #runs the gene through annotationdbi to get EntrezID conversion in order to conduct ortholog mapping
      output_data <- na.omit(output_data)  # omits the NA values
      id <- CoSIA::homolog(output_data$ENTREZID, species_number, ortholog_database)  #send the entrezid ids through the homolog conversion
      id <- data.frame(id)  # make these ortholog ids into a new dataframe
      names(output_data)[names(output_data) == "ENTREZID"] <- "species_one"  # make the EntrezID that were used for the conversion into species one
      merged_data <- merge.data.frame(data.frame(output_data), data.frame(id), by = "species_one")  # merge the two dataframes with the original ids and the ortholog ids together
      colnames(merged_data)[which(names(merged_data) == "ENTREZID")] <- paste(input_species, "Entrez_ID", sep = "_")  #changes name to more formal names
      colnames(merged_data)[which(names(merged_data) == "ENSEMBL")] <- paste(input_species, "Ensembl_ID", sep = "_")  #changes name to more formal names
      colnames(merged_data)[which(names(merged_data) == "SYMBOL")] <- paste(input_species, "Symbol", sep = "_")  #changes name to more formal names
      colnames(merged_data)[which(names(merged_data) == "ENSEMBLIDVERSION")] <- paste(input_species, "Ensembl_ID_with_Version_ID", sep = "_")  #changes name to more formal names
      ortho <- AnnotationDbi::select(output_org, keys = as.character(as.matrix(merged_data$species_two)), columns = output_ids, "ENTREZID",
                                     keytype = "ENTREZID")  # now that you have the entrezids of the converted species we can convert to different gene identifiers
      non_na <- which(is.na(ortho$ENTREZID) == FALSE)  # Determine the indices for the non-NA genes
      ortho <- ortho[non_na, ]  # Return only the genes with annotations using indices
      non_dups <- which(duplicated(ortho$ENTREZID) == FALSE)  # Determine the indices for the non-duplicated genes
      ortho <- ortho[non_dups, ]  # Return only the non-duplicated genes using indices
      ortho <- data.frame(ortho)  # make the ortholog conversion into a dataframe
      names(ortho)[names(ortho) == "ENTREZID"] <- "species_two"  #rename the species into species two in order to merge the orthologs and the original data
      merged_data <- merge.data.frame(data.frame(ortho), data.frame(merged_data), by = "species_two")  # merge
      colnames(merged_data)[names(merged_data) == "species_two"] <- paste(output_species, "Entrez_ID", sep = "_")  #changes name to more formal names
      # names(merged_data)[names(merged_data) == 'species_one'] <- paste(input_species,'Entrez_ID',sep='_') #changes name to more
      # formal names
      colnames(merged_data)[which(names(merged_data) == "ENTREZID")] <- paste(output_species, "Entrez_ID", sep = "_")  #changes name to more formal names
      colnames(merged_data)[which(names(merged_data) == "ENSEMBL")] <- paste(output_species, "Ensembl_ID", sep = "_")  #changes name to more formal names
      colnames(merged_data)[which(names(merged_data) == "SYMBOL")] <- paste(output_species, "Symbol", sep = "_")  #changes name to more formal names
      colnames(merged_data)[which(names(merged_data) == "ENSEMBLIDVERSION")] <- paste(output_species, "Ensembl_ID_with_Version_ID", sep = "_")  #changes name to more formal names
      merged_data <- merged_data[!duplicated(as.list(merged_data))]  # remove duplicates
      merged_data = subset(merged_data, select = -c(species_one))  #remove species one
      return(merged_data)  #works
    }
  }
}

#biomaRT

biomaRt <- function(input_id, input_dataset, output_ids, input_species, output_species, species_number, species_dataset, output_species_dataset,
                 ortholog_database) {
  # biomart
  input_id <- as.character(input_id)
  output_ids <- as.character(output_ids)
  ID_SWITCH <- Vectorize(vectorize.args = "ids", FUN = function(ids) {
    switch(as.character(ids), 
           Entrez_id = "entrezgene_id", 
           Ensembl_id = "ensembl_gene_id", 
           Ensembl_id_version = "ensembl_gene_id_version",
           Gene_name = "external_gene_name", 
           MGI_Symbol = "mgi_symbol", 
           HGNC_Symbol = "hgnc_symbol")
  })
  input_id <- ID_SWITCH(ids = input_id)
  output_ids <- ID_SWITCH(ids = output_ids)
  if (input_species == output_species) {
    # goes through this path if the input and output species are the same
    mart <- biomaRt::useMart("ensembl", dataset = species_dataset)  # pulls the biomaRt object for the species species that has been choosen
    attributes <- c(output_ids, input_id)  # sets the attributes that the user wants makes sure to set both the input and output values
    filters = input_id  # sets the filters as the input vales
    output_data <- biomaRt::getBM(attributes = attributes, filters = filters, values = input_dataset, mart = mart, uniqueRows = TRUE, bmHeader = FALSE)  # run conversion through biomaRt
    return(output_data)  # returnt he biomaRt output
  } else {
    # goes through this path if the input and output species are different
    mart <- biomaRt::useMart("ensembl", dataset = species_dataset)  # sets up biomart for species input
    attributes <- c("entrezgene_id", input_id)  # sets input ids and entrezids as attributes (output values)
    filters <- input_id  # set the input ids as the filter (input values)
    output_data <- biomaRt::getBM(attributes = attributes, filters = filters, values = input_dataset, mart = mart, uniqueRows = TRUE, bmHeader = FALSE)  # run the convesion through biomaRt
    output_data <- na.omit(output_data)  #omit the NA values
    id <- CoSIA::homolog(output_data$entrezgene_id, species_number, ortholog_database)  #now that we have the entrezids for the input species we can run it through the homolog function and then get the entrezids for the other species
    id <- data.frame(id)  #this makes sure that our homologous gene list is in a dataframe
    names(output_data)[names(output_data) == "entrezgene_id"] <- "species_one"  # change the name so we can merge the original and the homologous dataframes together
    merged_data <- merge.data.frame(data.frame(output_data), data.frame(id), by = "species_one")  #merge the two dataframes
    colnames(merged_data)[which(names(merged_data) == "entrezgene_id")] <- paste(input_species, "Entrez_ID", sep = "_")  #rename to a more formal name
    colnames(merged_data)[which(names(merged_data) == "ensembl_gene_id")] <- paste(input_species, "Ensembl_ID", sep = "_")  #rename to a more formal name
    colnames(merged_data)[which(names(merged_data) == "external_gene_name")] <- paste(input_species, "Gene_Name", sep = "_")  #rename to a more formal name
    colnames(merged_data)[which(names(merged_data) == "hgnc_symbol")] <- paste(input_species, "HGNC_Symbol", sep = "_")  #rename to a more formal name
    colnames(merged_data)[which(names(merged_data) == "mgi_symbol")] <- paste(input_species, "MGI_Symbol", sep = "_")  #rename to a more formal name
    colnames(merged_data)[which(names(merged_data) == "ensembl_gene_id_version")] <- paste(input_species, "Ensembl_ID_with_Version_ID", sep = "_")  #rename to a more formal name
    # names(output_data)[names(output_data) == 'entrezgene_id'] <- 'species_one'
    marts <- biomaRt::useMart("ensembl", dataset = output_species_dataset)  #set the biomart species to the new species
    ortho <- biomaRt::getBM(attributes = c(output_ids, "entrezgene_id"), filters = "entrezgene_id", values = id$species_two, mart = marts,
                            uniqueRows = TRUE, bmHeader = FALSE)  # run the homologous entrezids to the other gene identifiers
    non_na <- which(is.na(ortho$entrezgene_id) == FALSE)  #pick non na values from ortholog conversion
    ortho <- ortho[non_na, ]  #subset only nonna values
    non_dups <- which(duplicated(ortho$entrezgene_id) == FALSE)  #remove duplicates
    ortho <- ortho[non_dups, ]  #subset only duplicates
    ortho <- data.frame(ortho)  # set the ortho variable as a dataframe
    names(ortho)[names(ortho) == "entrezgene_id"] <- "species_two"  # rename the entrezid as species two
    merged_species_data <- merge.data.frame(data.frame(ortho), data.frame(merged_data), by = "species_two")  # merge the ortholog dataframe and the original data
    names(merged_species_data)[names(merged_species_data) == "species_two"] <- paste(output_species, "Entrez_ID", sep = "_")  # clean up names
    names(merged_species_data)[names(merged_species_data) == "species_one"] <- paste(input_species, "Entrez_ID", sep = "_")  #clean up names
    colnames(merged_species_data)[which(names(merged_species_data) == "entrezgene_id")] <- paste(output_species, "Entrez_ID", sep = "_")  #clean up names
    colnames(merged_species_data)[which(names(merged_species_data) == "ensembl_gene_id")] <- paste(output_species, "Ensembl_ID", sep = "_")  #clean up names
    colnames(merged_species_data)[which(names(merged_species_data) == "hgnc_symbol")] <- paste(output_species, "HGNC_Symbol", sep = "_")  #clean up names
    colnames(merged_species_data)[which(names(merged_species_data) == "mgi_symbol")] <- paste(output_species, "MGI_Symbol", sep = "_")  #clean up names
    colnames(merged_species_data)[which(names(merged_species_data) == "external_gene_name")] <- paste(output_species, "Gene_Name", sep = "_")  #clean up names
    colnames(merged_species_data)[which(names(merged_species_data) == "ensembl_gene_id_version")] <- paste(output_species, "Ensembl_ID_with_Version_ID",
                                                                                                           sep = "_")  #clean up names
    merged_species_data <- merged_species_data[!duplicated(as.list(merged_species_data))]  #remove duplicates
    return(merged_species_data)  # return the final table
  }
}

#homolog files

homolog <- function(entrez_data, species_number, ortholog_database) {
  if (ortholog_database == "HomoloGene") {
    myGenes <- as.character(entrez_data)
    species_one <- entrez_data
    homologene <- homologene::homologeneData2
    data <- annotationTools::getHOMOLOG(myGenes, species_number, homologene, noIDsymbol = NA, clusterCol = 1, speciesCol = 3, idCol = 4)
    data <- as.list(data)
    species_two <- as.character(data)
    homologs <- as.data.frame(species_two)
    homologs <- merge(data.frame(homologs, row.names = NULL), data.frame(species_one, row.names = NULL), by = 0, all = TRUE)[-1]
    manipulated_homologs <- data.frame(homologs)
    x <- 1:nrow(manipulated_homologs)
    for (i in seq_along(x)) {
      ## print(i)
      char <- manipulated_homologs[i, 1]
      ## print(char)
      if ((char == "NA") == TRUE) {
        manipulated_homologs[i, 1] <- NA
      }
    }
    species_two_entrez_ID <- na.omit(manipulated_homologs)
    return(species_two_entrez_ID)
  }
  if (ortholog_database == "NCBIOrtho") {
    myGenes <- as.vector(entrez_data)
    species_one <- entrez_data
    NCBIOrtho <- readr::read_table("https://ftp.ncbi.nih.gov/gene/DATA/gene_orthologs.gz", col_types = "ccccc", show_col_types = TRUE)
    NCBIOrtho <- as.data.frame(NCBIOrtho)
    data <- annotationTools::getHOMOLOG(myGenes, species_number, NCBIOrtho, tableType = "gene_orthologs")
    data <- as.list(data)
    species_two <- as.character(data)
    homologs <- as.data.frame(species_two)
    homologs <- merge(data.frame(homologs, row.names = NULL), data.frame(species_one, row.names = NULL), by = 0, all = TRUE)[-1]
    manipulated_homologs <- data.frame(homologs)
    x <- 1:nrow(manipulated_homologs)
    for (i in seq_along(x)) {
      ## print(i)
      char <- manipulated_homologs[i, 1]
      ## print(char)
      if ((char == "NA") == TRUE) {
        manipulated_homologs[i, 1] <- NA
      }
    }
    species_two_entrez_ID <- na.omit(manipulated_homologs)
    return(species_two_entrez_ID)
  }
}













