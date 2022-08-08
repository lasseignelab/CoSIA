#' getConversion Method
#'
#' @param object CoSIAn. 
#'
#' @return 
#' @export
#'
#' @examples
setMethod("getConversion", signature(object = "CoSIAn"), function(object) { # user's input of the function
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
    return(species_data)
})

# getConversion Functions

# Species Functions

c_elegans<-function(input_id,input_dataset,output_ids,output_species, tool, ortholog_database) { #input data funneling in from the cross species conversion function
  if (tool=="biomaRt"){ # if the user has choose the tool biomart this it the path the codes follows
    Filter_BO_CE <- switch(output_species,
                               c_elegans = {
                                 output_data<-CoSIA::BioM(input_id,input_dataset, output_ids, input_species="c_elegans", output_species, 6239, "celegans_gene_ensembl", "celegans_gene_ensembl", ortholog_database)
                               },
                               d_melanogaster = {
                                 output_data<-CoSIA::BioM(input_id,input_dataset, output_ids, input_species="c_elegans", output_species, 7227, "celegans_gene_ensembl", "dmelanogaster_gene_ensembl", ortholog_database)
                               }
                               m_musculus = {
                                 output_data<-CoSIA::BioM(input_id,input_dataset, output_ids,input_species="c_elegans", output_species, 10090, "celegans_gene_ensembl", "mmusculus_gene_ensembl", ortholog_database)
                               }
                               d_rerio = {
                                 output_data<-CoSIA::BioM(input_id,input_dataset, output_ids, input_species="c_elegans",output_species, 7955, "celegans_gene_ensembl", "drerio_gene_ensembl", ortholog_database)
                               }
                               h_sapien = {
                                 output_data<-CoSIA::BioM(input_id,input_dataset, output_ids,input_species="c_elegans", output_species, 9606, "celegans_gene_ensembl", "hsapiens_gene_ensembl", ortholog_database)
                               }
                               r_norvegicus ={
                                 output_data<-CoSIA::BioM(input_id,input_dataset, output_ids,input_species="c_elegans", output_species, 10116, "celegans_gene_ensembl", "rnorvegicus_gene_ensembl", ortholog_database)
                               },
                               stop("Error: Invalid o_species in CoSIAn Object. Make sure the species in the o_species slot is an avalible model 
                                    organism and is in the correct format.")
    )
  }
  if(tool=="annotationDBI"){# code follows this path if the user chooses annotationDBI as their tool of choice
    Filter_AO_CE <- switch(output_species,
                           h_sapiens = {
                             output_data<-CoSIA::AnnotateDBI(input_id,input_dataset,output_ids,input_species="c_elegans",output_species,9606, org.Ce.eg.db::org.Ce.eg.db , org.Hs.eg.db::org.Hs.eg.db, ortholog_database)
                             },
                           m_musculus = {
                             output_data<-CoSIA::AnnotateDBI(input_id,input_dataset,output_ids,input_species="c_elegans",output_species,10090, org.Ce.eg.db::org.Ce.eg.db , org.Mm.eg.db::org.Mm.eg.db, ortholog_database)
                           },
                           d_melanogaster = {
                             output_data<-CoSIA::AnnotateDBI(input_id,input_dataset,output_ids,input_species="c_elegans",output_species,7227, org.Ce.eg.db::org.Ce.eg.db , org.Dm.eg.db::org.Dm.eg.db, ortholog_database)
                           },
                           danio_rerio = {
                             output_data<-CoSIA::AnnotateDBI(input_id,input_dataset,output_ids,input_species="c_elegans",output_species,7955, org.Ce.eg.db::org.Ce.eg.db , org.Dr.eg.db::org.Dr.eg.db, ortholog_database)
                           },
                           c_elegans = {
                             output_data<-CoSIA::AnnotateDBI(input_id,input_dataset,output_ids,input_species="c_elegans",output_species,6239, org.Ce.eg.db::org.Ce.eg.db , org.Ce.eg.db::org.Ce.eg.db, ortholog_database)
                           }
                           r_norvegicus ={
                             output_data<-CoSIA::AnnotateDBI(input_id,input_dataset,output_ids,input_species="c_elegans",output_species,10116, org.Ce.eg.db::org.Ce.eg.db , org.Rn.eg.db::org.Rn.eg.db,ortholog_database)
                           },
                           stop("Error: Invalid o_species in CoSIAn Object. Make sure the species in the o_species slot is an avalible model 
                                    organism and is in the correct format.")
                           )
  }
  else{ # code follows this path if the tool that was inputted into the system does not match.
    stop("Error: Invalid tool in CoSIAn Object. Make sure the tool is either annotationDBI or biomaRt.")
  }
}

d_melanogaster<-function(input_id,input_dataset,output_ids,output_species, tool, ortholog_database) { #input data funneling in from the cross species conversion function
  if (tool=="biomaRt"){ # if the user has choose the tool biomart this it the path the codes follows
    Filter_BO_DM <- switch(output_species,
                           c_elegans = {
                             output_data<-CoSIA::BioM(input_id,input_dataset, output_ids, input_species="d_melanogaster",output_species, 6239, "dmelanogaster_gene_ensembl", "celegans_gene_ensembl", ortholog_database)
                           }
                           d_melanogaster = {
                             output_data<-CoSIA::BioM(input_id,input_dataset, output_ids,input_species="d_melanogaster", output_species, 7227, "dmelanogaster_gene_ensembl", "dmelanogaster_gene_ensembl", ortholog_database)
                           }
                           m_musculus = {
                             output_data<-CoSIA::BioM(input_id,input_dataset, output_ids,input_species="d_melanogaster", output_species, 10090, "dmelanogaster_gene_ensembl", "mmusculus_gene_ensembl", ortholog_database)
                           }
                           d_rerio = {
                             output_data<-CoSIA::BioM(input_id,input_dataset, output_ids,input_species="d_melanogaster", output_species, 7955, "dmelanogaster_gene_ensembl", "drerio_gene_ensembl", ortholog_database)
                           }
                           h_sapien = {
                             output_data<-CoSIA::BioM(input_id,input_dataset, output_ids,input_species="d_melanogaster", output_species, 9606, "dmelanogaster_gene_ensembl", "hsapiens_gene_ensembl", ortholog_database)
                           }
                           r_norvegicus ={
                             output_data<-CoSIA::BioM(input_id,input_dataset, output_ids,input_species="d_melanogaster", output_species, 10116, "dmelanogaster_gene_ensembl", "rnorvegicus_gene_ensembl", ortholog_database)
                           },
                           stop("Error: Invalid o_species in CoSIAn Object. Make sure the species in the o_species slot is an avalible model 
                                    organism and is in the correct format.")
    )
  }
  if(tool=="annotationDBI"){# code follows this path if the user chooses annotationDBI as their tool of choice
    Filter_AO_DM <- switch(output_species,
                           h_sapiens = {
                             output_data<-CoSIA::AnnotateDBI(input_id,input_dataset,output_ids,input_species="d_melanogaster",output_species,9606, org.Dm.eg.db::org.Dm.eg.db , org.Hs.eg.db::org.Hs.eg.db, ortholog_database)
                           },
                           m_musculus = {
                             output_data<-CoSIA::AnnotateDBI(input_id,input_dataset,output_ids,input_species="d_melanogaster",output_species,10090, org.Dm.eg.db::org.Dm.eg.db , org.Mm.eg.db::org.Mm.eg.db, ortholog_database)
                           },
                           d_melanogaster = {
                             output_data<-CoSIA::AnnotateDBI(input_id,input_dataset,output_ids,input_species="d_melanogaster",output_species,7227, org.Dm.eg.db::org.Dm.eg.db , org.Dm.eg.db::org.Dm.eg.db, ortholog_database)
                           },
                           danio_rerio = {
                             output_data<-CoSIA::AnnotateDBI(input_id,input_dataset,output_ids,input_species="d_melanogaster",output_species,7955, org.Dm.eg.db::org.Dm.eg.db , org.Dr.eg.db::org.Dr.eg.db, ortholog_database)
                           },
                           c_elegans = {
                             output_data<-CoSIA::AnnotateDBI(input_id,input_dataset,output_ids,input_species="d_melanogaster",output_species,6239, org.Dm.eg.db::org.Dm.eg.db , org.Ce.eg.db::org.Ce.eg.db, ortholog_database)
                           }
                           r_norvegicus ={
                             output_data<-CoSIA::AnnotateDBI(input_id,input_dataset,output_ids,input_species="d_melanogaster",output_species,10116, org.Dm.eg.db::org.Dm.eg.db , org.Rn.eg.db::org.Rn.eg.db, ortholog_database)
                           },
                           stop("Error: Invalid o_species in CoSIAn Object. Make sure the species in the o_species slot is an avalible model 
                                    organism and is in the correct format.")
    )
  }
  else{ # code follows this path if the tool that was inputted into the system does not match.
    stop("Error: Invalid tool in CoSIAn Object. Make sure the tool is either annotationDBI or biomaRt.")
  }
}





