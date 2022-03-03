c_elegans<-function(input_id,input_dataset,output_ids,output_species, tool, ortholog_database) { #input data funneling in from the cross species conversion function
  if (tool=="biomaRt"){ # if the user has choosen the tool biomart this it the path the codes follows
    if (output_species=="c_elegans"){ #code follows this path if the user chooses c_elegans as their output species
      output_data<-CoSIA::BioM(input_id,input_dataset, output_ids, input_species="c_elegans", output_species, 6239, "celegans_gene_ensembl", "celegans_gene_ensembl", ortholog_database)
      return(output_data)
    }
    if (output_species=="d_melanogaster"){ #code follows this path if the user chooses d_melanogaster as their output species
      output_data<-CoSIA::BioM(input_id,input_dataset, output_ids, input_species="c_elegans", output_species, 7227, "celegans_gene_ensembl", "dmelanogaster_gene_ensembl", ortholog_database)
      return(output_data)
    }
    if(output_species=="mus_musculus"){ #code follows this path if the user chooses mus_musculus as their output species
      output_data<-CoSIA::BioM(input_id,input_dataset, output_ids,input_species="c_elegans", output_species, 10090, "celegans_gene_ensembl", "mmusculus_gene_ensembl", ortholog_database)
      return(output_data)
    }
    if(output_species=="danio_rerio"){ # code follows this path if the user chooses danio_rerio
      output_data<-CoSIA::BioM(input_id,input_dataset, output_ids, input_species="c_elegans",output_species, 7955, "celegans_gene_ensembl", "drerio_gene_ensembl", ortholog_database)
      return(output_data)
    }
    if(output_species=="homo_sapien"){ # code follows this path if the user chooses homo_sapien
      output_data<-CoSIA::BioM(input_id,input_dataset, output_ids,input_species="c_elegans", output_species, 9606, "celegans_gene_ensembl", "hsapiens_gene_ensembl", ortholog_database)
      return(output_data)
    }
    if(output_species=="r_norvegicus"){ # code follows this path if the user chooses r_norvegicus
      output_data<-CoSIA::BioM(input_id,input_dataset, output_ids,input_species="c_elegans", output_species, 10116, "celegans_gene_ensembl", "rnorvegicus_gene_ensembl", ortholog_database)
      return(output_data)
    }
    else{ # code follows this path if the output species is written incorrectly
      stop("Error. Invalid output species. Make sure it matches the proper format")
    }
  }
  if(tool=="annotationDBI"){# code follows this path if the user chooses annotationDBI as their tool of choice
    if (output_species=="homo_sapiens"){ #code follows this path if the user chooses homo sapiens as their output species
      output_data<-CoSIA::AnnotateDBI(input_id,input_dataset,output_ids,input_species="c_elegans",output_species,9606, org.Ce.eg.db::org.Ce.eg.db , org.Hs.eg.db::org.Hs.eg.db, ortholog_database)
      return(output_data)
    }
    if (output_species=="mus_musculus"){ #code follows this path if the user chooses mus_musculus as their output species
      output_data<-CoSIA::AnnotateDBI(input_id,input_dataset,output_ids,input_species="c_elegans",output_species,10090, org.Ce.eg.db::org.Ce.eg.db , org.Mm.eg.db::org.Mm.eg.db, ortholog_database)
      return(output_data)
    }
    if (output_species=="d_melanogaster"){ #code follows this path if the user chooses d_melanogaster as their output species
      output_data<-CoSIA::AnnotateDBI(input_id,input_dataset,output_ids,input_species="c_elegans",output_species,7227, org.Ce.eg.db::org.Ce.eg.db , org.Dm.eg.db::org.Dm.eg.db, ortholog_database)
      return(output_data)
    }
    if (output_species=="danio_rerio"){ #code follows this path if the user chooses danio_rerio as their output species
      output_data<-CoSIA::AnnotateDBI(input_id,input_dataset,output_ids,input_species="c_elegans",output_species,7955, org.Ce.eg.db::org.Ce.eg.db , org.Dr.eg.db::org.Dr.eg.db, ortholog_database)
      return(output_data)
    }
    if (output_species=="c_elegans"){ #code follows this path if the user chooses c_elegans as their output species
      output_data<-CoSIA::AnnotateDBI(input_id,input_dataset,output_ids,input_species="c_elegans",output_species,6239, org.Ce.eg.db::org.Ce.eg.db , org.Ce.eg.db::org.Ce.eg.db, ortholog_database)
      return(output_data)
    }
    if (output_species=="r_norvegicus"){ #code follows this path if the user chooses r_norvegicus as their output species
      output_data<-CoSIA::AnnotateDBI(input_id,input_dataset,output_ids,input_species="c_elegans",output_species,10116, org.Ce.eg.db::org.Ce.eg.db , org.Rn.eg.db::org.Rn.eg.db,ortholog_database)
      return(output_data)
    }
    else{ #code follows this path if the output species that the user inputted does not match any of the species that are available
      stop("Error. Invalid output species. Make sure it matches the proper format.")
    }
  }
  else{ # code follows this path if the tool that was inputted into the system does not match.
    stop("Error. Invalid tool. Make sure it matches the proper format")
  }
}