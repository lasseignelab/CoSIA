mus_musculus<-function(input_id,input_dataset,output_ids,output_species, tool, ortholog_database) {
  if (tool=="biomaRt"){ # code follows this path if the user chooses biomaRt as their tool of choice{
    if (output_species=="danio_rerio"){ #code follows this path if the user chooses homo sapiens as their output species
      output_data<-CoSIA::BioM(input_id,input_dataset, output_ids, input_species="mus_musculus", output_species, 7955, "mmusculus_gene_ensembl", "drerio_gene_ensembl", ortholog_database)
      return(output_data)
    }
    if (output_species=="d_melanogaster"){ #code follows this path if the user chooses homo sapiens as their output species
      output_data<-CoSIA::BioM(input_id,input_dataset, output_ids, input_species="mus_musculus", output_species, 7227, "mmusculus_gene_ensembl", "dmelanogaster_gene_ensembl", ortholog_database)
      return(output_data)
    }
    if(output_species=="mus_musculus"){ #code follows this path if the user chooses mus musculus as their output species
      output_data<-CoSIA::BioM(input_id,input_dataset, output_ids,input_species="mus_musculus", output_species, 10090, "mmusculus_gene_ensembl", "mmusculus_gene_ensembl", ortholog_database)
      return(output_data)
    }
    if(output_species=="c_elegans"){
      output_data<-CoSIA::BioM(input_id,input_dataset, output_ids, input_species="mus_musculus",output_species, 6239, "mmusculus_gene_ensembl", "celegans_gene_ensembl", ortholog_database)
      return(output_data)
    }
    if(output_species=="homo_sapiens"){
      output_data<-CoSIA::BioM(input_id,input_dataset, output_ids,input_species="mus_musculus", output_species, 9606, "mmusculus_gene_ensembl", "hsapiens_gene_ensembl", ortholog_database)
      return(output_data)
    }
    if(output_species=="r_norvegicus"){
      output_data<-CoSIA::BioM(input_id,input_dataset, output_ids,input_species="mus_musculus", output_species, 10116, "mmusculus_gene_ensembl", "rnorvegicus_gene_ensembl", ortholog_database)
      return(output_data)
    }
    else{
      stop("Error. Invalid output species. Make sure it matches the proper format")
    }
  }
  if(tool=="annotationDBI"){# code follows this path if the user chooses annotationDBI as their tool of choice{
    if (output_species=="homo_sapiens"){ #code follows this path if the user chooses homo sapiens as their output species
      output_data<-CoSIA::AnnotateDBI(input_id,input_dataset,output_ids,input_species="mus_musculus",output_species,9606, org.Mm.eg.db::org.Mm.eg.db , org.Hs.eg.db::org.Hs.eg.db, ortholog_database)
      return(output_data)
    }
    if (output_species=="d_melanogaster"){ #code follows this path if the user chooses homo sapiens as their output species
      output_data<-CoSIA::AnnotateDBI(input_id,input_dataset,output_ids,input_species="mus_musculus",output_species,7227, org.Mm.eg.db::org.Mm.eg.db , org.Dm.eg.db::org.Dm.eg.db, ortholog_database)
      return(output_data)
    }
    if (output_species=="mus_musculus"){ #code follows this path if the user chooses homo sapiens as their output species
      output_data<-CoSIA::AnnotateDBI(input_id,input_dataset,output_ids,input_species="mus_musculus",output_species,10090, org.Mm.eg.db::org.Mm.eg.db , org.Mm.eg.db::org.Mm.eg.db, ortholog_database)
      return(output_data)
    }
    if (output_species=="danio_rerio"){ #code follows this path if the user chooses homo sapiens as their output species
      output_data<-CoSIA::AnnotateDBI(input_id,input_dataset,output_ids,input_species="mus_musculus",output_species,7955, org.Mm.eg.db::org.Mm.eg.db , org.Dr.eg.db::org.Dr.eg.db, ortholog_database)
      return(output_data)
    }
    if (output_species=="c_elegans"){ #code follows this path if the user chooses homo sapiens as their output species
      output_data<-CoSIA::AnnotateDBI(input_id,input_dataset,output_ids,input_species="mus_musculus",output_species,6239, org.Mm.eg.db::org.Mm.eg.db , org.Ce.eg.db::org.Ce.eg.db, ortholog_database)
      return(output_data)
    }
    if (output_species=="r_norvegicus"){ #code follows this path if the user chooses homo sapiens as their output species
      output_data<-CoSIA::AnnotateDBI(input_id,input_dataset,output_ids,input_species="mus_musculus",output_species,10116, org.Mm.eg.db::org.Mm.eg.db , org.Rn.eg.db::org.Rn.eg.db, ortholog_database)
      return(output_data)
    }
    else{
      stop("Error. Invalid output species. Make sure it matches the proper format")
    }
  }
  else{
    stop("Error. Invalid tool. Make sure it matches the proper format")
  }
}