BioM<-function(input_id,input_dataset,output_ids,input_species,output_species,species_number, species_dataset, output_species_dataset, ortholog_database){
  #biomart
  input_id <- switch(input_id,
                     "Entrez.id"="entrezgene_id",
                     "Ensembl.id"="ensembl_gene_id",
                     "Ensembl.id.version"="ensembl_gene_id_version",
                     "Gene.name"="external_gene_name",
                     "Symbol"="symbol")
  if(input_id=="Symbol"){
    if(input_species=="mus_musculus")
      input_id <- "mgi_symbol"
    if(input_species=="homo_sapiens")
      input_id <- "hgnc_symbol"
  }
  output_ids <- output_ids
  for(x in seq(length(output_ids))){
    output_ids[x] <- switch(output_ids[x],
                            "Entrez.id"="entrezgene_id",
                            "Ensembl.id"="ensembl_gene_id",
                            "Ensembl.id.version"="ensembl_gene_id_version",
                            "Gene.name"="external_gene_name",
                            "Symbol"="Symbol",
                            "DEFAULT")
    if(output_ids[x]=="Symbol"){
      if(input_species=="mus_musculus")
        output_ids[x] <- "mgi_symbol"
      if(input_species=="homo_sapiens")
        output_ids[x] <- "hgnc_symbol"
    }
  }
  if (input_species==output_species){ # goes through this path if the input and output species are the same
    mart<- biomaRt::useMart("ensembl", dataset= species_dataset) # pulls the biomaRt object for the species species that has been choosen
    attributes <- c(output_ids,input_id) # sets the attributes that the user wants makes sure to set both the input and output values
    filters = input_id # sets the filters as the input vales
    output_data <- biomaRt::getBM(attributes=attributes, filters= filters,values= input_dataset, mart=mart, uniqueRows=TRUE, bmHeader=FALSE) # run conversion through biomaRt
    return(output_data) # returnt he biomaRt output
  }
  else{ #goes through this path if the input and output species are different
    mart<- biomaRt::useMart("ensembl", dataset= species_dataset) # sets up biomart for species input
    attributes <- c("entrezgene_id",input_id) # sets input ids and entrezids as attributes (output values)
    filters <- input_id # set the input ids as the filter (input values)
    output_data <- biomaRt::getBM(attributes=attributes, filters= filters,values=input_dataset, mart=mart, uniqueRows=TRUE, bmHeader=FALSE)# run the convesion through biomaRt
    output_data<-na.omit(output_data) #omit the NA values
    id <- CoSIA::homolog(output_data$entrezgene_id, species_number, ortholog_database) #now that we have the entrezids for the input species we can run it through the homolog function and then get the entrezids for the other species
    id<-data.frame(id) #this makes sure that our homologous gene list is in a dataframe
    names(output_data)[names(output_data) == "entrezgene_id"] <- "species_one" # change the name so we can merge the original and the homologous dataframes together
    merged_data<-merge.data.frame(data.frame(output_data), data.frame(id), by = "species_one") #merge the two dataframes
    colnames(merged_data)[which(names(merged_data) == "entrezgene_id")] <- paste(input_species,"Entrez_ID",sep="_") #rename to a more formal name
    colnames(merged_data)[which(names(merged_data) == "ensembl_gene_id")] <- paste(input_species,"Ensembl_ID",sep="_")#rename to a more formal name
    colnames(merged_data)[which(names(merged_data) == "external_gene_name")] <- paste(input_species,"Gene_Name",sep="_")#rename to a more formal name
    colnames(merged_data)[which(names(merged_data) == "hgnc_symbol")] <- paste(input_species,"HGNC_Symbol",sep="_")#rename to a more formal name
    colnames(merged_data)[which(names(merged_data) == "mgi_symbol")] <- paste(input_species,"MGI_Symbol",sep="_")#rename to a more formal name
    colnames(merged_data)[which(names(merged_data) == "ensembl_gene_id_version")] <- paste(input_species,"Ensembl_ID_with_Version_ID",sep="_")#rename to a more formal name
    #names(output_data)[names(output_data) == "entrezgene_id"] <- "species_one"
    marts<-biomaRt::useMart("ensembl", dataset= output_species_dataset) #set the biomart species to the new species
    ortho <- biomaRt::getBM(attributes= c(output_ids,"entrezgene_id"), filters= "entrezgene_id" ,values= id$species_two, mart=marts,uniqueRows=TRUE,bmHeader=FALSE) # run the homologous entrezids to the other gene identifiers
    non_na <- which(is.na(ortho$entrezgene_id) == FALSE) #pick non na values from ortholog conversion
    ortho <- ortho[non_na,] #subset only nonna values
    non_dups <- which(duplicated(ortho$entrezgene_id) == FALSE) #remove duplicates
    ortho <- ortho[non_dups,] #subset only duplicates
    ortho<-data.frame(ortho) # set the ortho variable as a dataframe
    names(ortho)[names(ortho) == "entrezgene_id"] <- "species_two" # rename the entrezid as species two
    merged_species_data<-merge.data.frame(data.frame(ortho), data.frame(merged_data), by = "species_two") # merge the ortholog dataframe and the original data
    names(merged_species_data)[names(merged_species_data) == "species_two"] <- paste(output_species,"Entrez_ID",sep="_") # clean up names
    names(merged_species_data)[names(merged_species_data) == "species_one"] <- paste(input_species,"Entrez_ID",sep="_") #clean up names
    colnames(merged_species_data)[which(names(merged_species_data) == "entrezgene_id")] <- paste(output_species,"Entrez_ID",sep="_") #clean up names
    colnames(merged_species_data)[which(names(merged_species_data) == "ensembl_gene_id")] <- paste(output_species,"Ensembl_ID",sep="_") #clean up names
    colnames(merged_species_data)[which(names(merged_species_data) == "hgnc_symbol")] <- paste(output_species,"HGNC_Symbol",sep="_") #clean up names
    colnames(merged_species_data)[which(names(merged_species_data) == "mgi_symbol")] <- paste(output_species,"MGI_Symbol",sep="_") #clean up names
    colnames(merged_species_data)[which(names(merged_species_data) == "external_gene_name")] <- paste(output_species,"Gene_Name",sep="_") #clean up names
    colnames(merged_species_data)[which(names(merged_species_data) == "ensembl_gene_id_version")] <- paste(output_species,"Ensembl_ID_with_Version_ID",sep="_") #clean up names
    merged_species_data <- merged_species_data[!duplicated(as.list(merged_species_data))] #remove duplicates
    return(merged_species_data) # return the final table
  }
}