AnnotateDBI<-function(input_id,input_dataset,output_ids,input_species,output_species,species_number,input_org, output_org, ortholog_database){
  #AnnotationDBI
  input_id <- switch(input_id,
                     "Entrez.id"="ENTREZID",
                     "Ensembl.id"="ENSEMBL",
                     "Ensembl.id.version"="ENSEMBLIDVERSION",
                     "Gene.name"="GENENAME",
                     "Symbol"="SYMBOL")
  for(x in length(output_ids)){
    output_ids[x] <- switch(output_ids[x],
                            "Entrez.id"="ENTREZID",
                            "Ensembl.id"="ENSEMBL",
                            "Ensembl.id.version"="ENSEMBLIDVERSION",
                            "Gene.name"="GENENAME",
                            "Symbol"="SYMBOL")
  }
  if (output_species==input_species){ #code follows this path if the user chooses the same input species as their output species
    if (input_id=="ENSEMBLIDVERSION"){ #code follows this path if the user chooses ENSEMBLIDVERSION as their input id
      input_dataset_new<-CoSIA::remove_version_number(input_dataset,input_species) #puts the genes through the remove version number function to remove the version number (everything after the dot)
      #return(input_dataset_new)
      annotationDbi_data<- AnnotationDbi::select(input_org, keys=as.character(as.matrix(input_dataset_new$ENSEMBL)),columns=output_ids , keytype="ENSEMBL") # after the version number has been removed run the gene through the annotationDBI package using ensembl id as the gene identifer
      merge<-merge.data.frame(data.frame(annotationDbi_data), data.frame(input_dataset_new), by = "ENSEMBL") # merge the input and output data together
      return(merge) #works
    }
    else{ #code follows this path if the user chooses an ID other than ENSEMBLIDVERSION as their input id
      annotationDbi_data<- AnnotationDbi::select(input_org, keys= input_dataset, columns= output_ids , keytype= input_id) #run the code through annotationDBI
      return(annotationDbi_data) #works
    }
  }
  else{ # this code follows this path if the input gene is different than the output gene
    if (input_id=="ENSEMBLIDVERSION"){ #code follows this path if the user chooses ENSEMBLIDVERSION as their input id
      input_dataset_new<-CoSIA::remove_version_number(input_dataset,input_species) #puts the genes through the remove version number function to remove the version number (everything after the dot)
      input_dataset_new<- data.frame(input_dataset_new) #sets the new list of genes without versionIDs as a dataframe
      output_data<- AnnotationDbi::select(input_species, keys=as.character(as.matrix(input_dataset_new)), columns="ENTREZID","ENSEMBL", keytype="ENSEMBL") #runs the gene through annotationdbi to get EntrezID conversion in order to conduct ortholog mapping
      output_data<-na.omit(output_data) #omits the missing values that were not mapped
      merged_data<-merge.data.frame(data.frame(output_data), data.frame(input_dataset_new), by = "ENSEMBL") #merges the original df with the entrezid mapped data set
      id <- CoSIA::homolog(output_data$ENTREZID, species_number, ortholog_database) # run the EntrezID found in the new dataset through the homolog function
      id<-data.frame(id) #makes the homolog output into a dataframe
      names(merged_data)[names(merged_data) == "ENTREZID"] <- "species_one" #changes the name of the entrezid column to mathc the id dataframe as species one
      merged_data<-merge.data.frame(data.frame(merged_data), data.frame(id), by = "species_one") # merges the merged data and the homolog values that were found
      colnames(merged_data)[which(names(merged_data) == "ENTREZID")] <- paste(input_species,"Entrez_ID",sep="_") #change the names to more formal names for the columns
      colnames(merged_data)[which(names(merged_data) == "ENSEMBL")] <- paste(input_species,"Ensembl_ID",sep="_") #change the names to more formal names for the columns
      colnames(merged_data)[which(names(merged_data) == "SYMBOL")] <- paste(input_species,"Symbol",sep="_") #change the names to more formal names for the columns
      colnames(merged_data)[which(names(merged_data) == "ENSEMBLIDVERSION")] <- paste(input_species,"Ensembl_ID_with_Version_ID",sep="_") #change the names to more formal names for the columns
      #names(output_data)[names(output_data) == "ENTREZID"] <- "species_one" #change the name of the output dataframe to species_one for merging later
      ortho<- AnnotationDbi::select(output_org, keys=as.character(as.matrix(merged_data$species_two)), columns=c(output_ids,"ENTREZID"), keytype="ENTREZID") # next we are going to take the entrezids of the orthologs and return the output values of the gene
      non_na <- which(is.na(ortho$ENTREZID) == FALSE) # determine the indices for the non-NA genes
      ortho <- ortho[non_na,] #return only the genes with annotations using indices
      non_dups <- which(duplicated(ortho$ENTREZID) == FALSE) #determine the indices for the non-duplicated genes
      ortho<- ortho[non_dups,]#return only the non-duplicated genes using indices
      ortho<-data.frame(ortho) #make this into a dataframe
      names(ortho)[names(ortho) == "ENTREZID"] <- "species_two" # set the entrezID into a dataframe
      merged_data<-merge.data.frame(data.frame(ortho), data.frame(merged_data), by = "species_two") #merge the ortholog conversion dataframe with the merge dataframe
      #names(merged_data)[names(merged_data) == "species_two"] <- paste(output_species,"Entrez_ID",sep="_") #changes name to more formal names
      #names(merged_data)[names(merged_mouse_data) == "species_one"] <- paste(input_species,"Entrez_ID",sep="_") #changes name to more formal names
      colnames(merged_data)[which(names(merged_data) == "ENTREZID")] <- paste(output_species,"Entrez_ID",sep="_") #changes name to more formal names
      colnames(merged_data)[which(names(merged_data) == "ENSEMBL")] <- paste(output_species,"Ensembl_ID",sep="_") #changes name to more formal names
      colnames(merged_data)[which(names(merged_data) == "SYMBOL")] <- paste(output_species,"Symbol",sep="_") #changes name to more formal names
      colnames(merged_data)[which(names(merged_data) == "ENSEMBLIDVERSION")] <- paste(output_species,"Ensembl_ID_with_Version_ID",sep="_") #changes name to more formal names
      merged_data <- merged_data[!duplicated(as.list(merged_data))] # removes duplicates from the data
      merged_data = subset(merged_data, select = -c(species_two,species_one)) # removes the species two and species one parts of the data
      return(merged_data) #works
    }
    else{ #code follows this path if the user chooses an ID other than ENSEMBLIDVERSION as their input id
      output_data<- AnnotationDbi::select(input_org, keys=as.character(as.matrix(input_dataset)),columns="ENTREZID",input_id, keytype=input_id) #runs the gene through annotationdbi to get EntrezID conversion in order to conduct ortholog mapping
      output_data<-na.omit(output_data) # omits the NA values
      id <- CoSIA::homolog(output_data$ENTREZID, species_number, ortholog_database) #send the entrezid ids through the homolog conversion
      id<-data.frame(id) # make these ortholog ids into a new dataframe
      names(output_data)[names(output_data) == "ENTREZID"] <- "species_one" # make the EntrezID that were used for the conversion into species one
      merged_data<-merge.data.frame(data.frame(output_data), data.frame(id), by = "species_one") # merge the two dataframes with the original ids and the ortholog ids together
      colnames(merged_data)[which(names(merged_data) == "ENTREZID")] <- paste(input_species,"Entrez_ID",sep="_") #changes name to more formal names
      colnames(merged_data)[which(names(merged_data) == "ENSEMBL")] <- paste(input_species,"Ensembl_ID",sep="_")#changes name to more formal names
      colnames(merged_data)[which(names(merged_data) == "SYMBOL")] <- paste(input_species,"Symbol",sep="_")#changes name to more formal names
      colnames(merged_data)[which(names(merged_data) == "ENSEMBLIDVERSION")] <- paste(input_species,"Ensembl_ID_with_Version_ID",sep="_")#changes name to more formal names
      ortho <- AnnotationDbi::select(output_org, keys=as.character(as.matrix(merged_data$species_two)),columns=output_ids,"ENTREZID", keytype="ENTREZID") # now that you have the entrezids of the converted species we can convert to different gene identifiers
      non_na <- which(is.na(ortho$ENTREZID) == FALSE) # Determine the indices for the non-NA genes
      ortho <- ortho[non_na,]  # Return only the genes with annotations using indices
      non_dups <- which(duplicated(ortho$ENTREZID) == FALSE) # Determine the indices for the non-duplicated genes
      ortho <- ortho[non_dups,] # Return only the non-duplicated genes using indices
      ortho<-data.frame(ortho) # make the ortholog conversion into a dataframe
      names(ortho)[names(ortho) == "ENTREZID"] <- "species_two" #rename the species into species two in order to merge the orthologs and the original data
      merged_data<-merge.data.frame(data.frame(ortho), data.frame(merged_data), by = "species_two") # merge
      colnames(merged_data)[names(merged_data) == "species_two"] <- paste(output_species,"Entrez_ID",sep="_") #changes name to more formal names
      #names(merged_data)[names(merged_data) == "species_one"] <- paste(input_species,"Entrez_ID",sep="_") #changes name to more formal names
      colnames(merged_data)[which(names(merged_data) == "ENTREZID")] <- paste(output_species,"Entrez_ID",sep="_") #changes name to more formal names
      colnames(merged_data)[which(names(merged_data) == "ENSEMBL")] <- paste(output_species,"Ensembl_ID",sep="_")#changes name to more formal names
      colnames(merged_data)[which(names(merged_data) == "SYMBOL")] <- paste(output_species,"Symbol",sep="_") #changes name to more formal names
      colnames(merged_data)[which(names(merged_data) == "ENSEMBLIDVERSION")] <- paste(output_species,"Ensembl_ID_with_Version_ID",sep="_") #changes name to more formal names
      merged_data <- merged_data[!duplicated(as.list(merged_data))] # remove duplicates
      merged_data = subset(merged_data, select = -c(species_one)) #remove species one
      return(merged_data)#works
    }
  }
}