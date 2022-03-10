setMethod("getSpeciesExpression", signature(object = "CosiaExpressSpecies"), function(object) { # user's input of the function
  df<- data.frame("") #set up a blank dataframe to fill in using the filtered data that you are pulling
  x<-cbind(data.frame(object@list_of_ensembl_ids), data.frame(object@list_of_respective_species)) #combine the two arguments ensembl id list and the species that correlate with it into a data frame to parse through.
  for (val in 1:nrow(x)) { #takes about 10 mins to run this entire for loop if the data has already been downloaded to the location of your choice.
    single_gene<-x[val,1] #the for loop sets the ensembl id based on the number of rows of species data that is provided
    gene_species<-x[val,2] #the for loop sets the matching species to the ensembl id using the list that is provided by the user
    bgee_species <- BgeeDB::Bgee$new(species = gene_species, dataType= "rna_seq", pathToData = object@pathToData) # call Bgee data based on the species and download the data to the pathToData location
    bgee <- BgeeDB::getData(bgee_species) # get the data as a bgee object
    species_specific<- data.frame(dplyr::select(bgee, Gene.ID, Experiment.ID, Anatomical.entity.ID, Anatomical.entity.name, Read.count, TPM, FPKM,Detection.flag)) # filter the Bgee data by these columns
    species_specific$Anatomical.entity.name<-gsub('"',"",species_specific$Anatomical.entity.name) # remove the quotations from the anatomical entities
    gene_specific_data<-dplyr::filter(species_specific,Gene.ID== single_gene) # filter by the gene that is given
    tissue_specific_data<-dplyr::filter(gene_specific_data,Anatomical.entity.name %in% object@single_tissue) #filter by the specific anatomical entity name in the single tissue
    tissue_specific_data$Species <- gene_species #set a column with the species name so that can be sorted into columns
    df <-merge(x=df,y=tissue_specific_data,all=TRUE) #merge the bgee data for the species into the blank data frame so that it can be collected outside of the for loop
  }
  df <- subset( df, select = -c(X..)) # remove the X column that was blank and being used to start the blank data frame
  
  #start of the mediana nd sample size calculations
  sample_size <- data.frame(table(df$Species)) #find the sample sizes for the different species in order to put that in the legend
  colnames(sample_size)[which(names(sample_size) == "Var1")] <- "Species" # change name from Var 1 to species
  values<-aggregate(data = df,x = df$TPM,by = list(df$Species),FUN = median) # find the median TPM for each species
  colnames(values)[which(names(values) == "Group.1")] <- "Species" # change name
  value<- merge(values,sample_size,by="Species") #merge the median and sample size into the data
  value$S <- paste(value$Species,"(n=",value$Freq,")", sep = "") #set a new column with the new species column with sample size
  df_1 <- merge(value,df,by="Species") # merge the calculation to the main df
  object@dataframe<-df_1
    brewer.pal.info<-RColorBrewer::brewer.pal.info #set brewer's color
    palette5 <- RColorBrewer::brewer.pal.info[brewer.pal.info$category == "qual",] # set the palette
    color <- unlist(mapply(RColorBrewer::brewer.pal, palette5$maxcolors,rownames(palette5))) # set a list of unique colors to use
    #plotly figure making
    fig <-plotly::plot_ly(df_1, x = ~S, y = ~TPM, type = 'scatter', mode = 'markers', color = ~S, colors = color) %>%
      plotly::add_markers(x=~ S,y=~x ,marker = list(symbol = 'line-ew' ,size = 10,line = list(color = "grey",width = 2)))
    fig<-fig %>% plotly::add_trace(marker = list(size = 8,line = list(color = 'black',width = .75)),
                                   showlegend = F)
    fig <- fig %>%
      plotly::layout(title = stringr::str_wrap(paste("Gene Expression of the gene", object@list_of_ensembl_ids[1], "in", object@single_tissue, sep=" ")),
                     xaxis = list(title = "Anatomical Entity Name",size = 2),
                     yaxis = list(title = "TPM (transcript per million)",zeroline = F),
                     showlegend=FALSE)
    
    return(fig) # return the figure
}
)
