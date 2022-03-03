remove_version_number<- function(input_dataset, species){
  if (species=="mus_musculus"){ # if the user's input species is a mouse the code follows this path
    input_dataset<-data.frame(input_dataset) # create a new input dataset to manipulate
    colnames(input_dataset)[1] <- "ENSEMBLIDVERSION" # set the column name as ensemblid
    input_dataset = cbind(input_dataset, replicate(1,input_dataset$ENSEMBLIDVERSION)) # create a new identical column that will be used to maniputlate
    ensembl_id<- input_dataset  # set the input dataset into a new column
    ensembl_id<-data.frame(ensembl_id)
    x<-1:(nrow(ensembl_id))
    for (i in seq_along(x)) {
      char<-input_dataset[i,1]
      char1<-substr(char, start = 1, stop = 15)
      input_dataset[i,1]<-char1
    }
    input_dataset<-data.frame(input_dataset)
    colnames(input_dataset)[which(names(input_dataset) == "ENSEMBLIDVERSION")] <- "ENSEMBL"
    colnames(input_dataset)[which(names(input_dataset) == "replicate.1..input_dataset.ENSEMBLIDVERSION.")] <- "ENSEMBLIDVERSION"
    return (input_dataset)
  }
  
  if (species == "homo_sapiens"){
    input_dataset<-data.frame(input_dataset)
    colnames(input_dataset)[1] <- "ENSEMBLIDVERSION"
    input_dataset = cbind(input_dataset, replicate(1,input_dataset$ENSEMBLIDVERSION))
    ensembl_id<- input_dataset
    ensembl_id<-data.frame(ensembl_id)
    x<-1:(nrow(ensembl_id))
    for (i in seq_along(x)) {
      char<-input_dataset[i,1]
      char1<-substr(char, start = 1, stop = 18)
      input_dataset[i,1]<-char1
    }
    input_dataset<-data.frame(input_dataset)
    colnames(input_dataset)[which(names(input_dataset) == "ENSEMBLIDVERSION")] <- "ENSEMBL"
    colnames(input_dataset)[which(names(input_dataset) == "replicate.1..input_dataset.ENSEMBLIDVERSION.")] <- "ENSEMBLIDVERSION"
    return (input_dataset)
  }
  
  if (species=="danio_rerio"){
    input_dataset<-data.frame(input_dataset)
    colnames(input_dataset)[1] <- "ENSEMBLIDVERSION"
    input_dataset = cbind(input_dataset, replicate(1,input_dataset$ENSEMBLIDVERSION))
    ensembl_id<- input_dataset
    ensembl_id<-data.frame(ensembl_id)
    x<-1:(nrow(ensembl_id))
    for (i in seq_along(x)) {
      char<-input_dataset[i,1]
      char1<-substr(char, start = 1, stop = 18)
      input_dataset[i,1]<-char1
    }
    input_dataset<-data.frame(input_dataset)
    colnames(input_dataset)[which(names(input_dataset) == "ENSEMBLIDVERSION")] <- "ENSEMBL"
    colnames(input_dataset)[which(names(input_dataset) == "replicate.1..input_dataset.ENSEMBLIDVERSION.")] <- "ENSEMBLIDVERSION"
    return (input_dataset)
  }
  
  if (species == "c_elegans"){
    stop("Error. Incorrect ID evaluation. C. elegans do not have Ensembl Ids with Version.")
  }
  
  else{
    stop("Error. Invalid species. Make sure it matches the proper format.")
  }
  
}