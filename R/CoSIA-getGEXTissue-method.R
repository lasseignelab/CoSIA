#' getGEx Generic
#'
#' @param object
#'
#' @return
#' @export
#'
#' @examples

setGeneric("getGEx", function(object) standardGeneric("getGEx"))
#CoSIAn getGEx
setMethod("getGEx", signature(object = "CoSIAn"), function(object) {
  id_dataframe<-object@converted_id
  id_dataframe<- as.data.frame(id_dataframe)
  map_species<- object@map_species
  map_tissues<- object@map_tissues
  Species_SWITCH <- Vectorize(vectorize.args = "map_species", FUN = function(map_species) {
    switch(as.character(map_species), 
           h_sapiens ="h_sapiens_ensembl_id",
           m_musculus = "m_musculus_ensembl_id",
           r_norvegicus = "r_norvegicus_ensembl_id",
           d_rerio = "d_rerio_ensembl_id",
           d_melanogaster = "d_melanogaster_ensembl_id",
           c_elegans = "c_elegans_ensembl_id",
           stop("Error: Invalid map_species in CoSIAn Object. Make sure the species in the map_species slot are an avalible model 
                                    organism and are in the correct format.")
    )
  })
  map_species <- Species_SWITCH(map_species)
  id_dataframe<- dplyr::select(id_dataframe,matches(map_species))
  #load EH_Data here (get it off of cheaha soon the EH method will be used to pull the data here)
  return_filtered_Gex_data<- function(map_species){
    GEx_data<- data.frame(matrix(ncol = 12, nrow = 0))
    colnames(GEx_data) <- c('Anatomical_entity_name', 'Ensembl_ID', 'Sample_size', 'Minimum_TPM', 'First_Quartile_TPM', 'Median_TPM', 'Third_Quartile_TPM','Maximum_TPM','Standard_Deviation','Experiment_ID','Anatomical_entity_ID','Species')
    if (any(map_species == "m_musculus_ensembl_id")) {
      bgee_species <- dplyr::filter(Experimental_Hub_File, Species == "Mus_musculus")
      gene_specific_data <- dplyr::filter(bgee_species, Ensembl_ID == id_dataframe$m_musculus)
      GEx_data <- merge(GEx_data, gene_specific_data, by = "Ensembl_ID")
    }
    if (any(map_species == "r_norvegicus_ensembl_id")) {
      bgee_species <- dplyr::filter(Experimental_Hub_File, Species == "Rattus_norvegicus")
      gene_specific_data <- dplyr::filter(bgee_species, Ensembl_ID == id_dataframe$r_norvegicus)
      GEx_data <- merge(GEx_data, gene_specific_data, by = "Ensembl_ID")
    }
    if (any(map_species == "d_rerio_ensembl_id")) {
      bgee_species <- dplyr::filter(Experimental_Hub_File, Species == "Danio_rerio")
      gene_specific_data <- dplyr::filter(bgee_species, Ensembl_ID == id_dataframe$d_rerio)
      GEx_data <- merge(GEx_data, gene_specific_data, by = "Ensembl_ID")
    }
    if (any(map_species == "h_sapiens_ensembl_id")) {
      bgee_species <- dplyr::filter(Experimental_Hub_File, Species == "Homo_sapiens")
      gene_specific_data <- dplyr::filter(bgee_species, Ensembl_ID == id_dataframe$h_sapiens)
      GEx_data <- merge(GEx_data, gene_specific_data, by = "Ensembl_ID")
    }
    if (any(map_species == "c_elegans_ensembl_id")) {
      bgee_species <- dplyr::filter(Experimental_Hub_File, Species == "Caenorhabditis_elegans")
      gene_specific_data <- dplyr::filter(bgee_species, Ensembl_ID == id_dataframe$c_elegans)
      GEx_data <- merge(GEx_data, gene_specific_data, by = "Ensembl_ID")
    }
    if (any(map_species == "d_melanogaster_ensembl_id")) {
      bgee_species <- dplyr::filter(Experimental_Hub_File, Species == "Drosophila_melanogaster")
      gene_specific_data <- dplyr::filter(bgee_species, Ensembl_ID == id_dataframe$d_melanogaster)
      GEx_data <- merge(GEx_data, gene_specific_data, by = "Ensembl_ID")
    }
    return(GEx_data)
  }
  GEx_data<-sapply(map_species,return_filtered_Gex_data)
  combine_GEx<-data.frame(0)
  for (i in 1: length(GEx_data)){
    combine_GEx<-merge(combine_GEx,GEx_data[i])
  }
  return(combine_GEx)
  
  
  
  
  tissue_specific_data <- dplyr::filter(gene_specific_data, Anatomical.entity.name %in% object@tissues)
  sample_size <- data.frame(table(tissue_specific_data$Anatomical.entity.name))
  colnames(sample_size)[which(names(sample_size) == "Var1")] <- "Anatomical.entity.name"
  values <- aggregate(data = tissue_specific_data, x = tissue_specific_data$TPM, by = list(tissue_specific_data$Anatomical.entity.name),
                      FUN = median)
  colnames(values)[which(names(values) == "Group.1")] <- "Anatomical.entity.name"
  value <- merge(values, sample_size, by = "Anatomical.entity.name")
  value$AEN <- paste(value$Anatomical.entity.name, "(n=", value$Freq, ")", sep = "")
  tissue_specific_data <- merge(value, tissue_specific_data, by = "Anatomical.entity.name")
  return(tissue_specific_data)
}
)