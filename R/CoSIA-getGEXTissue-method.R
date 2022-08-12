#' getGExTissue Generic
#'
#' @param object
#'
#' @return
#' @export
#'
#' @examples

setGeneric("getGExTissue", function(object) standardGeneric("getGExTissue"))

#CoSIAn getGExTissue

setMethod("getGExTissue", signature(object = "CoSIAn"), function(object) {
  # user's input of the function
  if (object@map_tissues == "all_tissues") {
    gene_species <- object@map_species
    if(length(gene_species)> 1){
      stop("Error: map_species slot has more than 1 species. This method only works if the map_species slot has only 1 species identified.")
    }
    else{
      if (gene_species == "m_musculus") {
        bgee_species <- filter(Experimental_Hub_File, Species == "Mus_musculus")
        gene_set<-object@converted_id
        gene_list<- gene_set #call the specific column 
      }
      if (gene_species == "r_norvegicus") {
        bgee_species <- filter(Experimental_Hub_File, Species == "Rattus_norvegicus")
        gene_set<-object@converted_id
        gene_list<- gene_set
      }
      if (gene_species == "d_rerio") {
        bgee_species <- filter(Experimental_Hub_File, Species == "Danio_rerio")
        gene_set<-object@converted_id
        gene_list<- gene_set
      }
      if (gene_species == "h_sapiens") {
        bgee_species <- filter(Experimental_Hub_File, Species == "Homo_sapiens")
        gene_set<-object@converted_id
        gene_list<- gene_set
      }
    }
    #have to pull the conversion id here

    
    gene_specific_data <- dplyr::filter(bgee_species, Gene.ID == )
    sample_size <- data.frame(table(gene_specific_data$Anatomical.entity.name))
    colnames(sample_size)[which(names(sample_size) == "Var1")] <- "Anatomical.entity.name"
    values <- aggregate(data = gene_specific_data, x = gene_specific_data$TPM, by = list(gene_specific_data$Anatomical.entity.name), FUN = median)
    colnames(values)[which(names(values) == "Group.1")] <- "Anatomical.entity.name"
    value <- merge(values, sample_size, by = "Anatomical.entity.name")
    value$AEN <- paste(value$Anatomical.entity.name, "(n=", value$Freq, ")", sep = "")
    gene_specific_data <- merge(value, gene_specific_data, by = "Anatomical.entity.name")
    return(gene_specific_data)
    } 
  else {
    gene_species <- object@map_species
    #add an error here about gene species only being a length of 1
    if (gene_species == "mus_musculus") {
      bgee_species <- filter(Experimental_Hub_File, Species == "Mouse")
    }
    if (gene_species == "rattus_norvegicus") {
      bgee_species <- filter(Experimental_Hub_File, Species == "Rat")
    }
    if (gene_species == "danio_rerio") {
      bgee_species <- filter(Experimental_Hub_File, Species == "Zebrafish")
    }
    if (gene_species == "homo_sapiens") {
      bgee_species <- filter(Experimental_Hub_File, Species == "Human")
    }
    species_specific <- data.frame(dplyr::select(bgee_species, Gene.ID, Experiment.ID, Anatomical.entity.ID, Anatomical.entity.name, Read.count,
                                                 TPM, FPKM, Detection.flag))
    species_specific$Anatomical.entity.name <- gsub("\"", "", species_specific$Anatomical.entity.name)
    gene_specific_data <- dplyr::filter(species_specific, Gene.ID %in% object@single_gene)
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
})

