#' getGExMetrics Generic
#'
#' @param object
#'
#' @export
setGeneric("getGExMetrics", function(object) standardGeneric("getGExMetrics"))
#CoSIAn getGEx


#' getGExMetrics Method
#'
#' @param object CoSIAn. 
#'
#' @return CoSIAn Object with metric slot filled
#' @export
#' @importFrom dplyr group_by summarise
#' @importFrom tibble remove_rownames column_to_rownames
#'
#' @examples
#' Kidney_gene_metric<-getGExMetrics(Kidney_gene_gex)

setMethod("getGExMetrics", signature(object = "CoSIAn"), function(object) {
  id_dataframe<-object@converted_id
  id_dataframe<- as.data.frame(id_dataframe)
  map_species<- object@map_species
  map_species<- as.character(map_species)
  Converted_Species_SWITCH <- Vectorize(vectorize.args = "map_species", FUN = function(map_species) {
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
  map_species <- Converted_Species_SWITCH(map_species)
  id_dataframe<- dplyr::select(id_dataframe,tidyselect::matches(map_species))

  map_tissues<- object@map_tissues
  map_tissues<- as.character(map_tissues)
  Species_SWITCH <- Vectorize(vectorize.args = "map_species", FUN = function(map_species) {
    switch(as.character(map_species), 
           h_sapiens_ensembl_id ="Homo_sapiens",
           m_musculus_ensembl_id = "Mus_musculus",
           r_norvegicus_ensembl_id = "Rattus_norvegicus",
           d_rerio_ensembl_id = "Danio_rerio",
           d_melanogaster_ensembl_id = "Drosophila_melanogaster",
           c_elegans_ensembl_id = "Caenorhabditis_elegans",
           stop("Error: Invalid map_species in CoSIAn Object. Make sure the species in the map_species slot are an avalible model 
                                    organism and are in the correct format.")
    )
  })
  map_species <- Species_SWITCH(map_species)
  metric_type<- object@metric_type
  metric_type<- as.character(metric_type)
  
  CV_function <- function(x, na.rm = FALSE){
    stopifnot(is.numeric(x))
    if(na.rm) x <- x[!is.na(x)]
    if(any(x < 0)) #changed this from <= to < 
      stop("Your data must be greater than zero!")
    stats::sd(x, na.rm = FALSE)/stats::median(x)
  }
  
  
  CV_Tissue <- function(map_species, map_tissues){
    filter_species <- dplyr::filter(Experimental_Hub_File,Species %in% map_species)
    filter_tissue <- dplyr::filter(filter_species,Anatomical_entity_name %in% map_tissues)
    id<-as.vector(t(id_dataframe))
    filter_gene <- dplyr::filter(filter_tissue,Ensembl_ID %in% id)
    filter_gex<-tidyr::separate_rows(filter_gene, VST)
    filter_gex$VST <- as.numeric(filter_gex$VST)
    CV_Tissue<-filter_gex %>% group_by(Ensembl_ID,Anatomical_entity_name,Species) %>% summarise(CV_Tissue = CV_function(VST, na.rm=FALSE))
    return(CV_Tissue)
  }
  
  CV_Species <- function(map_species, map_tissues){
    filter_species <- dplyr::filter(Experimental_Hub_File,Species %in% map_species)
    id<-as.vector(t(id_dataframe))
    filter_gene <- dplyr::filter(filter_species,Ensembl_ID %in% id)
    filter_gex<-tidyr::separate_rows(filter_gene, VST)
    filter_gex$VST <- as.numeric(filter_gex$VST)
    CV_Species<-filter_gex %>% group_by(Ensembl_ID,Species) %>% summarise(CV_Species = CV_function(VST, na.rm=FALSE))
    CV_Tissue<-CV_Tissue %>% tidyr::spread(Species)
    return(CV_Species)
  }
  
  #DS_Gene : output genes restricted by mapped tissues and gene set 
  DS_Gene<- function(map_species, map_tissues){
    filter_species <- dplyr::filter(Experimental_Hub_File,Species %in% map_species)
    filter_tissue <- dplyr::filter(filter_species,Anatomical_entity_name %in% map_tissues)
    id<-as.vector(t(id_dataframe))
    filter_gene <- dplyr::filter(filter_tissue,Ensembl_ID %in% id)
    filter_gene$Scaled_Median_VST <- as.numeric(filter_gene$Scaled_Median_VST)
    filter_gex<- dplyr::select(filter_gene, Anatomical_entity_name, Scaled_Median_VST, Ensembl_ID)
    filter_gex_D<- filter_gex%>% tidyr::pivot_wider(names_from = Ensembl_ID, values_from = Scaled_Median_VST)
    filter_gex_D <- filter_gex_D %>% remove_rownames %>% tibble::column_to_rownames(var="Anatomical_entity_name")
    filter_gex_D<- data.matrix(filter_gex_D, )
    ENTROPY_DIVERSITY_G<-data.frame(BioQC::entropyDiversity(filter_gex_D,norm = TRUE)) # across genes
    colnames(ENTROPY_DIVERSITY_G)[which(names(ENTROPY_DIVERSITY_G) == "BioQC..entropyDiversity.filter_gex_D..norm...TRUE.")] <- "Diversity"
    
    filter_gex<- data.frame(filter_gex)
    filter_gex_S<- filter_gex%>% tidyr::pivot_wider(names_from = Anatomical_entity_name, values_from = Scaled_Median_VST)
    filter_gex_S <- filter_gex_S %>% remove_rownames %>% column_to_rownames(var="Ensembl_ID")
    filter_gex_S<- data.matrix(filter_gex_S, )
    ENTROPY_SPECIFITY_G<-data.frame(BioQC::entropySpecificity(filter_gex_S,norm = TRUE)) # across tissues
    colnames(ENTROPY_SPECIFITY_G)[which(names(ENTROPY_SPECIFITY_G) == "BioQC..entropySpecificity.filter_gex_S..norm...TRUE.")] <- "Specificity"
    
    DS <- merge(ENTROPY_SPECIFITY_G, ENTROPY_DIVERSITY_G, by = 'row.names')
    colnames(DS)[which(names(DS) == "Row.names")] <- "Ensembl_ID"
    DS$Ensembl_ID <- as.character(DS$Ensembl_ID)
    Species<-dplyr::select(filter_gene, Species, Ensembl_ID)
    DS <- merge(DS, Species, by = 'Ensembl_ID')
    DS<-data.frame(unique(DS))
    rownames(DS) <- NULL
    return(DS)
  }
  
  #DS_Tissues: output is tissues restricted to mapped tissues and gene set
  DS_Tissue<- function(map_species, map_tissues){
    DS<-data.frame(matrix(ncol = 4, nrow = 0))
    colnames(DS)[which(names(DS) == "map_tissues")] <- "Anatomical_entity_name"
    for (x in 1:length(map_species)){
      filter_species <- dplyr::filter(Experimental_Hub_File,Species == map_species[x])
      filter_tissue <- dplyr::filter(filter_species,Anatomical_entity_name %in% map_tissues)
      id<-as.vector(t(id_dataframe))
      filter_gene <- dplyr::filter(filter_tissue,Ensembl_ID %in% id)
      filter_gene$Scaled_Median_VST <- as.numeric(filter_gene$Scaled_Median_VST)
      filter_gex<- dplyr::select(filter_gene, Anatomical_entity_name, Scaled_Median_VST, Ensembl_ID)
      
      filter_gex_D<- filter_gex%>% tidyr::pivot_wider(names_from = Anatomical_entity_name, values_from = Scaled_Median_VST)
      filter_gex_D <- filter_gex_D %>% remove_rownames %>% tibble::column_to_rownames(var="Ensembl_ID")
      filter_gex_D<- data.matrix(filter_gex_D, )
      ENTROPY_DIVERSITY_T<-data.frame(BioQC::entropyDiversity(filter_gex_D,norm = TRUE)) # across genes
      colnames(ENTROPY_DIVERSITY_T)[which(names(ENTROPY_DIVERSITY_T) == "BioQC..entropyDiversity.filter_gex_D..norm...TRUE.")] <- "Diversity"
      
      filter_gex<- data.frame(filter_gex)
      filter_gex_S<- filter_gex%>% tidyr::pivot_wider(names_from = Ensembl_ID, values_from = Scaled_Median_VST)
      filter_gex_S <- filter_gex_S %>% remove_rownames %>% column_to_rownames(var="Anatomical_entity_name")
      filter_gex_S<- data.matrix(filter_gex_S, )
      ENTROPY_SPECIFITY_T<-data.frame(BioQC::entropySpecificity(filter_gex_S,norm = TRUE)) # across tissues
      colnames(ENTROPY_SPECIFITY_T)[which(names(ENTROPY_SPECIFITY_T) == "BioQC..entropySpecificity.filter_gex_S..norm...TRUE.")] <- "Specificity"
      
      SDS <- merge(ENTROPY_SPECIFITY_T, ENTROPY_DIVERSITY_T, by = 'row.names')
      colnames(SDS)[which(names(SDS) == "Row.names")] <- "Anatomical_entity_name"
      SDS$Anatomical_entity_name <- as.character(SDS$Anatomical_entity_name)
      Species<-dplyr::select(filter_gene, Species, Anatomical_entity_name)
      SDS <- merge(SDS, Species, by = 'Anatomical_entity_name')
      DS <- rbind(DS,SDS)
    }
    DS<-data.frame(unique(DS))
    rownames(DS) <- NULL
    return(DS)
    
  }
  #DS_Gene_all: outputs genes only restricted to selected genes across all tissues
  DS_Gene_all<- function(map_species, map_tissues){
    DS<-data.frame(matrix(ncol = 4, nrow = 0))
    colnames(DS)[1] <- "Ensembl_ID"
    for (x in 1:length(map_species)){
      filter_species <- dplyr::filter(Experimental_Hub_File,Species == map_species[x])
      id<-as.vector(t(id_dataframe))
      filter_gene <- dplyr::filter(filter_species,Ensembl_ID %in% id)
      filter_gene$Scaled_Median_VST <- as.numeric(filter_gene$Scaled_Median_VST)
      filter_gex<- dplyr::select(filter_gene, Anatomical_entity_name, Scaled_Median_VST, Ensembl_ID)
      
      filter_gex_D<- filter_gex%>% tidyr::pivot_wider(names_from = Ensembl_ID, values_from = Scaled_Median_VST)
      filter_gex_D <- filter_gex_D %>% remove_rownames %>% tibble::column_to_rownames(var="Anatomical_entity_name")
      filter_gex_D<- data.matrix(filter_gex_D, )
      ENTROPY_DIVERSITY_G<-data.frame(BioQC::entropyDiversity(filter_gex_D,norm = TRUE)) # across genes
      colnames(ENTROPY_DIVERSITY_G)[which(names(ENTROPY_DIVERSITY_G) == "BioQC..entropyDiversity.filter_gex_D..norm...TRUE.")] <- "Diversity"
      
      filter_gex<- data.frame(filter_gex)
      filter_gex_S<- filter_gex%>% tidyr::pivot_wider(names_from = Anatomical_entity_name, values_from = Scaled_Median_VST)
      filter_gex_S <- filter_gex_S %>% remove_rownames %>% column_to_rownames(var="Ensembl_ID")
      filter_gex_S<- data.matrix(filter_gex_S, )
      ENTROPY_SPECIFITY_G<-data.frame(BioQC::entropySpecificity(filter_gex_S,norm = TRUE)) # across tissues
      colnames(ENTROPY_SPECIFITY_G)[which(names(ENTROPY_SPECIFITY_G) == "BioQC..entropySpecificity.filter_gex_S..norm...TRUE.")] <- "Specificity"
      
      SDS <- merge(ENTROPY_SPECIFITY_G, ENTROPY_DIVERSITY_G, by = 'row.names')
      colnames(SDS)[which(names(SDS) == "Row.names")] <- "Ensembl_ID"
      SDS$Ensembl_ID <- as.character(SDS$Ensembl_ID)
      Species<-dplyr::select(filter_gene, Species, Ensembl_ID)
      SDS <- merge(SDS, Species, by = 'Ensembl_ID')
      DS <- rbind(DS,SDS)
    }
    DS<-data.frame(unique(DS))
    rownames(DS) <- NULL
    return(DS)
  }
  
  #DS_Tissues_All: output is tissues restricted to mapped tissues across all genes
  DS_Tissue_all<- function(map_species, map_tissues){
    DS<-data.frame(matrix(ncol = 4, nrow = 0))
    colnames(DS)[which(names(DS) == "map_tissues")] <- "Anatomical_entity_name"
    for (x in 1:length(map_species)){
      filter_species <- dplyr::filter(Experimental_Hub_File,Species == map_species[x])
      filter_tissue <- dplyr::filter(filter_species,Anatomical_entity_name %in% map_tissues)
      filter_tissue$Scaled_Median_VST <- as.numeric(filter_tissue$Scaled_Median_VST)
      filter_gex<- dplyr::select(filter_tissue, Anatomical_entity_name, Scaled_Median_VST, Ensembl_ID)
      
      filter_gex_D<- filter_gex%>% tidyr::pivot_wider(names_from = Anatomical_entity_name, values_from = Scaled_Median_VST)
      filter_gex_D <- filter_gex_D %>% remove_rownames %>% tibble::column_to_rownames(var="Ensembl_ID")
      filter_gex_D<- data.matrix(filter_gex_D, )
      ENTROPY_DIVERSITY_T<-data.frame(BioQC::entropyDiversity(filter_gex_D,norm = TRUE)) # across genes
      colnames(ENTROPY_DIVERSITY_T)[which(names(ENTROPY_DIVERSITY_T) == "BioQC..entropyDiversity.filter_gex_D..norm...TRUE.")] <- "Diversity"
      
      filter_gex<- data.frame(filter_gex)
      filter_gex_S<- filter_gex%>% tidyr::pivot_wider(names_from = Ensembl_ID, values_from = Scaled_Median_VST)
      filter_gex_S <- filter_gex_S %>% remove_rownames %>% column_to_rownames(var="Anatomical_entity_name")
      filter_gex_S<- data.matrix(filter_gex_S, )
      ENTROPY_SPECIFITY_T<-data.frame(BioQC::entropySpecificity(filter_gex_S,norm = TRUE)) # across tissues
      colnames(ENTROPY_SPECIFITY_T)[which(names(ENTROPY_SPECIFITY_T) == "BioQC..entropySpecificity.filter_gex_S..norm...TRUE.")] <- "Specificity"
      
      SDS <- merge(ENTROPY_SPECIFITY_T, ENTROPY_DIVERSITY_T, by = 'row.names')
      colnames(SDS)[which(names(SDS) == "Row.names")] <- "Anatomical_entity_name"
      SDS$Anatomical_entity_name <- as.character(SDS$Anatomical_entity_name)
      Species<-dplyr::select(filter_tissue, Species, Anatomical_entity_name)
      SDS <- merge(SDS, Species, by = 'Anatomical_entity_name')
      DS <- rbind(DS,SDS)
    }
    DS<-data.frame(unique(DS))
    rownames(DS) <- NULL
    return(DS)
  }
  
  
  if(metric_type == "CV_Tissue"){
    CV_Tissue<-CV_Tissue(map_species,map_tissues)
    object@metric<-CV_Tissue
  }
  else if(metric_type == "CV_Species"){
    CV_Species<-CV_Species(map_species,map_tissues)
    object@metric<-CV_Species
  }
  else if(metric_type == "DS_Gene"){
    DS_Gene<-DS_Gene(map_species,map_tissues)
    object@metric<-DS_Gene
  }
  else if(metric_type == "DS_Tissue"){
    DS_Tissue<-DS_Tissue(map_species,map_tissues)
    object@metric<-DS_Tissue
  }
  else if(metric_type == "DS_Gene_all"){
    DS_Gene_all<-DS_Gene_all(map_species,map_tissues)
    object@metric<-DS_Gene_all
  }
  else if(metric_type == "DS_Tissue_all"){
    DS_Tissue_all<-DS_Tissue_all(map_species,map_tissues)
    object@metric<-DS_Tissue_all
  }
  else(
    stop("Error: invalid metric type")
  )
  
return(object)
}) 
    
    

  
  