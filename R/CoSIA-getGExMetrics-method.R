#CoSIA getGExMetrics Function
#' getGExMetrics Generic
#'
#' @param object
#'
#' @return
#' @export
#'
#' @examples

setGeneric("getGExMetrics", function(object) standardGeneric("getGExMetrics"))
#CoSIAn getGEx
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
  id_dataframe<- dplyr::select(id_dataframe,matches(map_species))

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
    sd(x, na.rm = FALSE)/median(x)
  }
  
  CV_Tissue <- function(map_species, map_tissues){
    filter_species <- dplyr::filter(Experimental_Hub_File,Species %in% map_species)
    filter_tissue <- dplyr::filter(filter_species,Anatomical_entity_name %in% map_tissues)
    id<-as.vector(t(id_dataframe))
    filter_gene <- dplyr::filter(filter_tissue,Ensembl_ID %in% id)
    filter_gex<-tidyr::separate_rows(filter_gene, TPM)
    filter_gex$TPM <- as.numeric(filter_gex$TPM)
    filter_gex$Gene.ID_Tissue <- paste(filter_gex$Ensembl_ID, filter_gex$Anatomical_entity_name, sep = "_")
    CV_Tissue<-filter_gex %>% group_by(Gene.ID_Tissue) %>% summarise(CV = CV_function(TPM, na.rm=FALSE))
    dplyr::full_join(CV_Tissue,filter_gex)
    return(CV_Tissue)
    } 
  
  if (metric_type == "CV_Tissue"){
    CV_Tissue<-CV_Tissue(map_species,map_tissues)
    return(CV_Tissue)
  }
  else {
    stop("Error: Invalid metric_type in CoSIAn Object. Make sure the value given in the metric_type slot are avalible and in the correct format.")
  }
  
}) 
    
    

  
  