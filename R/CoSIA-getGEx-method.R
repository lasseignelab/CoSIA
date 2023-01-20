#' getGEx Generic
#'
#' @param object CoSIAn object with all user accessible slots filled with converted_id slot filled
#'
#' @export
#' @examples
#' Kidney_Genes<-CoSIAn(gene_set = c("ENSG00000008710","ENSG00000118762","ENSG00000152217"),
#' i_species = "h_sapiens",input_id = "Ensembl_id",o_species = c("d_melanogaster","m_musculus",
#' "h_sapiens", "d_rerio","c_elegans","r_norvegicus"),output_ids = c("Ensembl_id","Symbol"), 
#' mapping_tool = "annotationDBI",ortholog_database = "HomoloGene",map_tissues = "heart", 
#' map_species = c("m_musculus"),metric_type = "DS_Gene")
#' Kidney_gene_conversion<-CoSIA::getConversion(Kidney_Genes)
#' load("~/Desktop/EH_Data.RData")
#' Kidney_gene_gex<- getGEx(Kidney_gene_conversion)

setGeneric("getGEx", function(object) standardGeneric("getGEx"))
#CoSIAn getGEx


#' getGEx Method
#'
#' @param object CoSIAn object with all user accessible slots filled with converted_id slot filled 
#'
#' @return CoSIAn Object with gex slot filled
#' @export
#'
#' @examples
#' Kidney_Genes<-CoSIAn(gene_set = c("ENSG00000008710","ENSG00000118762","ENSG00000152217"),
#' i_species = "h_sapiens",input_id = "Ensembl_id",o_species = c("d_melanogaster","m_musculus",
#' "h_sapiens", "d_rerio","c_elegans","r_norvegicus"),output_ids = c("Ensembl_id","Symbol"), 
#' mapping_tool = "annotationDBI",ortholog_database = "HomoloGene",map_tissues = "heart", 
#' map_species = c("m_musculus"),metric_type = "DS_Gene")
#' Kidney_gene_conversion<-CoSIA::getConversion(Kidney_Genes)
#' load("~/Desktop/EH_Data.RData")
#' Kidney_gene_gex<- getGEx(Kidney_gene_conversion)


setMethod("getGEx", signature(object = "CoSIAn"), function(object) {
  #add validity that converted_id is filled in
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
  id_dataframe<- dplyr::select(id_dataframe,tidyselect::matches(map_species))
  #load EH_Data here (get it off of cheaha soon the EH method will be used to pull the data here)
  return_filtered_Gex_data<- function(map_species){
    GEx_data<- data.frame(matrix(ncol = 7, nrow = 0))
    colnames(GEx_data) <- c('Anatomical_entity_name', 'Ensembl_ID', 'Sample_size','VST','Experiment_ID','Anatomical_entity_ID','Species')
    
    if (any(map_species == "m_musculus_ensembl_id")) {
      bgee_species <- dplyr::filter(Experimental_Hub_File, Species == "Mus_musculus")
      m_ensembl_id<- id_dataframe$m_musculus
      gene_specific_data <- dplyr::filter(bgee_species, Ensembl_ID %in% m_ensembl_id)
      GEx_data <- rbind(GEx_data, gene_specific_data)
      GEx_data<-as.data.frame(GEx_data)
    }
    else if (any(map_species == "r_norvegicus_ensembl_id")) {
      bgee_species <- dplyr::filter(Experimental_Hub_File, Species == "Rattus_norvegicus")
      r_ensembl_id<- id_dataframe$r_norvegicus
      gene_specific_data <- dplyr::filter(bgee_species, Ensembl_ID %in% r_ensembl_id)
      GEx_data <- rbind(GEx_data, gene_specific_data)
      GEx_data<-as.data.frame(GEx_data)
    }
    else if (any(map_species == "d_rerio_ensembl_id")) {
      bgee_species <- dplyr::filter(Experimental_Hub_File, Species == "Danio_rerio")
      dr_ensembl_id<- id_dataframe$d_rerio
      gene_specific_data <- dplyr::filter(bgee_species, Ensembl_ID %in% dr_ensembl_id)
      GEx_data <- rbind(GEx_data, gene_specific_data)
      GEx_data<-as.data.frame(GEx_data)
    }
    else if (any(map_species == "h_sapiens_ensembl_id")) {
      bgee_species <- dplyr::filter(Experimental_Hub_File, Species == "Homo_sapiens")
      hs_ensembl_id<- id_dataframe$h_sapiens
      gene_specific_data <- dplyr::filter(bgee_species, Ensembl_ID %in% hs_ensembl_id)
      GEx_data <- rbind(GEx_data, gene_specific_data)
      GEx_data<-as.data.frame(GEx_data)
    }
    else if (any(map_species == "c_elegans_ensembl_id")) {
      bgee_species <- dplyr::filter(Experimental_Hub_File, Species == "Caenorhabditis_elegans")
      c_ensembl_id<- id_dataframe$c_elegans
      gene_specific_data <- dplyr::filter(bgee_species, Ensembl_ID %in% c_ensembl_id)
      GEx_data <- rbind(GEx_data, gene_specific_data)
      GEx_data<-as.data.frame(GEx_data)
    }
    else if (any(map_species == "d_melanogaster_ensembl_id")) {
      bgee_species <- dplyr::filter(Experimental_Hub_File, Species == "Drosophila_melanogaster")
      dm_ensembl_id<- id_dataframe$d_melanogaster
      gene_specific_data <- dplyr::filter(bgee_species, Ensembl_ID %in% dm_ensembl_id)
      GEx_data <- rbind(GEx_data, gene_specific_data)
      GEx_data<-as.data.frame(GEx_data)
    }
    else{
      stop("Error: map_species in CoSIAn Object. Make sure the species in the map_species slot are avalible model 
                                    organisms and are in the correct format.")
    }
    return(GEx_data)
  }
  GEx_data<-lapply(map_species,return_filtered_Gex_data)
  GEx_data<-as.data.frame(do.call(rbind, GEx_data))
  GEx_data <- dplyr::filter(GEx_data, Anatomical_entity_name %in% map_tissues)
  rownames(GEx_data) <- NULL
  object@gex <- data.frame(GEx_data)
  return(object)
})
