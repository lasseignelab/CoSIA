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
  expression_dataframe<-object@gex
  expression_dataframe<- as.data.frame(expression_dataframe)
  metric_type<- object@metric_type
  return_caluclated_metric_data<- function(metric_type){
    GEx_metric_data<- data.frame(matrix(0))
    if (any(metric_type == "CV")) {
      bgee_species <- dplyr::filter(Experimental_Hub_File, Species == "Mus_musculus")
      gene_specific_data <- dplyr::filter(bgee_species, Ensembl_ID == id_dataframe$m_musculus)
      GEx_data <- rbind(GEx_data, gene_specific_data)
      GEx_data<-as.data.frame(GEx_data)
    }
    else if (any(map_species == "r_norvegicus_ensembl_id")) {
      bgee_species <- dplyr::filter(Experimental_Hub_File, Species == "Rattus_norvegicus")
      gene_specific_data <- dplyr::filter(bgee_species, Ensembl_ID == id_dataframe$r_norvegicus)
      GEx_data <- rbind(GEx_data, gene_specific_data)
      GEx_data<-as.data.frame(GEx_data)
    }
    else()
    return(GEx_metric_data)
  }
  GEx_data<-lapply(map_species,return_filtered_Gex_data)
  GEx_data<-as.data.frame(do.call(rbind, GEx_data))
  object@gex <- data.frame(GEx_data)
  return(object)
}
)