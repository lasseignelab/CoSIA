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
  return_calculated_metric_data<- function(metric_type){
    GEx_metric_data<- data.frame(matrix(0))
    if (any(metric_type == "CV")) {
    #   CV_function <- function(x, na.rm = FALSE){
    #     stopifnot(is.numeric(x))
    #     if(na.rm) x <- x[!is.na(x)]
    #     if(any(x < 0)) #changed this from <= to < 
    #       stop("Your data must be greater than zero!")
    #     sd(x, na.rm = FALSE)/mean(x)
    #   }
    #   CV<-expression_dataframe %>% group_by(Gene.ID_Tissue) %>% summarise(SUM = sum(TPM),Sample_Size = n(),Median = median(TPM), CV = CV_function(TPM, na.rm=FALSE))
     }
    else if (any(metric_type == "Entropy_Specificity")) {
      bgee_species <- dplyr::filter(Experimental_Hub_File, Species == "Rattus_norvegicus")
      gene_specific_data <- dplyr::filter(bgee_species, Ensembl_ID == id_dataframe$r_norvegicus)
      GEx_data <- rbind(GEx_data, gene_specific_data)
      GEx_data<-as.data.frame(GEx_data)
    }
    else{
    stop("Error: metric_type in CoSIAn Object. Make sure the metric type in the metric_type slot are avalible metrics and are written in the correct format.")
    }
    return(GEx_metric_data)
  }
  GEx_data<-lapply(map_species,return_filtered_Gex_data)
  GEx_data<-as.data.frame(do.call(rbind, GEx_data))
  object@gex <- data.frame(GEx_data)
  return(object)
}
)