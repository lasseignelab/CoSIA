#' plotSpeciesGEx Generic
#'
#' @param object
#'
#' @return
#' @export
#'
#' @examples

setGeneric("plotSpeciesGEx", function(object, single_tissue, single_gene) standardGeneric("plotSpeciesGEx"))

#' plotSpeciesGEx Method
#'
#' @param object CoSIAn. 
#'
#' @return 
#' @export
#'
#' @examples

setMethod("plotSpeciesGEx", signature(object = "CoSIAn"), function(object, single_tissue, single_gene) {
  gex_dataframe<-object@gex
  gex_dataframe<- as.data.frame(gex_dataframe)
  converted_id<-object@converted_id
  converted_id<- as.data.frame(converted_id)
  filter_ids<-dplyr::select(converted_id, ends_with("ensembl_id"))
  ids <- filter_ids[Reduce(`|`, lapply(filter_ids, grepl, pattern = single_gene)),]
  ids<- as.character(ids)
  filter_gex<-gex_dataframe[gex_dataframe$Ensembl_ID %in% ids, ]
  filter_gex<-gex_dataframe[gex_dataframe$Anatomical_entity_name %in% single_tissue, ]
  filter_gex<-separate_rows(filter_gex, TPM)
  #add some validation methods here : check and make sure that the tissue is in gex the species in mapped species are in gex and that the single gene has been converted
  #color palette for plot (we can make this more modifiable)
  brewer.pal.info <- RColorBrewer::brewer.pal.info
  palette5 <- RColorBrewer::brewer.pal.info[brewer.pal.info$category == "qual", ]
  color <- unlist(mapply(RColorBrewer::brewer.pal, palette5$maxcolors, rownames(palette5)))
  
  fig <- df %>%
    plot_ly(
      y = ~total_bill,
      type = 'violin',
      box = list(
        visible = T
      ),
      meanline = list(
        visible = T
      ),
      x0 = 'Total Bill'
    ) 
  
  fig <- fig %>%
    layout(
      yaxis = list(
        title = "",
        zeroline = F
      )
    )
  
  fig
  
  
  fig <- plotly::plot_ly(filter_gex, x = ~Species, y = ~TPM, type = "scatter", mode = "markers", color = ~Species, colors = color)
  fig <- fig %>%
    plotly::add_trace(marker = list(size = 8, line = list(color = "black", width = 0.75)), showlegend = F)
  fig <- fig %>%
    plotly::layout(title = stringr::str_wrap(paste("Gene Expression of the gene", single_gene, "in", single_tissue, "across Species",
                                                   sep = " ")), xaxis = list(title = "Anatomical Entity Name", size = 2), yaxis = list(title = "TPM (transcript per million)", zeroline = F),
                   showlegend = FALSE)
  return(fig)
})
  
  
    
