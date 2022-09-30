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
  filter_ids %>% tibble::as_tibble() %>% dplyr::filter_all(any_vars(. %in% single_gene))
  
  return(filter_ids)
  
  
  map_species<- object@map_species
  single_tissue<- single_tissue
  single_gene<- single_gene
  if (converted_id %in% single_gene ){
    stop("Single gene is not in converted_id. Use getConversions to get the ensembl_id and ortholog_ensembl ids for species in the mapped_species slot")
  }
  converted_id<-converted_id %>% dplyr::filter(stringr::str_detect(single_gene))
  return(converted_id)
  
  #filter gex_dataframe for single gene and single tissue
  filtered_gex_dataframe<- dplyr::filter(gex_dataframe, Anatomical_entity_name == single_tissue & Ensembl_ID == list_of_ensembl_ids)
  #color palette for plot (we can make this more modifiable)
  brewer.pal.info <- RColorBrewer::brewer.pal.info
  palette5 <- RColorBrewer::brewer.pal.info[brewer.pal.info$category == "qual", ]
  color <- unlist(mapply(RColorBrewer::brewer.pal, palette5$maxcolors, rownames(palette5)))
  #build plot
  
  
  
  fig <- plotly::plot_ly(filtered_gex_dataframe, x = ~Species, y = ~Median_TPM, type = "scatter", mode = "markers", color = ~Anatomical_entity_name, colors = color) %>%
  
  plotly::add_markers(x = ~Anatomical_entity_name, y = ~Species, marker = list(symbol = "line-ew", size = 10, line = list(color = "grey", width = 2)) %>%
                        plotly::add_trace(marker = list(size = 8, line = list(color = "black", width = 0.5)), showlegend = F))
fig <- fig %>%
  plotly::layout(title = stringr::str_wrap(paste("Gene Expression of the gene amoung",
                                                 single_tissue, sep = " ")), xaxis = list(title = "Anatomical Entity Name", size = 2), yaxis = list(title = "TPM (transcript per million)",
                                                                                                                                                     zeroline = F), showlegend = FALSE)

  return(fig)


brewer.pal.info <- RColorBrewer::brewer.pal.info
palette5 <- RColorBrewer::brewer.pal.info[brewer.pal.info$category == "qual", ]
color <- unlist(mapply(RColorBrewer::brewer.pal, palette5$maxcolors, rownames(palette5)))

fig <- plotly::plot_ly(tissue_specific_data, x = ~AEN, y = ~TPM, type = "scatter", mode = "markers", color = ~AEN, colors = color) %>%
  plotly::add_markers(x = ~AEN, y = ~x, marker = list(symbol = "line-ew", size = 10, line = list(color = "grey", width = 2)) %>%
                        plotly::add_trace(marker = list(size = 8, line = list(color = "black", width = 0.5)), showlegend = F))
fig <- fig %>%
  plotly::layout(title = stringr::str_wrap(paste("Gene Expression of the gene", object@single_gene, "in", object@gene_species, "amoung",
                                                 object@tissues, sep = " ")), xaxis = list(title = "Anatomical Entity Name", size = 2), yaxis = list(title = "TPM (transcript per million)",
                                                                                                                                                     zeroline = F), showlegend = FALSE)

return(fig)

})
