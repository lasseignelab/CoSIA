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
  filter_gex<-filter_gex[filter_gex$Anatomical_entity_name %in% single_tissue, ]
  filter_gex<-tidyr::separate_rows(filter_gex, TPM)
  filter_gex$TPM <- as.numeric(filter_gex$TPM)
  #add some validation methods here : check and make sure that the tissue is in gex the species in mapped species are in gex and that the single gene has been converted
  #color palette for plot (we can make this more modifiable)
  brewer.pal.info <- RColorBrewer::brewer.pal.info
  palette5 <- RColorBrewer::brewer.pal.info[brewer.pal.info$category == "qual", ]
  color <- unlist(mapply(RColorBrewer::brewer.pal, palette5$maxcolors, rownames(palette5)))
  fig <- filter_gex %>%
    plotly::plot_ly(
      x = ~Species,
      y = ~TPM,
      type = 'scatter',
      mode = "markers", 
      color = ~Species,
      colors = color
    ) %>%
    plotly::add_markers(x = ~Species, 
                        y = ~median(TPM), 
                        marker = list(symbol = "line-ew", 
                                      size = 20, 
                                      line = list(color = "grey"
                                                  , width = 2)))

  fig <- fig %>%
    plotly::add_trace(marker = list(size = 8, line = list(color = "black", width = 0.75)), 
                      showlegend = F)
  
  fig <- fig %>%
    plotly::add_trace(type = "violin", 
                      showlegend = F)
  fig <- fig %>%
    plotly::layout(xaxis = list(title = "Anatomical Entity Name", size = 2), 
                   yaxis = list(title = "TPM (transcript per million)",zeroline = F),
                   title = stringr::str_wrap(paste("Gene Expression of the gene", single_gene, "in", single_tissue , "across species" , sep = " ")),
                   showlegend = FALSE)

  return(fig)
})

#' plotTissueGEx Generic
#'
#' @param object
#'
#' @return
#' @export
#'
#' @examples

setGeneric("plotTissueGEx", function(object, single_species, single_gene) standardGeneric("plotTissueGEx"))

#' plotTissueGEx Method
#'
#' @param object CoSIAn. 
#'
#' @return 
#' @export
#'
#' @examples

setMethod("plotTissueGEx", signature(object = "CoSIAn"), function(object, single_species, single_gene) {
  gex_dataframe<-object@gex
  gex_dataframe<- as.data.frame(gex_dataframe)
  converted_id<-object@converted_id
  converted_id<- as.data.frame(converted_id)
  Species_SWITCH <- Vectorize(vectorize.args = "single_species", FUN = function(single_species) {
    switch(as.character(single_species), 
           h_sapiens ="Homo_sapiens",
           m_musculus = "Mus_musculus",
           r_norvegicus = "Rattus_norvegicus",
           d_rerio = "Danio_rerio",
           d_melanogaster = "Drosophila_melanogaster",
           c_elegans = "Caenorhabditis_elegans",
           stop("Error: Invalid map_species in CoSIAn Object. Make sure the species in the map_species slot are an avalible model 
                                    organism and are in the correct format.")
    )
  })
  
  single_species <- Species_SWITCH(single_species)
  filter_ids<-dplyr::select(converted_id, ends_with("ensembl_id"))
  ids <- filter_ids[Reduce(`|`, lapply(filter_ids, grepl, pattern = single_gene)),]
  ids<- as.character(ids)
  filter_gex<-gex_dataframe[gex_dataframe$Ensembl_ID %in% ids, ]
  filter_gex<-filter_gex[filter_gex$Species %in% single_species, ]
  filter_gex<-tidyr::separate_rows(filter_gex, TPM)
  filter_gex$TPM <- as.numeric(filter_gex$TPM)
  #add some validation methods here : check and make sure that the tissue is in gex the species in mapped species are in gex and that the single gene has been converted
  #color palette for plot (we can make this more modifiable)
  brewer.pal.info <- RColorBrewer::brewer.pal.info
  palette5 <- RColorBrewer::brewer.pal.info[brewer.pal.info$category == "qual", ]
  color <- unlist(mapply(RColorBrewer::brewer.pal, palette5$maxcolors, rownames(palette5)))
  fig <- filter_gex %>%
    plotly::plot_ly(
      x = ~Anatomical_entity_name,
      y = ~TPM,
      type = 'scatter',
      mode = "markers", 
      color = ~Anatomical_entity_name,
      colors = color
    ) %>%
    plotly::add_markers(x = ~Anatomical_entity_name, 
                        y = ~median(TPM), 
                        marker = list(symbol = "line-ew", 
                                      size = 20, 
                                      line = list(color = "grey"
                                                  , width = 2)))
  
  fig <- fig %>%
    plotly::add_trace(marker = list(size = 8, line = list(color = "black", width = 0.75)), 
                      showlegend = F)
  
  fig <- fig %>%
    plotly::add_trace(type = "violin", 
                      showlegend = F)
  fig <- fig %>%
    plotly::layout(xaxis = list(title = "Anatomical Entity Name", size = 2), 
                   yaxis = list(title = "TPM (transcript per million)",zeroline = F),
                   title = stringr::str_wrap(paste("Gene Expression of the gene", single_gene, "in", single_species , "across species" , sep = " ")),
                   showlegend = FALSE)
  
  return(fig)
})

  
    
