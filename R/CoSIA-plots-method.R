#' plotSpeciesGEx Generic
#'
#' @param object
#'
#' @export
#'

setGeneric("plotSpeciesGEx", function(object, single_tissue, single_gene) standardGeneric("plotSpeciesGEx"))

#' plotSpeciesGEx Method
#'
#' @param object CoSIAn. 
#'
#' @return plot object
#' @export
#'
#' @examples
#' plotSpeciesGEx(Kidney_gene_gex,"liver","ENSG00000008710")
#' plotSpeciesGEx(Kidney_gene_gex,"brain","ENSG00000118762")

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
  filter_gex<-tidyr::separate_rows(filter_gex, VST)
  filter_gex$VST <- as.numeric(filter_gex$VST)
  filter_gex<-as.data.frame(filter_gex)
  #add some validation methods here : check and make sure that the tissue is in gex the species in mapped species are in gex and that the single gene has been converted
  #color palette for plot (we can make this more modifiable)
  brewer.pal.info <- RColorBrewer::brewer.pal.info
  palette5 <- RColorBrewer::brewer.pal.info[brewer.pal.info$category == "qual", ]
  color <- unlist(mapply(RColorBrewer::brewer.pal, palette5$maxcolors, rownames(palette5)))
  fig <- filter_gex %>%
    plotly::plot_ly(
      x = ~Species,
      y = ~VST,
      type = 'scatter',
      mode = "markers", 
      color = ~Species,
      colors = color
    ) %>%
    plotly::add_markers(x = ~Species, 
                        y = ~Median_VST, 
                        marker = list(symbol = "line-ew", 
                                      size = 20, 
                                      line = list(color = "grey"
                                                  , width = 2)))
  fig <- fig %>%
    plotly::add_trace(marker = list(size = 8, line = list(color = "black", width = 0.75)), 
                      showlegend = F)
  
  fig <- fig %>%
    plotly::add_trace(type = "violin", spanmode = "hard",
                      showlegend = F)
  fig <- fig %>%
    plotly::layout(xaxis = list(title = "Species", size = 2), 
                   yaxis = list(title = "VST (Variance Stabilized Transformation of Read Counts)",zeroline = F),
                   title = stringr::str_wrap(paste("Gene Expression of the gene", single_gene, "in", single_tissue , "across species" , sep = " ")),
                   showlegend = FALSE)
  return(fig)
})

#' plotTissueGEx Generic
#'
#' @param object
#'
#' @export

setGeneric("plotTissueGEx", function(object, single_species, single_gene) standardGeneric("plotTissueGEx"))

#' plotTissueGEx Method
#'
#' @param object CoSIAn. 
#'
#' @return plot object
#' @export
#'
#' @examples
#' plotTissueGEx(Kidney_gene_gex,"h_sapiens","ENSG00000008710")
#' plotTissueGEx(Kidney_gene_gex,"m_musculus","ENSG00000008710")


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
  filter_gex<-tidyr::separate_rows(filter_gex, VST)
  filter_gex$VST <- as.numeric(filter_gex$VST)

  #add some validation methods here : check and make sure that the tissue is in gex the species in mapped species are in gex and that the single gene has been converted
  #color palette for plot (we can make this more modifiable)
  brewer.pal.info <- RColorBrewer::brewer.pal.info
  palette5 <- RColorBrewer::brewer.pal.info[brewer.pal.info$category == "qual", ]
  color <- unlist(mapply(RColorBrewer::brewer.pal, palette5$maxcolors, rownames(palette5)))
  fig <- filter_gex %>%
    plotly::plot_ly(
      x = ~Anatomical_entity_name,
      y = ~VST,
      type = 'scatter',
      mode = "markers", 
      color = ~Anatomical_entity_name,
      colors = color
    ) %>%
    plotly::add_markers(x = ~Anatomical_entity_name, 
                        y = ~Median_VST, 
                        marker = list(symbol = "line-ew", 
                                      size = 20, 
                                      line = list(color = "grey"
                                                  , width = 2)))
  
  fig <- fig %>%
    plotly::add_trace(marker = list(size = 8, line = list(color = "black", width = 0.75)), 
                      showlegend = F)
  
  fig <- fig %>%
    plotly::add_trace(type = "violin", spanmode = "hard",
                      showlegend = F)
  fig <- fig %>%
    plotly::layout(xaxis = list(title = "Anatomical Entity Name", size = 2), 
                   yaxis = list(title = "VST (Variance Stabilized Transformation of Read Counts)",zeroline = F),
                   title = stringr::str_wrap(paste("Gene Expression of the gene", single_gene, "in", single_species , "across tissues" , sep = " ")),
                   showlegend = FALSE)
  
  return(fig)
})


#' plotDSGEx Generic
#'
#' @param object
#'
#' @export
setGeneric("plotDSGEx", function(object) standardGeneric("plotDSGEx"))


#' plotDSGEx Method
#'
#' @param object CoSIAn. 
#'
#' @return plot object
#' @export
#'
#' @examples
#' plotDSGEx(Kidney_gene_metric)


setMethod("plotDSGEx", signature(object = "CoSIAn"), function(object) {
  metric_type<- object@metric_type
  df_metric<-object@metric
  if (metric_type == "DS_Gene"){
    cols <- c("#88CCEE", "#CC6677", "#DDCC77",  "#332288", "#AA4499","#44AA99", "#999933", "#882255", "#661100", "#117733", "#6699CC", "#888888") #make these colors color blind friendly
    DS_plot<-ggplot2::ggplot(df_metric, ggplot2::aes(x = Specificity, y = Diversity, color = Species))
    DS_plot<-DS_plot+
      ggplot2::geom_point(size =3, ggplot2::aes())+
      ggplot2::scale_color_manual(values = cols)+
      ggplot2::ggtitle("Diversity versus Specificity of Genes in Geneset \nAcross Mapped Tissues in a Species")+
      ggplot2::theme_classic()
  }
  else if (metric_type == "DS_Tissue"){
    cols <- c("#88CCEE", "#CC6677", "#DDCC77",  "#332288", "#AA4499","#44AA99", "#999933", "#117733","#882255", "#661100", "#6699CC", "#888888") #make these colors color blind friendly
    DS_plot<-ggplot2::ggplot(df_metric, ggplot2::aes(x = Specificity, y = Diversity, color = Anatomical_entity_name))
    DS_plot<-DS_plot+
      ggplot2::geom_point(size =3, ggplot2::aes(shape=Species))+
      ggplot2::scale_color_manual(values = cols)+
      ggplot2::ggtitle("Diversity versus Specificity of Genes in Geneset \nAcross Anatomical Entity Names")+
      ggplot2::theme_classic()
  }
  else if (metric_type == "DS_Tissue_all"){
    cols <- c("#88CCEE", "#CC6677", "#DDCC77",  "#332288", "#AA4499","#44AA99", "#999933", "#117733","#882255", "#661100", "#6699CC", "#888888") #make these colors color blind friendly
    DS_plot<-ggplot2::ggplot(df_metric, ggplot2::aes(x = Specificity, y = Diversity, color = Anatomical_entity_name))
    DS_plot<-DS_plot+
      ggplot2::geom_point(size =3, ggplot2::aes(shape=Species))+
      ggplot2::scale_color_manual(values = cols)+
      ggplot2::ggtitle("Diversity versus Specificity of All Genes \nAcross Anatomical Entity Names")+
      ggplot2::theme_classic()
  }
  else if (metric_type == "DS_Gene_all"){
    cols <- c("#88CCEE", "#CC6677", "#DDCC77",  "#332288", "#AA4499","#44AA99", "#999933", "#882255", "#661100", "#117733", "#6699CC", "#888888") #make these colors color blind friendly
    DS_plot<-ggplot2::ggplot(df_metric, ggplot2::aes(x = Specificity, y = Diversity, color = Species))
    DS_plot<-DS_plot+
      ggplot2::geom_point(size =3, ggplot2::aes())+
      ggplot2::scale_color_manual(values = cols)+
      ggplot2::ggtitle("Diversity versus Specificity of Genes in Geneset \nAcross all Tissues in a Species")+
      ggplot2::theme_classic()
  }
  else{
    stop("Error: Invalid metric type for plotDS make sure you have a DS argument as the metric type and the values are saved in the metric slot before proceeding. ")
  }
  return(DS_plot)
  
})

#' plotCVGEx Generic
#'
#' @param object
#'
#' @export
setGeneric("plotCVGEx", function(object) standardGeneric("plotCVGEx"))


#' plotCVGEx Method
#'
#' @param object CoSIAn. 
#'
#' @return plot object
#' @export
#'
#' @examples
#' plotCVGEx(Kidney_gene_metric)


setMethod("plotCVGEx", signature(object = "CoSIAn"), function(object) { #make multiple heat maps per species and then put them together
  metric_type<- object@metric_type
  df_metric<-object@metric
  if (metric_type == "CV_Species"){
    df_metric.wide <- tidyr::pivot_wider(df_metric, names_from = Ensembl_ID, values_from = 'CV_Species')
    df_metric.wide <- df_metric.wide %>% remove_rownames %>% tibble::column_to_rownames(var="Species")
    
    CV_plot<-heatmaply::heatmaply(df_metric.wide,
                         na.value = "white",
                         dendrogram = "none",
                         xlab = "Ensembl ID",
                         ylab = "Species",
                         main = "The Coeffecient of Variation of Gene Expression of a set of Genes across Species",
                         fontsize_row = 8,
                         fontsize_col = 7,
                        
        )
  }
  else if (metric_type == "CV_Tissue"){
    df_metric_remove_s<-subset(df_metric, select = -c(Species) )
    df_metric.wide <- tidyr::pivot_wider(df_metric_remove_s, names_from = Ensembl_ID, values_from = 'CV_Tissue')
    df_metric.wide <- df_metric.wide %>% remove_rownames %>% tibble::column_to_rownames(var="Anatomical_entity_name")
    CV_plot<-heatmaply::heatmaply(df_metric.wide,
                         na.value = "white",
                         dendrogram = "none",
                         xlab = "Ensembl ID",
                         ylab = "Anatomical_entity_name",
                         main = "The Coeffecient of Variation of Gene Expression of a set of Genes across Tissues",
                         fontsize_row = 8,
                         fontsize_col = 7,
    )
    
  }
  else{
    stop("Error: Invalid metric type for plotCV make sure you have a CV argument as the metric type and the values are saved in the metric slot before proceeding. ")
  }
  return(CV_plot)
  
})
