#' plotSpeciesGEx Generic
#'
#' @param object CoSIAn object with all user accessible slots filled in as well as the converted_id and gex slot filled
#' @param single_tissue one tissue that the user wants to investigate across the mapped species
#' @param single_gene one ensembl id that the user wants to investigate across the mapped species
#' @return initializes a generic function for plotSpeciesGEx as preparation for defining the plotSpeciesGEx Method
#'
#' @export
#' @examples
#' Kidney_Genes <- CoSIAn(
#'   gene_set = c("ENSG00000008710", "ENSG00000118762", "ENSG00000152217"),
#'   i_species = "h_sapiens", input_id = "Ensembl_id", o_species = c(
#'     "d_melanogaster", "m_musculus",
#'     "h_sapiens", "d_rerio", "c_elegans", "r_norvegicus"
#'   ), output_ids = c("Ensembl_id", "Symbol"),
#'   mapping_tool = "annotationDBI", ortholog_database = "HomoloGene", map_tissues = "heart",
#'   map_species = c("m_musculus"), metric_type = "DS_Gene"
#' )
#' Kidney_gene_conversion <- CoSIA::getConversion(Kidney_Genes)
#' Kidney_gene_gex <- getGEx(Kidney_gene_conversion)
#' plotSpeciesGEx(Kidney_gene_gex, "liver", "ENSG00000008710")
#' plotSpeciesGEx(Kidney_gene_gex, "brain", "ENSG00000118762")
setGeneric("plotSpeciesGEx", function(object, single_tissue, single_gene) standardGeneric("plotSpeciesGEx"))

#' plotSpeciesGEx Method
#'
#' @param object CoSIAn object with all user accessible slots filled in as well as the converted_id and gex slot filled
#' @param single_tissue one tissue that the user wants to investigate across the mapped species
#' @param single_gene one ensembl id that the user wants to investigate across the mapped species
#'
#' @return plot object
#' @export
#' @importFrom tidyselect ends_with
#'
#' @examples
#' Kidney_Genes <- CoSIAn(
#'   gene_set = c("ENSG00000008710", "ENSG00000118762", "ENSG00000152217"),
#'   i_species = "h_sapiens", input_id = "Ensembl_id", o_species = c(
#'     "d_melanogaster", "m_musculus",
#'     "h_sapiens", "d_rerio", "c_elegans", "r_norvegicus"
#'   ), output_ids = c("Ensembl_id", "Symbol"),
#'   mapping_tool = "annotationDBI", ortholog_database = "HomoloGene", map_tissues = "heart",
#'   map_species = c("m_musculus"), metric_type = "DS_Gene"
#' )
#' Kidney_gene_conversion <- CoSIA::getConversion(Kidney_Genes)
#' Kidney_gene_gex <- getGEx(Kidney_gene_conversion)
#' plotSpeciesGEx(Kidney_gene_gex, "liver", "ENSG00000008710")
#' plotSpeciesGEx(Kidney_gene_gex, "brain", "ENSG00000118762")
setMethod("plotSpeciesGEx", signature(object = "CoSIAn"), function(object, single_tissue, single_gene) {
  gex_dataframe <- object@gex
  gex_dataframe <- as.data.frame(gex_dataframe)
  converted_id <- object@converted_id
  converted_id <- as.data.frame(converted_id)
  filter_ids <- dplyr::select(converted_id, ends_with("ensembl_id"))
  ids <- filter_ids[Reduce(`|`, lapply(filter_ids, grepl, pattern = single_gene)), ]
  ids <- as.character(ids)
  filter_gex <- gex_dataframe[gex_dataframe$Ensembl_ID %in% ids, ]
  filter_gex <- filter_gex[filter_gex$Anatomical_entity_name %in% single_tissue, ]
  filter_gex <- tidyr::separate_rows(filter_gex, VST)
  filter_gex$VST <- as.numeric(filter_gex$VST)
  filter_gex <- as.data.frame(filter_gex)
  # add some validation methods here : check and make sure that the tissue is in gex the species in mapped species are
  # in gex and that the single gene has been converted color palette for plot (we can make this more modifiable)
  brewer.pal.info <- RColorBrewer::brewer.pal.info
  palette5 <- RColorBrewer::brewer.pal.info[brewer.pal.info$category == "qual", ]
  color <- unlist(mapply(RColorBrewer::brewer.pal, palette5$maxcolors, rownames(palette5)))
  fig <- filter_gex %>%
    plotly::plot_ly(x = ~Species, y = ~VST, type = "scatter", mode = "markers", color = ~Species, colors = color)
  fig <- fig %>%
    plotly::add_trace(marker = list(size = 8, line = list(color = "black", width = 0.75)), showlegend = FALSE)

  fig <- fig %>%
    plotly::add_trace(type = "violin", spanmode = "hard", showlegend = FALSE)
  fig <- fig %>%
    plotly::layout(xaxis = list(title = "Species", size = 2), yaxis = list(
      title = "VST (Variance Stabilized Transformation of Read Counts)",
      zeroline = FALSE
    ), title = stringr::str_wrap(paste("Gene Expression of the gene", single_gene, "in", single_tissue,
      "across species",
      sep = " "
    )), showlegend = FALSE)
  return(fig)
})

#' plotTissueGEx Generic
#'
#' @param object CoSIAn object with all user accessible slots filled in as well as the converted_id and gex slot filled
#' @param single_species one species that the user wants to investigate across the mapped tissues
#' @param single_gene one ensembl id that the user wants to investigate across the mapped tissues
#' @return initializes a generic function for plotTissueGEx as preparation for defining the plotTissueGEx Method
#'
#' @export
#' @examples
#' Kidney_Genes <- CoSIAn(
#'   gene_set = c("ENSG00000008710", "ENSG00000118762", "ENSG00000152217"),
#'   i_species = "h_sapiens", input_id = "Ensembl_id", o_species = c(
#'     "d_melanogaster", "m_musculus",
#'     "h_sapiens", "d_rerio", "c_elegans", "r_norvegicus"
#'   ), output_ids = c("Ensembl_id", "Symbol"),
#'   mapping_tool = "annotationDBI", ortholog_database = "HomoloGene", map_tissues = "heart",
#'   map_species = c("m_musculus"), metric_type = "DS_Gene"
#' )
#' Kidney_gene_conversion <- CoSIA::getConversion(Kidney_Genes)
#' Kidney_gene_gex <- getGEx(Kidney_gene_conversion)
#' plotTissueGEx(Kidney_gene_gex, "h_sapiens", "ENSG00000008710")
#' plotTissueGEx(Kidney_gene_gex, "m_musculus", "ENSG00000008710")
setGeneric("plotTissueGEx", function(object, single_species, single_gene) standardGeneric("plotTissueGEx"))

#' plotTissueGEx Method
#'
#' @param object CoSIAn object with all user accessible slots filled in as well as the converted_id and gex slot filled
#' @param single_species one species that the user wants to investigate across the mapped tissues
#' @param single_gene one ensembl id that the user wants to investigate across the mapped tissues
#'
#' @return plot object
#' @export
#' @importFrom tidyselect ends_with

#'
#' @examples
#' Kidney_Genes <- CoSIAn(
#'   gene_set = c("ENSG00000008710", "ENSG00000118762", "ENSG00000152217"),
#'   i_species = "h_sapiens", input_id = "Ensembl_id", o_species = c(
#'     "d_melanogaster", "m_musculus",
#'     "h_sapiens", "d_rerio", "c_elegans", "r_norvegicus"
#'   ), output_ids = c("Ensembl_id", "Symbol"),
#'   mapping_tool = "annotationDBI", ortholog_database = "HomoloGene", map_tissues = "heart",
#'   map_species = c("m_musculus"), metric_type = "DS_Gene"
#' )
#' Kidney_gene_conversion <- CoSIA::getConversion(Kidney_Genes)
#' Kidney_gene_gex <- getGEx(Kidney_gene_conversion)
#' plotTissueGEx(Kidney_gene_gex, "h_sapiens", "ENSG00000008710")
#' plotTissueGEx(Kidney_gene_gex, "m_musculus", "ENSG00000008710")
setMethod("plotTissueGEx", signature(object = "CoSIAn"), function(object, single_species, single_gene) {
  gex_dataframe <- object@gex
  gex_dataframe <- as.data.frame(gex_dataframe)
  converted_id <- object@converted_id
  converted_id <- as.data.frame(converted_id)
  Species_SWITCH <- Vectorize(vectorize.args = "single_species", FUN = function(single_species) {
    switch(as.character(single_species),
      h_sapiens = "Homo_sapiens",
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
  filter_ids <- dplyr::select(converted_id, ends_with("ensembl_id"))
  ids <- filter_ids[Reduce(`|`, lapply(filter_ids, grepl, pattern = single_gene)), ]
  ids <- as.character(ids)
  filter_gex <- gex_dataframe[gex_dataframe$Ensembl_ID %in% ids, ]
  filter_gex <- filter_gex[filter_gex$Species %in% single_species, ]
  filter_gex <- tidyr::separate_rows(filter_gex, VST)
  filter_gex$VST <- as.numeric(filter_gex$VST)

  # add some validation methods here : check and make sure that the tissue is in gex the species in mapped species are
  # in gex and that the single gene has been converted color palette for plot (we can make this more modifiable)
  brewer.pal.info <- RColorBrewer::brewer.pal.info
  palette5 <- RColorBrewer::brewer.pal.info[brewer.pal.info$category == "qual", ]
  color <- unlist(mapply(RColorBrewer::brewer.pal, palette5$maxcolors, rownames(palette5)))
  fig <- filter_gex %>%
    plotly::plot_ly(
      x = ~Anatomical_entity_name, y = ~VST, type = "scatter", mode = "markers", color = ~Anatomical_entity_name,
      colors = color
    )

  fig <- fig %>%
    plotly::add_trace(marker = list(size = 8, line = list(color = "black", width = 0.75)), showlegend = FALSE)

  fig <- fig %>%
    plotly::add_trace(type = "violin", spanmode = "hard", showlegend = FALSE)
  fig <- fig %>%
    plotly::layout(xaxis = list(title = "Anatomical Entity Name", size = 2), yaxis = list(
      title = "VST (Variance Stabilized Transformation of Read Counts)",
      zeroline = FALSE
    ), title = stringr::str_wrap(paste("Gene Expression of the gene", single_gene, "in", single_species,
      "across tissues",
      sep = " "
    )), showlegend = FALSE)

  return(fig)
})


#' plotDSGEx Generic
#'
#' @param object CoSIAn object with all user accessible slots filled in as well as the converted_id and metric slot filled
#' @return initializes a generic function for plotDSGEx as preparation for defining the plotDSGEx Method
#' @export
#' @examples
#' Kidney_Genes <- CoSIAn(
#'   gene_set = c("ENSG00000008710", "ENSG00000118762", "ENSG00000152217"),
#'   i_species = "h_sapiens", input_id = "Ensembl_id", o_species = c(
#'     "d_melanogaster", "m_musculus",
#'     "h_sapiens", "d_rerio", "c_elegans", "r_norvegicus"
#'   ), output_ids = c("Ensembl_id", "Symbol"),
#'   mapping_tool = "annotationDBI", ortholog_database = "HomoloGene", map_tissues = "heart",
#'   map_species = c("m_musculus"), metric_type = "DS_Gene"
#' )
#' Kidney_gene_conversion <- CoSIA::getConversion(Kidney_Genes)
#' Kidney_gene_metric <- getGExMetrics(Kidney_gene_conversion)
#' plotDSGEx(Kidney_gene_metric)
setGeneric("plotDSGEx", function(object) standardGeneric("plotDSGEx"))


#' plotDSGEx Method
#'
#' @param object CoSIAn object with all user accessible slots filled in as well as the converted_id and metric slot filled
#'
#' @return plot object
#' @export
#'
#' @examples
#' Kidney_Genes <- CoSIAn(
#'   gene_set = c("ENSG00000008710", "ENSG00000118762", "ENSG00000152217"),
#'   i_species = "h_sapiens", input_id = "Ensembl_id", o_species = c(
#'     "d_melanogaster", "m_musculus",
#'     "h_sapiens", "d_rerio", "c_elegans", "r_norvegicus"
#'   ), output_ids = c("Ensembl_id", "Symbol"),
#'   mapping_tool = "annotationDBI", ortholog_database = "HomoloGene", map_tissues = "heart",
#'   map_species = c("m_musculus"), metric_type = "DS_Gene"
#' )
#' Kidney_gene_conversion <- CoSIA::getConversion(Kidney_Genes)
#' Kidney_gene_metric <- getGExMetrics(Kidney_gene_conversion)
#' plotDSGEx(Kidney_gene_metric)
setMethod("plotDSGEx", signature(object = "CoSIAn"), function(object) {
  metric_type <- object@metric_type
  df_metric <- object@metric
  if (metric_type == "DS_Gene") {
    cols <- c(
      "#88CCEE", "#CC6677", "#DDCC77", "#332288", "#AA4499", "#44AA99", "#999933", "#882255", "#661100", "#117733",
      "#6699CC", "#888888"
    ) # make these colors color blind friendly
    DS_plot <- ggplot2::ggplot(df_metric, ggplot2::aes(x = Specificity, y = Diversity, color = Species))
    DS_plot <- DS_plot + ggplot2::geom_point(size = 3, ggplot2::aes()) + ggplot2::scale_color_manual(values = cols) + ggplot2::ggtitle("Diversity versus Specificity of Genes in Geneset \nAcross Mapped Tissues in a Species") +
      ggplot2::theme_classic()
  } else if (metric_type == "DS_Tissue") {
    cols <- c(
      "#88CCEE", "#CC6677", "#DDCC77", "#332288", "#AA4499", "#44AA99", "#999933", "#117733", "#882255", "#661100",
      "#6699CC", "#888888"
    ) # make these colors color blind friendly
    DS_plot <- ggplot2::ggplot(df_metric, ggplot2::aes(x = Specificity, y = Diversity, color = Anatomical_entity_name))
    DS_plot <- DS_plot + ggplot2::geom_point(size = 3, ggplot2::aes(shape = Species)) + ggplot2::scale_color_manual(values = cols) +
      ggplot2::ggtitle("Diversity versus Specificity of Genes in Geneset \nAcross Anatomical Entity Names") + ggplot2::theme_classic()
  } else if (metric_type == "DS_Tissue_all") {
    cols <- c(
      "#88CCEE", "#CC6677", "#DDCC77", "#332288", "#AA4499", "#44AA99", "#999933", "#117733", "#882255", "#661100",
      "#6699CC", "#888888"
    ) # make these colors color blind friendly
    DS_plot <- ggplot2::ggplot(df_metric, ggplot2::aes(x = Specificity, y = Diversity, color = Anatomical_entity_name))
    DS_plot <- DS_plot + ggplot2::geom_point(size = 3, ggplot2::aes(shape = Species)) + ggplot2::scale_color_manual(values = cols) +
      ggplot2::ggtitle("Diversity versus Specificity of All Genes \nAcross Anatomical Entity Names") + ggplot2::theme_classic()
  } else if (metric_type == "DS_Gene_all") {
    cols <- c(
      "#88CCEE", "#CC6677", "#DDCC77", "#332288", "#AA4499", "#44AA99", "#999933", "#882255", "#661100", "#117733",
      "#6699CC", "#888888"
    ) # make these colors color blind friendly
    DS_plot <- ggplot2::ggplot(df_metric, ggplot2::aes(x = Specificity, y = Diversity, color = Species))
    DS_plot <- DS_plot + ggplot2::geom_point(size = 3, ggplot2::aes()) + ggplot2::scale_color_manual(values = cols) + ggplot2::ggtitle("Diversity versus Specificity of Genes in Geneset \nAcross all Tissues in a Species") +
      ggplot2::theme_classic()
  } else {
    stop("Error: Invalid metric type for plotDS make sure you have a DS argument as the metric type and the values are saved in the metric slot before proceeding. ")
  }
  return(DS_plot)
})

#' plotCVGEx Generic
#'
#' @param object CoSIAn object with all user accessible slots filled in as well as the converted_id and metric slot filled
#'
#' @return initializes a generic function for plotCVGEx as preparation for defining the plotCVGEx Method
#' @export
#' @examples
#' Kidney_Genes <- CoSIAn(
#'   gene_set = c("ENSG00000008710", "ENSG00000118762", "ENSG00000152217"),
#'   i_species = "h_sapiens", input_id = "Ensembl_id", o_species = c(
#'     "d_melanogaster", "m_musculus",
#'     "h_sapiens", "d_rerio", "c_elegans", "r_norvegicus"
#'   ), output_ids = c("Ensembl_id", "Symbol"),
#'   mapping_tool = "annotationDBI", ortholog_database = "HomoloGene", map_tissues = "heart",
#'   map_species = c("h_sapiens", "m_musculus"), metric_type = "CV_Species"
#' )
#' Kidney_gene_conversion <- CoSIA::getConversion(Kidney_Genes)
#' Kidney_gene_metric <- getGExMetrics(Kidney_gene_conversion)
#' plotCVGEx(Kidney_gene_metric)
setGeneric("plotCVGEx", function(object) standardGeneric("plotCVGEx"))


#' plotCVGEx Method
#'
#' @param object CoSIAn object with all user accessible slots filled in as well as the converted_id and metric slot filled
#'
#' @return plot object
#' @export
#' @importFrom tibble remove_rownames
#'
#' @examples
#' Kidney_Genes <- CoSIAn(
#'   gene_set = c("ENSG00000008710", "ENSG00000118762", "ENSG00000152217"),
#'   i_species = "h_sapiens", input_id = "Ensembl_id", o_species = c(
#'     "d_melanogaster", "m_musculus",
#'     "h_sapiens", "d_rerio", "c_elegans", "r_norvegicus"
#'   ), output_ids = c("Ensembl_id", "Symbol"),
#'   mapping_tool = "annotationDBI", ortholog_database = "HomoloGene", map_tissues = "heart",
#'   map_species = c("h_sapiens", "m_musculus"), metric_type = "CV_Species"
#' )
#' Kidney_gene_conversion <- CoSIA::getConversion(Kidney_Genes)
#' Kidney_gene_metric <- getGExMetrics(Kidney_gene_conversion)
#' plotCVGEx(Kidney_gene_metric)
setMethod("plotCVGEx", signature(object = "CoSIAn"), function(object) {
  # set object slots into variables
  metric_type <- object@metric_type
  df_metric <- object@metric

  # Dumbbell Plot Function Plot df: CV_Species_Plot segments: whether to add segments (TRUE) or not (FALSE) text:
  # whether to add text (TRUE) or not (FALSE) pch: symbol pt.cex: size of the points segcol: color of the segment lwd:
  # width of the segment ... : additional arguments to be passed to dotchart function
  dumbbell <- function(df, segments = TRUE, text = FALSE, pch = 19, pt.cex = 1, segcol = 1, lwd = 1, main = "", metric = "",
                       ...) {
    n_col <- ncol(df) - 1
    n_row <- nrow(df)
    df_subset <- df[, -1]
    group <- rep(1, n_row)

    o <- sort.list(as.numeric(group), decreasing = TRUE)
    group <- group[o]
    offset <- cumsum(c(0, diff(as.numeric(group)) != 0))
    y <- 1L:n_row + 2 * offset

    graphics::dotchart(df[, 2],
      labels = df$ensembl_id, color = 1, xlim = range(0, 1) + c(0, 0.1), groups = group, pch = NA,
      pt.cex = pt.cex, main = main, xlab = "Coefficient of Variation (CV)", ylab = col1
    )

    if (segments == TRUE) {
      for (i in seq_len(n_row)) {
        segments(min(df_subset[i, ]), y[i], max(df_subset[i, ]), y[i], lwd = lwd, col = segcol)
      }
    }
    if (metric == "CV_Tissue") {
      col <- c("plum", "tomato", "paleturquoise", "lightcoral")
      pch <- c(19, 18, 17, 16, 15, 3, 4, 8, 7, 9, 10, 13)
      count <- 1
      for (i in seq_len(n_row)) {
        if (count / length(map_tissues) == 1) {
          count <- 1
        } else {
          count <- count + 1
        }
        for (j in seq_len(n_col)) {
          graphics::points(df_subset[i, j], y[i], pch = pch[count], cex = pt.cex, col = col[j])
        }
      }
    }
    if (metric == "CV_Species") {
      col <- c("plum", "tomato", "paleturquoise", "lightcoral")
      for (i in seq_len(n_row)) {
        for (j in seq_len(n_col)) {
          graphics::points(df_subset[i, j], y[i], pch = 19, cex = pt.cex, col = col[j])
        }
      }
    }

    # if(text == TRUE) { for(i in 1:length(v1)) { text(min(v2[i ], v1[i]) - 1.5, y[i], labels = min(v2[i], v1[i]))
    # text(max(v2[i], v1[i]) + 1.5, y[i], labels = max(v2[i], v1[i])) } }

    graphics::legend("topright",
      legend = colnames(df_subset), col = col, pch = 19, bty = "n", pt.cex = 1, cex = 1, text.col = "black",
      horiz = FALSE, inset = c(-0.1, -0.04)
    )
  }

  if (metric_type == "CV_Species") {
    # CV Species
    input_species <- object@i_species
    col1 <- paste(input_species, "ensembl_id", sep = "_")
    CV_Species_Plot <- df_metric %>%
      remove_rownames() %>%
      column_to_rownames(var = col1)
    CV_Species_Plot <- dplyr::select(CV_Species_Plot, contains("CV_species"))
    CV_Species_Plot <- tibble::rownames_to_column(CV_Species_Plot, "ensembl_id")
    colnames(CV_Species_Plot) <- gsub("_CV_species", "", colnames(CV_Species_Plot))
    CV_plot <- dumbbell(CV_Species_Plot,
      text = FALSE, segments = TRUE, pt.cex = 1, main = "CV across Species in Geneset",
      metric = "CV_Species"
    )
  } else if (metric_type == "CV_Tissue") {
    map_tissues <- object@map_tissues
    input_species <- object@i_species
    col1 <- paste(input_species, "ensembl_id", sep = "_")
    CV_Tissue_Plot <- df_metric %>%
      tidyr::unite("ensembl_id_AE_names", tidyselect::all_of(col1), Anatomical_entity_name, sep = "_", remove = FALSE)
    CV_Tissue_Plot <- CV_Tissue_Plot %>%
      remove_rownames() %>%
      column_to_rownames(var = "ensembl_id_AE_names")
    CV_Tissue_Plot <- dplyr::select(CV_Tissue_Plot, contains("CV_tissue"))
    CV_Tissue_Plot <- tibble::rownames_to_column(CV_Tissue_Plot, "ensembl_id_AE_names")
    colnames(CV_Tissue_Plot) <- gsub("_CV_tissue", "", colnames(CV_Tissue_Plot))
    CV_Tissue_Plot <- CV_Tissue_Plot[with(CV_Tissue_Plot, order(ensembl_id_AE_names)), ]
    CV_plot <- dumbbell(CV_Tissue_Plot,
      text = FALSE, segments = TRUE, pt.cex = 1, main = "CV across Tissues in Geneset",
      metric = "CV_Tissue", map_tissues = map_tissues
    )
  } else {
    stop("Error: Invalid metric type for plotCV make sure you have a CV argument as the metric type and the values are saved in the metric slot before proceeding. ")
  }
  return(CV_plot)
})
