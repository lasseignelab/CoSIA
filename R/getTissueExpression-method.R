setMethod("getTissueExpression", signature(object = "CosiaExpressTissue"), function(object) {
    # user's input of the function
    if (object@tissues == "all tissues") {
        gene_species <- object@gene_species
        if (gene_species == "mus_musculus") {
          bgee_species <- filter_mouse
        }
        if (gene_species == "rattus_norvegicus") {
          bgee_species <- filter_rat
        }
        if (gene_species == "danio_rerio") {
          bgee_species <- filter_zebrafish
        }
        if (gene_species == "homo_sapiens") {
          bgee_species <- filter_human
        }
        species_specific <- data.frame(dplyr::select(bgee_species, Gene.ID, Experiment.ID, Anatomical.entity.ID, Anatomical.entity.name, Read.count,
            TPM, FPKM, Detection.flag))
        species_specific$Anatomical.entity.name <- gsub("\"", "", species_specific$Anatomical.entity.name)
        gene_specific_data <- dplyr::filter(species_specific, Gene.ID == object@single_gene)
        sample_size <- data.frame(table(gene_specific_data$Anatomical.entity.name))
        colnames(sample_size)[which(names(sample_size) == "Var1")] <- "Anatomical.entity.name"
        values <- aggregate(data = gene_specific_data, x = gene_specific_data$TPM, by = list(gene_specific_data$Anatomical.entity.name), FUN = median)
        colnames(values)[which(names(values) == "Group.1")] <- "Anatomical.entity.name"
        value <- merge(values, sample_size, by = "Anatomical.entity.name")
        value$AEN <- paste(value$Anatomical.entity.name, "(n=", value$Freq, ")", sep = "")
        gene_specific_data <- merge(value, gene_specific_data, by = "Anatomical.entity.name")

        brewer.pal.info <- RColorBrewer::brewer.pal.info
        palette5 <- RColorBrewer::brewer.pal.info[brewer.pal.info$category == "qual", ]
        color <- unlist(mapply(RColorBrewer::brewer.pal, palette5$maxcolors, rownames(palette5)))

        fig <- plotly::plot_ly(gene_specific_data, x = ~AEN, y = ~TPM, type = "scatter", mode = "markers", color = ~AEN, colors = color) %>%
            plotly::add_markers(x = ~AEN, y = ~x, marker = list(symbol = "line-ew", size = 10, line = list(color = "grey", width = 2)) %>%
                plotly::add_trace(marker = list(size = 8, line = list(color = "black", width = 0.5)), showlegend = F))
        fig <- fig %>%
            plotly::layout(title = stringr::str_wrap(paste("Gene Expression of the gene", object@single_gene, "in", object@gene_species, "amoung",
                object@tissues, sep = " ")), xaxis = list(title = "Anatomical Entity Name", size = 2), yaxis = list(title = "TPM (transcript per million)",
                zeroline = F), showlegend = FALSE)

        return(fig)
    } else {
        gene_species <- object@gene_species
        if (gene_species == "mus_musculus") {
          bgee_species <- mouse_specific
        }
        if (gene_species == "rattus_norvegicus") {
          bgee_species <- rat_specific
        }
        if (gene_species == "danio_rerio") {
          bgee_species <- zebrafish_specific
        }
        if (gene_species == "homo_sapiens") {
          bgee_species <- human_specific
        }
        gene_specific_data <- dplyr::filter(bgee_species, Gene.ID %in% object@single_gene)
        tissue_specific_data <- dplyr::filter(gene_specific_data, Anatomical.entity.name %in% object@tissues)
        sample_size <- data.frame(table(tissue_specific_data$Anatomical.entity.name))
        colnames(sample_size)[which(names(sample_size) == "Var1")] <- "Anatomical.entity.name"
        values <- aggregate(data = tissue_specific_data, x = tissue_specific_data$TPM, by = list(tissue_specific_data$Anatomical.entity.name),
            FUN = median)
        colnames(values)[which(names(values) == "Group.1")] <- "Anatomical.entity.name"
        value <- merge(values, sample_size, by = "Anatomical.entity.name")
        value$AEN <- paste(value$Anatomical.entity.name, "(n=", value$Freq, ")", sep = "")
        tissue_specific_data <- merge(value, tissue_specific_data, by = "Anatomical.entity.name")

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
    }
})

