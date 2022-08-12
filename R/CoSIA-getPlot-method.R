# brewer.pal.info <- RColorBrewer::brewer.pal.info
# palette5 <- RColorBrewer::brewer.pal.info[brewer.pal.info$category == "qual", ]
# color <- unlist(mapply(RColorBrewer::brewer.pal, palette5$maxcolors, rownames(palette5)))
# 
# fig <- plotly::plot_ly(gene_specific_data, x = ~AEN, y = ~TPM, type = "scatter", mode = "markers", color = ~AEN, colors = color) %>%
#   plotly::add_markers(x = ~AEN, y = ~x, marker = list(symbol = "line-ew", size = 10, line = list(color = "grey", width = 2)) %>%
#                         plotly::add_trace(marker = list(size = 8, line = list(color = "black", width = 0.5)), showlegend = F))
# fig <- fig %>%
#   plotly::layout(title = stringr::str_wrap(paste("Gene Expression of the gene", object@single_gene, "in", object@gene_species, "amoung",
#                                                  object@tissues, sep = " ")), xaxis = list(title = "Anatomical Entity Name", size = 2), yaxis = list(title = "TPM (transcript per million)",
#                                                                                                                                                      zeroline = F), showlegend = FALSE)
# 
# return(fig)
# 
# 
# brewer.pal.info <- RColorBrewer::brewer.pal.info
# palette5 <- RColorBrewer::brewer.pal.info[brewer.pal.info$category == "qual", ]
# color <- unlist(mapply(RColorBrewer::brewer.pal, palette5$maxcolors, rownames(palette5)))
# 
# fig <- plotly::plot_ly(tissue_specific_data, x = ~AEN, y = ~TPM, type = "scatter", mode = "markers", color = ~AEN, colors = color) %>%
#   plotly::add_markers(x = ~AEN, y = ~x, marker = list(symbol = "line-ew", size = 10, line = list(color = "grey", width = 2)) %>%
#                         plotly::add_trace(marker = list(size = 8, line = list(color = "black", width = 0.5)), showlegend = F))
# fig <- fig %>%
#   plotly::layout(title = stringr::str_wrap(paste("Gene Expression of the gene", object@single_gene, "in", object@gene_species, "amoung",
#                                                  object@tissues, sep = " ")), xaxis = list(title = "Anatomical Entity Name", size = 2), yaxis = list(title = "TPM (transcript per million)",
#                                                                                                                                                      zeroline = F), showlegend = FALSE)
# 
# return(fig)