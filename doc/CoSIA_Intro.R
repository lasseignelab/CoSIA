## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>"
)

## ---- eval=FALSE--------------------------------------------------------------
#  library(devtools)
#  install_github("lasseignelab/CoSIA", ref = "main", auth_token = "<PAT>")

## ----setup, warning=FALSE, message=FALSE--------------------------------------
library(CoSIA)
load("../inst/extdata/proccessed/monogenic_kidney_genes.rda")

## ----getTissues_1, warning=FALSE, message=FALSE-------------------------------
CoSIA::getTissues("d_rerio")

## ----getTissues_2, warning=FALSE, message=FALSE-------------------------------
CoSIA::getTissues(c("h_sapiens", "m_musculus", "r_norvegicus"))

## ----CoSIAnObj----------------------------------------------------------------
CoSIAn_Obj <- CoSIA::CoSIAn(
    gene_set = unique(monogenic_kidney_genes$Gene),
    i_species = "h_sapiens",
    o_species = c(
        "h_sapiens",
        "m_musculus",
        "r_norvegicus"
    ),
    input_id = "Symbol",
    output_ids = "Ensembl_id",
    map_species = c(
        "h_sapiens",
        "m_musculus",
        "r_norvegicus"
    ),
    map_tissues = c(
        "adult mammalian kidney",
        "heart"
    ),
    mapping_tool = "annotationDBI",
    ortholog_database = "HomoloGene",
    metric_type = "CV_Species"
)

str(CoSIAn_Obj)

## ----use1,  message=FALSE, warning=FALSE--------------------------------------
CoSIAn_Obj_convert <- CoSIA::getConversion(CoSIAn_Obj)

str(CoSIAn_Obj_convert)

## ----use2_1,  message=FALSE, warning=FALSE------------------------------------
CoSIAn_Obj_gex <- CoSIA::getGEx(CoSIAn_Obj_convert)

str(CoSIAn_Obj_gex)

## ----plotSpeciesGEx, fig.small=TRUE, fig.cap = "Gene Expression of GATM in Kidney Across Species", message=FALSE, warning=FALSE----
#CoSIAn_Obj_gexplot <- CoSIA::plotSpeciesGEx(CoSIAn_Obj_gex, "adult mammalian kidney", "ENSG00000171766")

#CoSIAn_Obj_gexplot

## ----plotCVGEx, dpi=200, fig.height=15, fig.width=8, fig.cap = "Gene Expression Variability Across Species in Kidney Tissue", fig.wide=TRUE, message=FALSE, warning=FALSE----
# downsampling data for figure
monogenic_kid_sample <- dplyr::sample_n(monogenic_kidney_genes, 50)

CoSIAn_Obj_CV <- CoSIA::CoSIAn(
    gene_set = unique(monogenic_kid_sample$Gene),
    i_species = "h_sapiens",
    o_species = c(
        "h_sapiens",
        "m_musculus",
        "r_norvegicus"
    ),
    input_id = "Symbol",
    output_ids = "Ensembl_id",
    map_species = c(
        "h_sapiens",
        "m_musculus",
        "r_norvegicus"
    ),
    map_tissues = c(
        "adult mammalian kidney",
        "heart"
    ),
    mapping_tool = "annotationDBI",
    ortholog_database = "HomoloGene",
    metric_type = "CV_Species"
)

CoSIAn_Obj_CV <- CoSIA::getConversion(CoSIAn_Obj_CV)

CoSIAn_Obj_CV <- CoSIA::getGExMetrics(CoSIAn_Obj_CV)

CoSIAn_Obj_CVplot <- CoSIA::plotCVGEx(CoSIAn_Obj_CV)

CoSIAn_Obj_CVplot

## ----use4, message=FALSE, warning=FALSE---------------------------------------
CoSIAn_Obj_DS <- CoSIA::CoSIAn(
    gene_set = unique(monogenic_kidney_genes$Gene),
    i_species = "h_sapiens",
    o_species = c("h_sapiens", "m_musculus", "r_norvegicus"),
    input_id = "Symbol",
    output_ids = "Ensembl_id",
    map_species = c("h_sapiens", "m_musculus", "r_norvegicus"),
    map_tissues = c("adult mammalian kidney", "heart"),
    mapping_tool = "annotationDBI",
    ortholog_database = "HomoloGene",
    metric_type = "DS_Tissue"
)

CoSIAn_Obj_DS <- CoSIA::getConversion(CoSIAn_Obj_DS)

CoSIAn_Obj_DS <- CoSIA::getGExMetrics(CoSIAn_Obj_DS)

CoSIAn_Obj_DSplot <- CoSIA::plotDSGEx(CoSIAn_Obj_DS)

CoSIAn_Obj_DSplot

## -----------------------------------------------------------------------------
sessionInfo()

