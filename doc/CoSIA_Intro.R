## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval=FALSE--------------------------------------------------------------
#
#  #if (!requireNamespace("BiocManager", quietly=TRUE))
#      #install.packages("BiocManager")
#  #BiocManager::install("CoSIA")
#
#  #library(devtools)
#  #install_github("lasseignelab/CoSIA", ref= "main")
#

## ----setup, warning=FALSE, message=FALSE--------------------------------------

library(CoSIA)

load("~/Desktop/EH_Data.RData")
# load("../../Learning/useCoSIA/data/EH_Data.RData")
load("../inst/extdata/proccessed/monogenic_kidney_genes.rda")



## ----getTissues_1-------------------------------------------------------------

CoSIA::getTissues("d_rerio")


## ----getTissues_2-------------------------------------------------------------

CoSIA::getTissues(c("h_sapiens", "m_musculus", "r_norvegicus"))


## ----CoSIAnObj----------------------------------------------------------------

CoSIAn_Obj <- CoSIA::CoSIAn(
  gene_set = unique(monogenic_kidney_genes$Gene),
  i_species = "h_sapiens",
  o_species = c("h_sapiens", "m_musculus", "r_norvegicus"),
  input_id = "Symbol",
  output_ids = "Ensembl_id",
  map_species = c("h_sapiens", "m_musculus", "r_norvegicus"),
  map_tissues = c("adult mammalian kidney", "heart"),
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


## ----use2_2, message=FALSE, warning=FALSE-------------------------------------

CoSIAn_Obj_gexplot <- CoSIA::plotSpeciesGEx(CoSIAn_Obj_gex, "adult mammalian kidney", "ENSG00000171766")

CoSIAn_Obj_gexplot


## ----use3, message=FALSE, warning=FALSE---------------------------------------

CoSIAn_Obj_CV <- CoSIA::getGExMetrics(CoSIAn_Obj_gex)

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
