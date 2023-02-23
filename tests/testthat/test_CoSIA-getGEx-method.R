#library(testthat)
#library(covr)
#test_file("/Users/Nathan2/Desktop/lasseigneLab/Cosia_Testing/tests/testthat/test_CoSIA-getGEx-method.R")
#covr <- file_coverage("R/CoSIA-getGEx-method.R","tests/testthat/test_CoSIA-getGEx-method.R")
#covr
#report(covr)
source("/Users/Nathan2/Desktop/lasseigneLab/Cosia_Testing/R/CoSIA-class.R")
source("/Users/Nathan2/Desktop/lasseigneLab/Cosia_Testing/R/CoSIA-getConversions-method.R")
source("/Users/Nathan2/Desktop/lasseigneLab/Cosia_Testing/R/CoSIA-getGEx-method.R")
#load("/Users/Nathan2/Downloads/EH_Data.RData")
cosia_object <- CoSIAn(
  gene_set = "ENSG00000040531",
  i_species="h_sapiens",
  o_species="h_sapiens",
  input_id="Ensembl_id",
  o_ids="Ensembl_id",
  mapping_tool="annotationDBI",
  ortholog_database = "HomoloGene",
  map_tissues = "",
  map_species = "h_sapiens",
  metric_type=""
)
expect_true(all.equal(getGEx(getConversions(cosia_object))@gex$Median_VST, 9.243398,tolerance=.001))

cosia_object <- CoSIAn(
  gene_set = "ENSMUSG00000089948",
  i_species="m_musculus",
  o_species="m_musculus",
  input_id="Ensembl_id",
  o_ids="Ensembl_id",
  mapping_tool="annotationDBI",
  ortholog_database = "HomoloGene",
  map_tissues = "",
  map_species = "m_musculus",
  metric_type=""
)
getGEx(getConversions(cosia_object))@gex

cosia_object <- CoSIAn(
  gene_set = "ENSRNOG00000024972",
  i_species="r_norvegicus",
  o_species="r_norvegicus",
  input_id="Ensembl_id",
  o_ids="Ensembl_id",
  mapping_tool="annotationDBI",
  ortholog_database = "HomoloGene",
  map_tissues = "",
  map_species = "r_norvegicus",
  metric_type=""
)
getGEx(getConversions(cosia_object))@gex

