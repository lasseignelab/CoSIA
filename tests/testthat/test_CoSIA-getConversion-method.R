#library(testthat)
#library(covr)
#test_file("/Users/Nathan2/Desktop/lasseigneLab/Cosia_Testing/tests/testthat/test_CoSIA-getConversions-method.R")
#covr <- file_coverage("R/CoSIA-getConversions-method.R","tests/testthat/test_CoSIA-getConversions-method.R")
#covr
#report(covr)
source("/Users/Nathan2/Desktop/lasseigneLab/Cosia_Testing/R/CoSIA-class.R")
source("/Users/Nathan2/Desktop/lasseigneLab/Cosia_Testing/R/CoSIA-getConversions-method.R")
input_id <- "Ensembl_id"#CHECK THIS LATER for some reason I can't instantiate the CoSIAn without it



test_that("no orthologs",{
  cosia_object <- CoSIAn(
    gene_set = "ENSG00000265301",
    i_species="h_sapiens",
    o_species="m_musculus",
    input_id="Ensembl_id",
    o_ids= "Ensembl_id",
    mapping_tool="annotationDBI",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  expect_warning(getConversions(cosia_object)@converted_id,"No orthologs were found for all genes in one of the species.")
  
  cosia_object <- CoSIAn(
    gene_set = "ENSMUSG00000075968",
    i_species="m_musculus",
    o_species="h_sapiens",
    input_id="Ensembl_id",
    o_ids= "Ensembl_id",
    mapping_tool="annotationDBI",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  expect_warning(getConversions(cosia_object)@converted_id,"No orthologs were found for all genes in one of the species.")
  
  cosia_object <- CoSIAn(
    gene_set = "ENSRNOG00000014951",
    i_species="r_norvegicus",
    o_species="d_rerio",
    input_id="Ensembl_id",
    o_ids= "Ensembl_id",
    mapping_tool="annotationDBI",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )

  expect_warning(getConversions(cosia_object)@converted_id,"No orthologs were found for all genes in one of the species.")
 
  
  cosia_object <- CoSIAn(
    gene_set = "ENSDARG00000023290",
    i_species="d_rerio",
    o_species="r_norvegicus",
    input_id="Ensembl_id",
    o_ids= "Ensembl_id",
    mapping_tool="annotationDBI",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )

    expect_warning(getConversions(cosia_object)@converted_id,"No orthologs were found for all genes in one of the species.")
 
  
  cosia_object <- CoSIAn(
    gene_set = "WBGene00017651",
    i_species="c_elegans",
    o_species="r_norvegicus",
    input_id="Ensembl_id",
    o_ids= "Ensembl_id",
    mapping_tool="annotationDBI",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
 
    expect_warning(getConversions(cosia_object)@converted_id,"No orthologs were found for all genes in one of the species.")
  
  
  cosia_object <- CoSIAn(
    gene_set = "FBgn0003975",
    i_species="d_melanogaster",
    o_species="r_norvegicus",
    input_id="Ensembl_id",
    o_ids= "Ensembl_id",
    mapping_tool="annotationDBI",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
 
    expect_warning(getConversions(cosia_object)@converted_id,"No orthologs were found for all genes in one of the species.")
 
})

test_that("clutter cleanup",{
  cosia_object <- CoSIAn(
    gene_set = "ENSG00000040531",
    i_species="h_sapiens",
    o_species="m_musculus",
    input_id="Ensembl_id",
    o_ids="MGI_Symbol",
    mapping_tool="biomaRt",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  h_sapiens_ensembl_id <- "ENSG00000040531"
  m_musculus_mgi_symbol <- "Ctns"
  expect_identical(getConversions(cosia_object)@converted_id,data.frame(h_sapiens_ensembl_id,m_musculus_mgi_symbol))
  
  cosia_object <- CoSIAn(
    gene_set = "ENSRNOG00000060666",
    i_species="r_norvegicus",
    o_species="d_melanogaster",
    input_id="Ensembl_id",
    o_ids= "Ensembl_id",
    mapping_tool="biomaRt",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  
  expect_warning(getConversions(cosia_object)@converted_id,"No orthologs were found for all genes in one of the species.")

  
  cosia_object <- CoSIAn(
    gene_set = c("ENSG00000265301","ENSG00000040531"),
    i_species="h_sapiens",
    o_species="m_musculus",
    input_id="Ensembl_id",
    o_ids= "Ensembl_id",
    mapping_tool="annotationDBI",
    ortholog_database = "NCBIOrtho",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  h_sapiens_ensembl_id <- c("ENSG00000265301","ENSG00000040531")
  m_musculus_ensembl_id <- c(NA,"ENSMUSG00000005949")
  expect_identical( getConversions(cosia_object)@converted_id,data.frame(h_sapiens_ensembl_id,m_musculus_ensembl_id))
  
  expect_error(remove_version_numbers("",""),"Error. Invalid species. Make sure it matches the proper format.")
})

##Annotation DBI
test_that("Check annotationDBI",{
  #humans -> other species
  cosia_object <- CoSIAn(
    gene_set = "ENSG00000040531",
    i_species="h_sapiens",
    o_species="h_sapiens",
    input_id="Ensembl_id",
    o_ids="Ensembl_id",
    mapping_tool="annotationDBI",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  h_sapiens_ensembl_id <- "ENSG00000040531"
  expect_identical(getConversions(cosia_object)@converted_id,data.frame(h_sapiens_ensembl_id))

  cosia_object <- CoSIAn(
    gene_set = "ENSG00000040531",
    i_species="h_sapiens",
    o_species="m_musculus",
    input_id="Ensembl_id",
    o_ids="Ensembl_id",
    mapping_tool="annotationDBI",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  h_sapiens_ensembl_id <- "ENSG00000040531"
  m_musculus_ensembl_id <- "ENSMUSG00000005949"
  expect_identical(getConversions(cosia_object)@converted_id,data.frame(h_sapiens_ensembl_id, m_musculus_ensembl_id))

  cosia_object <- CoSIAn(
    gene_set = "ENSG00000040531",
    i_species="h_sapiens",
    o_species="r_norvegicus",
    input_id="Ensembl_id",
    o_ids="Ensembl_id",
    mapping_tool="annotationDBI",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  h_sapiens_ensembl_id <- "ENSG00000040531"
  r_norvegicus_ensembl_id <- "ENSRNOG00000028688"
  expect_identical(getConversions(cosia_object)@converted_id,data.frame(h_sapiens_ensembl_id,r_norvegicus_ensembl_id))

  cosia_object <- CoSIAn(
    gene_set = "ENSG00000040531",
    i_species="h_sapiens",
    o_species="d_rerio",
    input_id="Ensembl_id",
    o_ids="Ensembl_id",
    mapping_tool="annotationDBI",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  h_sapiens_ensembl_id <- "ENSG00000040531"
  d_rerio_ensembl_id <- "ENSDARG00000008890"
  expect_identical(getConversions(cosia_object)@converted_id,data.frame(h_sapiens_ensembl_id,d_rerio_ensembl_id))

  cosia_object <- CoSIAn(
    gene_set = "ENSG00000040531",
    i_species="h_sapiens",
    o_species="c_elegans",
    input_id="Ensembl_id",
    o_ids="Ensembl_id",
    mapping_tool="annotationDBI",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  h_sapiens_ensembl_id <- "ENSG00000040531"
  c_elegans_ensembl_id <- "WBGene00008052"
  expect_identical(getConversions(cosia_object)@converted_id,data.frame(h_sapiens_ensembl_id,c_elegans_ensembl_id))

  cosia_object <- CoSIAn(
    gene_set = "ENSG00000040531",
    i_species="h_sapiens",
    o_species="d_melanogaster",
    input_id="Ensembl_id",
    o_ids="Ensembl_id",
    mapping_tool="annotationDBI",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  h_sapiens_ensembl_id <- "ENSG00000040531"
  d_melanogaster_ensembl_id <- "FBgn0039045"
  expect_identical(getConversions(cosia_object)@converted_id,data.frame(h_sapiens_ensembl_id, d_melanogaster_ensembl_id ))


  #mouse -> other species
  cosia_object <- CoSIAn(
    gene_set = "ENSMUSG00000002413",
    i_species="m_musculus",
    o_species="h_sapiens",
    input_id="Ensembl_id",
    o_ids="Ensembl_id",
    mapping_tool="annotationDBI",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  m_musculus_ensembl_id <- "ENSMUSG00000002413"
  h_sapiens_ensembl_id <- "ENSG00000157764"
  expect_identical(getConversions(cosia_object)@converted_id,data.frame(m_musculus_ensembl_id,h_sapiens_ensembl_id))

  cosia_object <- CoSIAn(
    gene_set = "ENSMUSG00000002413",
    i_species="m_musculus",
    o_species="m_musculus",
    input_id="Ensembl_id",
    o_ids="Ensembl_id",
    mapping_tool="annotationDBI",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  m_musculus_ensembl_id <- "ENSMUSG00000002413"
  expect_identical(getConversions(cosia_object)@converted_id,data.frame(m_musculus_ensembl_id))

  cosia_object <- CoSIAn(
    gene_set = "ENSMUSG00000002413",
    i_species="m_musculus",
    o_species="r_norvegicus",
    input_id="Ensembl_id",
    o_ids="Ensembl_id",
    mapping_tool="annotationDBI",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  m_musculus_ensembl_id <- "ENSMUSG00000002413"
  r_norvegicus_ensembl_id <- "ENSRNOG00000010957"
  expect_identical(getConversions(cosia_object)@converted_id,data.frame(m_musculus_ensembl_id,r_norvegicus_ensembl_id))

  cosia_object <- CoSIAn(
    gene_set = "ENSMUSG00000002413",
    i_species="m_musculus",
    o_species="d_rerio",
    input_id="Ensembl_id",
    o_ids="Ensembl_id",
    mapping_tool="annotationDBI",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  m_musculus_ensembl_id <- "ENSMUSG00000002413"
  d_rerio_ensembl_id <- "ENSDARG00000017661"
  expect_identical(getConversions(cosia_object)@converted_id,data.frame(m_musculus_ensembl_id,d_rerio_ensembl_id))

  cosia_object <- CoSIAn(
    gene_set = "ENSMUSG00000002413",
    i_species="m_musculus",
    o_species="c_elegans",
    input_id="Ensembl_id",
    o_ids="Ensembl_id",
    mapping_tool="annotationDBI",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  m_musculus_ensembl_id <- "ENSMUSG00000002413"
  c_elegans_ensembl_id <- "WBGene00003030"
  expect_identical(getConversions(cosia_object)@converted_id,data.frame(m_musculus_ensembl_id,c_elegans_ensembl_id))

  cosia_object <- CoSIAn(
    gene_set = "ENSMUSG00000002413",
    i_species="m_musculus",
    o_species="d_melanogaster",
    input_id="Ensembl_id",
    o_ids="Ensembl_id",
    mapping_tool="annotationDBI",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  m_musculus_ensembl_id <- "ENSMUSG00000002413"
  d_melanogaster_ensembl_id <- "FBgn0003079"
  expect_identical(getConversions(cosia_object)@converted_id,data.frame(m_musculus_ensembl_id,d_melanogaster_ensembl_id))


  #rat -> other species species
  cosia_object <- CoSIAn(
    gene_set = "ENSRNOG00000024972",
    i_species="r_norvegicus",
    o_species="h_sapiens",
    input_id="Ensembl_id",
    o_ids="Ensembl_id",
    mapping_tool="annotationDBI",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  r_norvegicus_ensembl_id <- "ENSRNOG00000024972"
  h_sapiens_ensembl_id <- "ENSG00000006695"
  expect_identical(getConversions(cosia_object)@converted_id,data.frame(r_norvegicus_ensembl_id,h_sapiens_ensembl_id))

  cosia_object <- CoSIAn(
    gene_set = "ENSRNOG00000024972",
    i_species="r_norvegicus",
    o_species="m_musculus",
    input_id="Ensembl_id",
    o_ids="Ensembl_id",
    mapping_tool="annotationDBI",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  r_norvegicus_ensembl_id <- "ENSRNOG00000024972"
  m_musculus_ensembl_id <- "ENSMUSG00000042148"
  expect_identical(getConversions(cosia_object)@converted_id,data.frame(r_norvegicus_ensembl_id,m_musculus_ensembl_id))

  cosia_object <- CoSIAn(
    gene_set = "ENSRNOG00000024972",
    i_species="r_norvegicus",
    o_species="r_norvegicus",
    input_id="Ensembl_id",
    o_ids="Ensembl_id",
    mapping_tool="annotationDBI",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  r_norvegicus_ensembl_id <- "ENSRNOG00000024972"
  expect_identical(getConversions(cosia_object)@converted_id,data.frame(r_norvegicus_ensembl_id))

  cosia_object <- CoSIAn(
    gene_set = "ENSRNOG00000024972",
    i_species="r_norvegicus",
    o_species="d_rerio",
    input_id="Ensembl_id",
    o_ids="Ensembl_id",
    mapping_tool="annotationDBI",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  r_norvegicus_ensembl_id <-c("ENSRNOG00000024972","ENSRNOG00000024972") #"ENSRNOG00000024972"
  d_rerio_ensembl_id <-c("ENSDARG00000034309","ENSDARG00000113352") #"ENSDARG00000034309"
  expect_identical(getConversions(cosia_object)@converted_id,data.frame(r_norvegicus_ensembl_id,d_rerio_ensembl_id))

  cosia_object <- CoSIAn(
    gene_set = "ENSRNOG00000024972",
    i_species="r_norvegicus",
    o_species="c_elegans",
    input_id="Ensembl_id",
    o_ids="Ensembl_id",
    mapping_tool="annotationDBI",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  r_norvegicus_ensembl_id <- "ENSRNOG00000024972"
  c_elegans_ensembl_id <- "WBGene00012895"
  expect_identical(getConversions(cosia_object)@converted_id,data.frame(r_norvegicus_ensembl_id,c_elegans_ensembl_id))

  cosia_object <- CoSIAn(
    gene_set = "ENSRNOG00000024972",
    i_species="r_norvegicus",
    o_species="d_melanogaster",
    input_id="Ensembl_id",
    o_ids="Ensembl_id",
    mapping_tool="annotationDBI",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  r_norvegicus_ensembl_id <- "ENSRNOG00000024972"
  d_melanogaster_ensembl_id <- "FBgn0032222"
  expect_identical(getConversions(cosia_object)@converted_id,data.frame(r_norvegicus_ensembl_id,d_melanogaster_ensembl_id))


  #d_rerio -> other species
  cosia_object <- CoSIAn(
    gene_set = "ENSDARG00000060006",
    i_species="d_rerio",
    o_species="h_sapiens",
    input_id="Ensembl_id",
    o_ids="Ensembl_id",
    mapping_tool="annotationDBI",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  d_rerio_ensembl_id <- "ENSDARG00000060006"
  h_sapiens_ensembl_id <- "ENSG00000049759"
  expect_identical(getConversions(cosia_object)@converted_id,data.frame(d_rerio_ensembl_id,h_sapiens_ensembl_id))

  cosia_object <- CoSIAn(
    gene_set = "ENSDARG00000060006",
    i_species="d_rerio",
    o_species="m_musculus",
    input_id="Ensembl_id",
    o_ids="Ensembl_id",
    mapping_tool="annotationDBI",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  d_rerio_ensembl_id <- "ENSDARG00000060006"
  m_musculus_ensembl_id <- "ENSMUSG00000024589"
  expect_identical(getConversions(cosia_object)@converted_id,data.frame(d_rerio_ensembl_id,m_musculus_ensembl_id))

  cosia_object <- CoSIAn(
    gene_set = "ENSDARG00000060006",
    i_species="d_rerio",
    o_species="r_norvegicus",
    input_id="Ensembl_id",
    o_ids="Ensembl_id",
    mapping_tool="annotationDBI",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  d_rerio_ensembl_id <- "ENSDARG00000060006"
  r_norvegicus_ensembl_id <- "ENSRNOG00000017610"
  expect_identical(getConversions(cosia_object)@converted_id,data.frame(d_rerio_ensembl_id,r_norvegicus_ensembl_id))


  cosia_object <- CoSIAn(
    gene_set = "ENSDARG00000060006",
    i_species="d_rerio",
    o_species="d_rerio",
    input_id="Ensembl_id",
    o_ids="Ensembl_id",
    mapping_tool="annotationDBI",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  d_rerio_ensembl_id <- "ENSDARG00000060006"
  expect_identical(getConversions(cosia_object)@converted_id,data.frame(d_rerio_ensembl_id))

  cosia_object <- CoSIAn(
    gene_set = "ENSDARG00000060006",
    i_species="d_rerio",
    o_species="c_elegans",
    input_id="Ensembl_id",
    o_ids="Ensembl_id",
    mapping_tool="annotationDBI",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  d_rerio_ensembl_id <- "ENSDARG00000060006"
  c_elegans_ensembl_id <- "WBGene00022358"
  expect_identical(getConversions(cosia_object)@converted_id,data.frame(d_rerio_ensembl_id,c_elegans_ensembl_id))

  cosia_object <- CoSIAn(
    gene_set = "ENSDARG00000060006",
    i_species="d_rerio",
    o_species="d_melanogaster",
    input_id="Ensembl_id",
    o_ids="Ensembl_id",
    mapping_tool="annotationDBI",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  d_rerio_ensembl_id <- "ENSDARG00000060006"
  d_melanogaster_ensembl_id <- "FBgn0259174"
  expect_identical(getConversions(cosia_object)@converted_id,data.frame(d_rerio_ensembl_id,d_melanogaster_ensembl_id ))


  #c_elegans -> other species
  cosia_object <- CoSIAn(
    gene_set = "WBGene00000838",
    i_species="c_elegans",
    o_species="h_sapiens",
    input_id="Ensembl_id",
    o_ids="Ensembl_id",
    mapping_tool="annotationDBI",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  c_elegans_ensembl_id <- "WBGene00000838"
  h_sapiens_ensembl_id <- "ENSG00000036257"
  expect_identical(getConversions(cosia_object)@converted_id,data.frame(c_elegans_ensembl_id,h_sapiens_ensembl_id))

  cosia_object <- CoSIAn(
    gene_set = "WBGene00000838",
    i_species="c_elegans",
    o_species="m_musculus",
    input_id="Ensembl_id",
    o_ids="Ensembl_id",
    mapping_tool="annotationDBI",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  c_elegans_ensembl_id <- "WBGene00000838"
  m_musculus_ensembl_id <- "ENSMUSG00000004364"
  expect_identical(getConversions(cosia_object)@converted_id,data.frame(c_elegans_ensembl_id,m_musculus_ensembl_id))

  cosia_object <- CoSIAn(
    gene_set = "WBGene00000838",
    i_species="c_elegans",
    o_species="r_norvegicus",
    input_id="Ensembl_id",
    o_ids="Ensembl_id",
    mapping_tool="annotationDBI",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  c_elegans_ensembl_id <- "WBGene00000838"
  r_norvegicus_ensembl_id <- "ENSRNOG00000015633"
  expect_identical(getConversions(cosia_object)@converted_id,data.frame(c_elegans_ensembl_id,r_norvegicus_ensembl_id))

  cosia_object <- CoSIAn(
    gene_set = "WBGene00000838",
    i_species="c_elegans",
    o_species="d_rerio",
    input_id="Ensembl_id",
    o_ids="Ensembl_id",
    mapping_tool="annotationDBI",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  c_elegans_ensembl_id <- "WBGene00000838"
  d_rerio_ensembl_id <- "ENSDARG00000100108"
  expect_identical(getConversions(cosia_object)@converted_id,data.frame(c_elegans_ensembl_id,d_rerio_ensembl_id))

  cosia_object <- CoSIAn(
    gene_set = "WBGene00000838",
    i_species="c_elegans",
    o_species="c_elegans",
    input_id="Ensembl_id",
    o_ids="Ensembl_id",
    mapping_tool="annotationDBI",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  c_elegans_ensembl_id <- "WBGene00000838"
  expect_identical(getConversions(cosia_object)@converted_id,data.frame(c_elegans_ensembl_id))

  cosia_object <- CoSIAn(
    gene_set = "WBGene00000838",
    i_species="c_elegans",
    o_species="d_melanogaster",
    input_id="Ensembl_id",
    o_ids="Ensembl_id",
    mapping_tool="annotationDBI",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  c_elegans_ensembl_id <- "WBGene00000838"
  d_melanogaster_ensembl_id <- "FBgn0261268"
  expect_identical(getConversions(cosia_object)@converted_id,data.frame(c_elegans_ensembl_id,d_melanogaster_ensembl_id))


  #d melanogaster -> other species
  cosia_object <- CoSIAn(
    gene_set = "FBgn0261268",
    i_species="d_melanogaster",
    o_species="h_sapiens",
    input_id="Ensembl_id",
    o_ids="Ensembl_id",
    mapping_tool="annotationDBI",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  d_melanogaster_ensembl_id <- "FBgn0261268"
  h_sapiens_ensembl_id <- "ENSG00000036257"
  expect_identical(getConversions(cosia_object)@converted_id,data.frame(d_melanogaster_ensembl_id,h_sapiens_ensembl_id))

  cosia_object <- CoSIAn(
    gene_set = "FBgn0261268",
    i_species="d_melanogaster",
    o_species="m_musculus",
    input_id="Ensembl_id",
    o_ids="Ensembl_id",
    mapping_tool="annotationDBI",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  d_melanogaster_ensembl_id <- "FBgn0261268"
  m_musculus_ensembl_id <- "ENSMUSG00000004364"
  expect_identical(getConversions(cosia_object)@converted_id,data.frame(d_melanogaster_ensembl_id,m_musculus_ensembl_id))

  cosia_object <- CoSIAn(
    gene_set = "FBgn0261268",
    i_species="d_melanogaster",
    o_species="r_norvegicus",
    input_id="Ensembl_id",
    o_ids="Ensembl_id",
    mapping_tool="annotationDBI",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  d_melanogaster_ensembl_id <- "FBgn0261268"
  r_norvegicus_ensembl_id <- "ENSRNOG00000015633"
  expect_identical(getConversions(cosia_object)@converted_id,data.frame(d_melanogaster_ensembl_id,r_norvegicus_ensembl_id))

  cosia_object <- CoSIAn(
    gene_set = "FBgn0261268",
    i_species="d_melanogaster",
    o_species="d_rerio",
    input_id="Ensembl_id",
    o_ids="Ensembl_id",
    mapping_tool="annotationDBI",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  d_melanogaster_ensembl_id <- "FBgn0261268"
  d_rerio_ensembl_id <- "ENSDARG00000100108"
  expect_identical(getConversions(cosia_object)@converted_id,data.frame(d_melanogaster_ensembl_id,d_rerio_ensembl_id))

  cosia_object <- CoSIAn(
    gene_set = "FBgn0261268",
    i_species="d_melanogaster",
    o_species="c_elegans",
    input_id="Ensembl_id",
    o_ids="Ensembl_id",
    mapping_tool="annotationDBI",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  d_melanogaster_ensembl_id <- "FBgn0261268"
  c_elegans_ensembl_id <- "WBGene00000838"
  expect_identical(getConversions(cosia_object)@converted_id,data.frame(d_melanogaster_ensembl_id,c_elegans_ensembl_id))

  cosia_object <- CoSIAn(
    gene_set = "FBgn0261268",
    i_species="d_melanogaster",
    o_species="d_melanogaster",
    input_id="Ensembl_id",
    o_ids="Ensembl_id",
    mapping_tool="annotationDBI",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  d_melanogaster_ensembl_id <- "FBgn0261268"
  expect_identical(getConversions(cosia_object)@converted_id,data.frame(d_melanogaster_ensembl_id))

})

##No Tool error
test_that("no tool",{
  cosia_object <- CoSIAn(
    gene_set = "ENSG00000040531",
    i_species="h_sapiens",
    o_species="h_sapiens",
    input_id="Ensembl_id",
    o_ids="Ensembl_id",
    mapping_tool="FALSE_TOOL",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  expect_error(getConversions(cosia_object)@converted_id,"Error: Invalid tool in CoSIAn Object. Make sure the tool is either annotationDBI or biomaRt.")

  cosia_object <- CoSIAn(
    gene_set = "ENSG00000040531",
    i_species="m_musculus",
    o_species="h_sapiens",
    input_id="Ensembl_id",
    o_ids="Ensembl_id",
    mapping_tool="FALSE_TOOL",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  expect_error(getConversions(cosia_object)@converted_id,"Error: Invalid tool in CoSIAn Object. Make sure the tool is either annotationDBI or biomaRt.")

  cosia_object <- CoSIAn(
    gene_set = "ENSG00000040531",
    i_species="r_norvegicus",
    o_species="h_sapiens",
    input_id="Ensembl_id",
    o_ids="Ensembl_id",
    mapping_tool="FALSE_TOOL",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  expect_error(getConversions(cosia_object)@converted_id,"Error: Invalid tool in CoSIAn Object. Make sure the tool is either annotationDBI or biomaRt.")

  cosia_object <- CoSIAn(
    gene_set = "ENSG00000040531",
    i_species="d_rerio",
    o_species="h_sapiens",
    input_id="Ensembl_id",
    o_ids="Ensembl_id",
    mapping_tool="FALSE_TOOL",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  expect_error(getConversions(cosia_object)@converted_id,"Error: Invalid tool in CoSIAn Object. Make sure the tool is either annotationDBI or biomaRt.")

  cosia_object <- CoSIAn(
    gene_set = "ENSG00000040531",
    i_species="d_rerio",
    o_species="h_sapiens",
    input_id="Ensembl_id",
    o_ids="Ensembl_id",
    mapping_tool="FALSE_TOOL",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  expect_error(getConversions(cosia_object)@converted_id,"Error: Invalid tool in CoSIAn Object. Make sure the tool is either annotationDBI or biomaRt.")

  cosia_object <- CoSIAn(
    gene_set = "ENSG00000040531",
    i_species="c_elegans",
    o_species="h_sapiens",
    input_id="Ensembl_id",
    o_ids="Ensembl_id",
    mapping_tool="FALSE_TOOL",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  expect_error(getConversions(cosia_object)@converted_id,"Error: Invalid tool in CoSIAn Object. Make sure the tool is either annotationDBI or biomaRt.")

  cosia_object <- CoSIAn(
    gene_set = "ENSG00000040531",
    i_species="d_melanogaster",
    o_species="h_sapiens",
    input_id="Ensembl_id",
    o_ids="Ensembl_id",
    mapping_tool="FALSE_TOOL",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  expect_error(getConversions(cosia_object)@converted_id,"Error: Invalid tool in CoSIAn Object. Make sure the tool is either annotationDBI or biomaRt.")
})

#across id_types
test_that("across id types",{
  cosia_object <- CoSIAn(
    gene_set = "ENSG00000040531",
    i_species="h_sapiens",
    o_species="h_sapiens",
    input_id="Ensembl_id",
    o_ids=c("Entrez_id","Ensembl_id","Gene_name","Symbol"),
    mapping_tool="annotationDBI",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  h_sapiens_ensembl_id <- "ENSG00000040531"
  h_sapiens_entrez_id <- "1497"
  h_sapiens_gene_name <- "cystinosin, lysosomal cystine transporter"
  h_sapiens_symbol <- "CTNS"
  expect_identical(getConversions(cosia_object)@converted_id,data.frame(h_sapiens_ensembl_id,h_sapiens_entrez_id,h_sapiens_gene_name,h_sapiens_symbol))
 

  cosia_object <- CoSIAn(
    gene_set = "ENSG00000040531",
    i_species="h_sapiens",
    o_species="h_sapiens",
    input_id="Ensembl_id",
    o_ids=c("Ensembl_id","Entrez_id","Ensembl_id_version","Gene_name","HGNC_Symbol"),
    mapping_tool="biomaRt",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  h_sapiens_ensembl_id <- "ENSG00000040531"
  h_sapiens_entrez_id <- "1497"
  h_sapiens_ensembl_id_version <- "ENSG00000040531.16"
  h_sapiens_gene_name <- "CTNS"
  h_sapiens_hgnc_symbol <- "CTNS"
  expect_identical(getConversions(cosia_object)@converted_id,data.frame(h_sapiens_ensembl_id,h_sapiens_entrez_id,h_sapiens_ensembl_id_version,h_sapiens_gene_name,h_sapiens_hgnc_symbol))
 
})

test_that("NCBI",{
  cosia_object <- CoSIAn(
    gene_set = "ENSG00000040531",
    i_species="h_sapiens",
    o_species="m_musculus",
    input_id="Ensembl_id",
    o_ids="Ensembl_id",
    mapping_tool="annotationDBI",
    ortholog_database = "NCBIOrtho",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  h_sapiens_ensembl_id <- "ENSG00000040531"
  m_musculus_ensembl_id <- "ENSMUSG00000005949"
  expect_identical(getConversions(cosia_object)@converted_id,data.frame(h_sapiens_ensembl_id,m_musculus_ensembl_id))
})

#version numbers annotationDBI
test_that("version_numbers",{
  cosia_object <- CoSIAn(
    gene_set = "ENSG00000040531.16",
    i_species="h_sapiens",
    o_species="h_sapiens",
    input_id="Ensembl_id_version",
    o_ids=c("Ensembl_id","Entrez_id"),
    mapping_tool="annotationDBI",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  h_sapiens_entrez_id <- "1497"
  h_sapiens_ensembl_id <- "ENSG00000040531"
  h_sapiens_ensembl_id_version <- "ENSG00000040531.16"
  expect_identical(getConversions(cosia_object)@converted_id,data.frame(h_sapiens_ensembl_id_version,h_sapiens_ensembl_id,h_sapiens_entrez_id))

  cosia_object <- CoSIAn(
    gene_set = "ENSMUSG00000005949.10",
    i_species="m_musculus",
    o_species="m_musculus",
    input_id="Ensembl_id_version",
    o_ids="Ensembl_id",
    mapping_tool="annotationDBI",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  m_musculus_ensembl_id <- "ENSMUSG00000005949"
  m_musculus_ensembl_id_version <- "ENSMUSG00000005949.10"
  expect_identical(getConversions(cosia_object)@converted_id,data.frame(m_musculus_ensembl_id_version,m_musculus_ensembl_id))

  cosia_object <- CoSIAn(
    gene_set = "ENSRNOG00000028688.8",
    i_species="r_norvegicus",
    o_species="r_norvegicus",
    input_id="Ensembl_id_version",
    o_ids="Ensembl_id",
    mapping_tool="annotationDBI",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  r_norvegicus_ensembl_id <- "ENSRNOG00000028688"
  r_norvegicus_ensembl_id_version <- "ENSRNOG00000028688.8"
  expect_identical(getConversions(cosia_object)@converted_id,data.frame(r_norvegicus_ensembl_id_version,r_norvegicus_ensembl_id))

  cosia_object <- CoSIAn(
    gene_set = "ENSDARG00000008890.8",
    i_species="d_rerio",
    o_species="d_rerio",
    input_id="Ensembl_id_version",
    o_ids="Ensembl_id",
    mapping_tool="annotationDBI",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  d_rerio_ensembl_id <- "ENSDARG00000008890"
  d_rerio_ensembl_id_version <- "ENSDARG00000008890.8"
  expect_identical(getConversions(cosia_object)@converted_id,data.frame(d_rerio_ensembl_id_version,d_rerio_ensembl_id))

  cosia_object <- CoSIAn(
    gene_set = "ENSDARG00000008890.8",
    i_species="c_elegans",
    o_species="d_rerio",
    input_id="Ensembl_id_version",
    o_ids="Ensembl_id",
    mapping_tool="annotationDBI",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  expect_error(getConversions(cosia_object),"Error. Incorrect ID evaluation. C. elegans do not have Ensembl Ids with Version.")

  cosia_object <- CoSIAn(
    gene_set = "ENSDARG00000008890.8",
    i_species="d_melanogaster",
    o_species="d_rerio",
    input_id="Ensembl_id_version",
    o_ids="Ensembl_id",
    mapping_tool="annotationDBI",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  expect_error(getConversions(cosia_object),"Error. Incorrect ID evaluation. D. melanogaster do not have Ensembl Ids with Version.")

  cosia_object <- CoSIAn(
    gene_set = "ENSG00000040531.16",
    i_species="h_sapiens",
    o_species="m_musculus",
    input_id="Ensembl_id_version",
    o_ids=c("Ensembl_id","Entrez_id"),
    mapping_tool="annotationDBI",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  m_musculus_entrez_id <- "83429"
  m_musculus_ensembl_id <- "ENSMUSG00000005949"
  h_sapiens_ensembl_id_version <- "ENSG00000040531.16"
  expect_identical(getConversions(cosia_object)@converted_id,data.frame(h_sapiens_ensembl_id_version,m_musculus_entrez_id,m_musculus_ensembl_id))
 

  cosia_object <- CoSIAn(
    gene_set = "ENSG00000040531.16",
    i_species="h_sapiens",
    o_species="m_musculus",
    input_id="Ensembl_id_version",
    o_ids=c("Ensembl_id"),
    mapping_tool="annotationDBI",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  m_musculus_ensembl_id <- "ENSMUSG00000005949"
  h_sapiens_ensembl_id_version <- "ENSG00000040531.16"
  expect_identical(getConversions(cosia_object)@converted_id,data.frame(h_sapiens_ensembl_id_version,m_musculus_ensembl_id))
 
})

##Biomart (about 30%)
test_that("Check biomaRt",{
  #so that this section isn't skipped
  expect_identical(1,1)

  cosia_object <- CoSIAn(
    gene_set = "ENSG00000040531",
    i_species="h_sapiens",
    o_species=c("h_sapiens","m_musculus","r_norvegicus","d_rerio","c_elegans","d_melanogaster"),
    input_id="Ensembl_id",
    o_ids="Ensembl_id",
    mapping_tool="biomaRt",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  h_sapiens_ensembl_id <- "ENSG00000040531"
  m_musculus_ensembl_id <- "ENSMUSG00000005949"
  r_norvegicus_ensembl_id <- "ENSRNOG00000028688"
  d_rerio_ensembl_id <- "ENSDARG00000008890"
  c_elegans_ensembl_id <- "WBGene00008052"
  d_melanogaster_ensembl_id <- "FBgn0039045"
  expect_identical(getConversions(cosia_object)@converted_id,data.frame(h_sapiens_ensembl_id,m_musculus_ensembl_id,r_norvegicus_ensembl_id,d_rerio_ensembl_id,c_elegans_ensembl_id,d_melanogaster_ensembl_id))

  cosia_object <- CoSIAn(
    gene_set = "ENSMUSG00000005949",
    i_species="m_musculus",
    o_species=c("h_sapiens","m_musculus","r_norvegicus","d_rerio","c_elegans","d_melanogaster"),
    input_id="Ensembl_id",
    o_ids="Ensembl_id",
    mapping_tool="biomaRt",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  h_sapiens_ensembl_id <- "ENSG00000040531"
  m_musculus_ensembl_id <- "ENSMUSG00000005949"
  r_norvegicus_ensembl_id <- "ENSRNOG00000028688"
  d_rerio_ensembl_id <- "ENSDARG00000008890"
  c_elegans_ensembl_id <- "WBGene00008052"
  d_melanogaster_ensembl_id <- "FBgn0039045"
  expect_identical(getConversions(cosia_object)@converted_id,data.frame(m_musculus_ensembl_id,h_sapiens_ensembl_id,r_norvegicus_ensembl_id,d_rerio_ensembl_id,c_elegans_ensembl_id,d_melanogaster_ensembl_id))

  cosia_object <- CoSIAn(
    gene_set = "ENSRNOG00000028688",
    i_species="r_norvegicus",
    o_species=c("h_sapiens","m_musculus","r_norvegicus","d_rerio","c_elegans","d_melanogaster"),
    input_id="Ensembl_id",
    o_ids="Ensembl_id",
    mapping_tool="biomaRt",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  h_sapiens_ensembl_id <- "ENSG00000040531"
  m_musculus_ensembl_id <- "ENSMUSG00000005949"
  r_norvegicus_ensembl_id <- "ENSRNOG00000028688"
  d_rerio_ensembl_id <- "ENSDARG00000008890"
  c_elegans_ensembl_id <- "WBGene00008052"
  d_melanogaster_ensembl_id <- "FBgn0039045"
  expect_identical(getConversions(cosia_object)@converted_id,data.frame(r_norvegicus_ensembl_id,h_sapiens_ensembl_id,m_musculus_ensembl_id,d_rerio_ensembl_id,c_elegans_ensembl_id,d_melanogaster_ensembl_id))

  cosia_object <- CoSIAn(
    gene_set = "ENSDARG00000008890",
    i_species="d_rerio",
    o_species=c("h_sapiens","m_musculus","r_norvegicus","d_rerio","c_elegans","d_melanogaster"),
    input_id="Ensembl_id",
    o_ids="Ensembl_id",
    mapping_tool="biomaRt",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  h_sapiens_ensembl_id <- "ENSG00000040531"
  m_musculus_ensembl_id <- "ENSMUSG00000005949"
  r_norvegicus_ensembl_id <- "ENSRNOG00000028688"
  d_rerio_ensembl_id <- "ENSDARG00000008890"
  c_elegans_ensembl_id <- "WBGene00008052"
  d_melanogaster_ensembl_id <- "FBgn0039045"
  expect_identical(getConversions(cosia_object)@converted_id,data.frame(d_rerio_ensembl_id,h_sapiens_ensembl_id,m_musculus_ensembl_id,r_norvegicus_ensembl_id,c_elegans_ensembl_id,d_melanogaster_ensembl_id))

  cosia_object <- CoSIAn(
    gene_set = "WBGene00008052",
    i_species="c_elegans",
    o_species=c("h_sapiens","m_musculus","r_norvegicus","d_rerio","c_elegans","d_melanogaster"),
    input_id="Ensembl_id",
    o_ids="Ensembl_id",
    mapping_tool="biomaRt",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  h_sapiens_ensembl_id <- "ENSG00000040531"
  m_musculus_ensembl_id <- "ENSMUSG00000005949"
  r_norvegicus_ensembl_id <- "ENSRNOG00000028688"
  d_rerio_ensembl_id <- "ENSDARG00000008890"
  c_elegans_ensembl_id <- "WBGene00008052"
  d_melanogaster_ensembl_id <- "FBgn0039045"
  expect_identical(getConversions(cosia_object)@converted_id,data.frame(c_elegans_ensembl_id,h_sapiens_ensembl_id,m_musculus_ensembl_id,r_norvegicus_ensembl_id,d_rerio_ensembl_id,d_melanogaster_ensembl_id))

  cosia_object <- CoSIAn(
    gene_set = "FBgn0039045",
    i_species="d_melanogaster",
    o_species=c("h_sapiens","m_musculus","r_norvegicus","d_rerio","c_elegans","d_melanogaster"),
    input_id="Ensembl_id",
    o_ids="Ensembl_id",
    mapping_tool="biomaRt",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  h_sapiens_ensembl_id <- "ENSG00000040531"
  m_musculus_ensembl_id <- "ENSMUSG00000005949"
  r_norvegicus_ensembl_id <- "ENSRNOG00000028688"
  d_rerio_ensembl_id <- "ENSDARG00000008890"
  c_elegans_ensembl_id <- "WBGene00008052"
  d_melanogaster_ensembl_id <- "FBgn0039045"
  expect_identical(getConversions(cosia_object)@converted_id,data.frame(d_melanogaster_ensembl_id,h_sapiens_ensembl_id,m_musculus_ensembl_id,r_norvegicus_ensembl_id,d_rerio_ensembl_id,c_elegans_ensembl_id))


})

test_that("edge_cases",{
  cosia_object <- CoSIAn(
    gene_set = c("ENSG00000008710","ENSG00000118762","ENSG00000152217"),
    i_species="h_sapiens",
    o_species=c("h_sapiens","m_musculus","r_norvegicus","d_rerio","c_elegans","d_melanogaster"),
    input_id="Ensembl_id",
    o_ids=c("Ensembl_id","Gene_name","Entrez_id"),
    mapping_tool="biomaRt",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  getConversions(cosia_object)@converted_id
  
  cosia_object <- CoSIAn(
    gene_set = c("ENSG00000008710","ENSG00000118762","ENSG00000152217"),
    i_species="h_sapiens",
    o_species=c("h_sapiens","m_musculus","r_norvegicus","d_rerio","c_elegans","d_melanogaster"),
    input_id="Ensembl_id",
    o_ids=c("Ensembl_id","Gene_name","Entrez_id", "Symbol"),
    mapping_tool="annotationDBI",
    ortholog_database = "HomoloGene",
    map_tissues = "",
    map_species = "",
    metric_type=""
  )
  getConversions(cosia_object)@converted_id
})