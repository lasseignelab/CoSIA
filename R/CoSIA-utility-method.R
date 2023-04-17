# CoSIA Utility Functions

# Find the Available/Common Tissues in CoSIA for a given combination of species

#' getTissues
#'
#' @param species name of a species or multiple species that you want to get
#' available tissue list for
#'
#' @return list of tissues that are common/available among the species or
#' multiple species inputted
#' @export
#'
#' @examples
#' tissue <- getTissues(c("m_musculus"))
getTissues <- function(species) {
    # Loading CoSIAdata
    CoSIAdata_load <- function(species) {
        eh <- ExperimentHub::ExperimentHub()
        merged_CoSIAdata <- data.frame(matrix(ncol = 7, nrow = 0))
        colnames(merged_CoSIAdata) <- c(
            "Anatomical_entity_name", "Ensembl_ID",
            "Sample_size", "VST", "Experiment_ID",
            "Anatomical_entity_ID", "Species"
        )
        if (any(species == "m_musculus")) {
            mm_EH_File <- eh[["EH7859"]]
            merged_CoSIAdata <- rbind(merged_CoSIAdata, mm_EH_File)
            merged_CoSIAdata <- as.data.frame(merged_CoSIAdata)
        } else if (any(species == "r_norvegicus")) {
            rn_EH_File <- eh[["EH7860"]]
            merged_CoSIAdata <- rbind(merged_CoSIAdata, rn_EH_File)
            merged_CoSIAdata <- as.data.frame(merged_CoSIAdata)
        } else if (any(species == "d_rerio")) {
            dr_EH_File <- eh[["EH7861"]]
            merged_CoSIAdata <- rbind(merged_CoSIAdata, dr_EH_File)
            merged_CoSIAdata <- as.data.frame(merged_CoSIAdata)
        } else if (any(species == "h_sapiens")) {
            hs_EH_File <- eh[["EH7858"]]
            merged_CoSIAdata <- rbind(merged_CoSIAdata, hs_EH_File)
            merged_CoSIAdata <- as.data.frame(merged_CoSIAdata)
        } else if (any(species == "c_elegans")) {
            ce_EH_File <- eh[["EH7863"]]
            merged_CoSIAdata <- rbind(merged_CoSIAdata, ce_EH_File)
            merged_CoSIAdata <- as.data.frame(merged_CoSIAdata)
        } else if (any(species == "d_melanogaster")) {
            dm_EH_File <- eh[["EH7862"]]
            merged_CoSIAdata <- rbind(merged_CoSIAdata, dm_EH_File)
            merged_CoSIAdata <- as.data.frame(merged_CoSIAdata)
        } else {
            stop("Issue with species in CoSIAn Object. Make sure the species in the
          species argument are avalible organisms through CoSIA and are in the
          correct format.")
        }
        return(merged_CoSIAdata)
    }

    merged_CoSIAdata <- lapply(species, CoSIAdata_load)
    merged_CoSIAdata <- as.data.frame(do.call(rbind, merged_CoSIAdata))

    List_of_Tissues <- merged_CoSIAdata %>%
        dplyr::group_by(Anatomical_entity_name) %>%
        dplyr::reframe(
            Anatomical_entity_ID = unique(Anatomical_entity_ID),
            Species = unique(Species)
        )
    L <- List_of_Tissues %>% dplyr::group_by(Anatomical_entity_name)%>%
        dplyr::reframe(Frequency = table(Anatomical_entity_name))
    frequency_value <- length(species)
    common_tissue <- dplyr::filter(L, Frequency == frequency_value)
    common_tissue <- subset(common_tissue, select = -c(Frequency))
    colnames(common_tissue)[which(names(common_tissue) ==
        "Anatomical_entity_name")] <- "Common_Anatomical_Entity_Name"
    return(common_tissue)
}

#CoSIA Accessor Functions

#' viewConvertId Generics
#'
#' @param object CoSIAn object with all user accessible slots filled
#'
#' @return initializes a generic function for getConvertId as preparation for
#' defining the getConvertId Method
#' @export
#'
#' @examples
#' Kidney_Genes <- CoSIAn(
#'     gene_set = c("ENSG00000008710", "ENSG00000118762", "ENSG00000152217"),
#'     i_species = "h_sapiens", input_id = "Ensembl_id",
#'     o_species = c(
#'         "d_melanogaster", "m_musculus", "h_sapiens", "d_rerio",
#'         "c_elegans", "r_norvegicus"
#'     ),
#'     output_ids = c("Ensembl_id", "Symbol"), mapping_tool = "annotationDBI",
#'     ortholog_database = "HomoloGene", map_tissues = "heart",
#'     map_species = c("m_musculus"), metric_type = "DS_Gene"
#' )
#' Kidney_gene_conversion <- CoSIA::getConversion(Kidney_Genes)
#' viewConvertId(Kidney_gene_conversion)
setGeneric(
  "viewConvertId", 
  function(object) standardGeneric("viewConvertId"))

#' viewConvertId
#'
#' @param object CoSIAn object with all user accessible slots filled
#' @return converted_id slot in CoSIAn Object
#' @export
#'
#' @examples
#' Kidney_Genes <- CoSIAn(
#'     gene_set = c("ENSG00000008710", "ENSG00000118762", "ENSG00000152217"),
#'     i_species = "h_sapiens", input_id = "Ensembl_id",
#'     o_species = c(
#'         "d_melanogaster", "m_musculus", "h_sapiens", "d_rerio",
#'         "c_elegans", "r_norvegicus"
#'     ),
#'     output_ids = c("Ensembl_id", "Symbol"), mapping_tool = "annotationDBI",
#'     ortholog_database = "HomoloGene", map_tissues = "heart",
#'     map_species = c("m_musculus"), metric_type = "DS_Gene"
#' )
#' Kidney_gene_conversion <- CoSIA::getConversion(Kidney_Genes)
#' viewConvertId(Kidney_gene_conversion)
setMethod(
  "viewConvertId", 
  signature(object = "CoSIAn"), 
  function(object) object@converted_id)

#' viewGEx Generic
#'
#' @param object CoSIAn object with all user accessible slots filled with
#' converted_id slot filled
#'
#' @return initializes a generic function for viewGEx as preparation for
#' defining the viewGEx Method
#' @export
#'
#' @examples
#' Kidney_Genes <- CoSIAn(
#'     gene_set = c("ENSG00000008710", "ENSG00000118762", "ENSG00000152217"),
#'     i_species = "h_sapiens", input_id = "Ensembl_id",
#'     o_species = c(
#'         "d_melanogaster", "m_musculus",
#'         "h_sapiens", "d_rerio", "c_elegans", "r_norvegicus"
#'     ), output_ids = c("Ensembl_id", "Symbol"),
#'     mapping_tool = "annotationDBI", ortholog_database = "HomoloGene",
#'     map_tissues = "heart", map_species = c("m_musculus"),
#'     metric_type = "DS_Gene"
#' )
#' Kidney_gene_conversion <- CoSIA::getConversion(Kidney_Genes)
#' Kidney_gene_gex <- getGEx(Kidney_gene_conversion)
#' viewGEx(Kidney_gene_gex)
setGeneric(
  "viewGEx",
  function(object) standardGeneric("viewGEx"))
#' viewGEx
#'
#' @param object CoSIAn object with all user accessible slots filled with
#' converted_id slot filled
#'
#' @return gex slot in CoSIAn Object
#' @export
#'
#' @examples
#' Kidney_Genes <- CoSIAn(
#'     gene_set = c("ENSG00000008710", "ENSG00000118762", "ENSG00000152217"),
#'     i_species = "h_sapiens", input_id = "Ensembl_id",
#'     o_species = c(
#'         "d_melanogaster", "m_musculus",
#'         "h_sapiens", "d_rerio", "c_elegans", "r_norvegicus"
#'     ), output_ids = c("Ensembl_id", "Symbol"),
#'     mapping_tool = "annotationDBI", ortholog_database = "HomoloGene",
#'     map_tissues = "heart", map_species = c("m_musculus"),
#'     metric_type = "DS_Gene"
#' )
#' Kidney_gene_conversion <- CoSIA::getConversion(Kidney_Genes)
#' Kidney_gene_gex <- getGEx(Kidney_gene_conversion)
#' viewGEx(Kidney_gene_gex)
setMethod( "viewGEx",
           signature(object = "CoSIAn"),
           function(object) object@gex)