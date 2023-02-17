#' getGExMetrics Generic
#'
#' @param object CoSIAn object with all user accessible slots filled with converted_id slot filled
#'
#' @return initializes a generic function for getGExMetrics as preparation for defining the getGExMetrics Method
#' @export
#' @examples
#' Kidney_Genes<-CoSIAn(gene_set = c('ENSG00000008710','ENSG00000118762','ENSG00000152217'),
#' i_species = 'h_sapiens',input_id = 'Ensembl_id',o_species = c('d_melanogaster','m_musculus',
#' 'h_sapiens', 'd_rerio','c_elegans','r_norvegicus'),output_ids = c('Ensembl_id','Symbol'), 
#' mapping_tool = 'annotationDBI',ortholog_database = 'HomoloGene',map_tissues = 'heart', 
#' map_species = c('m_musculus'),metric_type = 'DS_Gene')
#' Kidney_gene_conversion<-CoSIA::getConversion(Kidney_Genes)
#' load('~/Desktop/EH_Data.RData')
#' Kidney_gene_metric<-getGExMetrics(Kidney_gene_conversion)

setGeneric("getGExMetrics", function(object) standardGeneric("getGExMetrics"))
# CoSIAn getGEx


#' getGExMetrics Method
#'
#' @param object CoSIAn object with all user accessible slots filled with converted_id slot filled 
#'
#' @return CoSIAn Object with metric slot filled
#' @export
#' @importFrom dplyr group_by summarise
#' @importFrom tibble remove_rownames column_to_rownames
#'
#' @examples
#' Kidney_Genes<-CoSIAn(gene_set = c('ENSG00000008710','ENSG00000118762','ENSG00000152217'),
#' i_species = 'h_sapiens',input_id = 'Ensembl_id',o_species = c('d_melanogaster','m_musculus',
#' 'h_sapiens', 'd_rerio','c_elegans','r_norvegicus'),output_ids = c('Ensembl_id','Symbol'), 
#' mapping_tool = 'annotationDBI',ortholog_database = 'HomoloGene',map_tissues = 'heart', 
#' map_species = c('m_musculus'),metric_type = 'DS_Gene')
#' Kidney_gene_conversion<-CoSIA::getConversion(Kidney_Genes)
#' load('~/Desktop/EH_Data.RData')
#' Kidney_gene_metric<-getGExMetrics(Kidney_gene_conversion)

setMethod("getGExMetrics", signature(object = "CoSIAn"), function(object) {
    id_dataframe <- object@converted_id
    id_dataframe <- as.data.frame(id_dataframe)
    map_species <- object@map_species
    map_species <- as.character(map_species)
    Converted_Species_SWITCH <- Vectorize(
      vectorize.args = "map_species", 
      FUN = function(map_species) {
        switch(as.character(map_species), 
               h_sapiens = "h_sapiens_ensembl_id", 
               m_musculus = "m_musculus_ensembl_id", 
               r_norvegicus = "r_norvegicus_ensembl_id",
               d_rerio = "d_rerio_ensembl_id", 
               d_melanogaster = "d_melanogaster_ensembl_id", 
               c_elegans = "c_elegans_ensembl_id",
            stop("Error: Invalid map_species in CoSIAn Object. Make sure the 
            species in the map_species slot are an avalible model organism and 
            are in the correct format."))
    })
    map_species <- Converted_Species_SWITCH(map_species)
    id_dataframe <- dplyr::select(id_dataframe, tidyselect::matches(map_species))
    colnames(id_dataframe)[which(names(id_dataframe) == "h_sapiens")] <- "h_sapiens_ensembl_id"
    colnames(id_dataframe)[which(names(id_dataframe) == "m_musculus")] <- "m_musculus_ensembl_id"
    colnames(id_dataframe)[which(names(id_dataframe) == "r_norvegicus")] <- "r_norvegicus_ensembl_id"
    colnames(id_dataframe)[which(names(id_dataframe) == "d_rerio")] <- "d_rerio_ensembl_id"
    colnames(id_dataframe)[which(names(id_dataframe) == "d_melanogaster")] <- "d_melanogaster_ensembl_id"
    colnames(id_dataframe)[which(names(id_dataframe) == "c_elegans")] <- "c_elegans_ensembl_id"
    map_tissues <- object@map_tissues
    map_tissues <- as.character(map_tissues)
    Species_SWITCH <- Vectorize(
      vectorize.args = "map_species", 
      FUN = function(map_species) {
        switch(as.character(map_species), 
               h_sapiens_ensembl_id = "Homo_sapiens",
               m_musculus_ensembl_id = "Mus_musculus", 
               r_norvegicus_ensembl_id = "Rattus_norvegicus",
               d_rerio_ensembl_id = "Danio_rerio", 
               d_melanogaster_ensembl_id = "Drosophila_melanogaster", 
               c_elegans_ensembl_id = "Caenorhabditis_elegans",
            stop("Error: Invalid map_species in CoSIAn Object. Make sure the 
            species in the map_species slot are an avalible model organism and 
                 are in the correct format."))
    })
    map_species <- Species_SWITCH(map_species)
    metric_type <- object@metric_type
    metric_type <- as.character(metric_type)
    
    #Load CoSIAdata and merge into merged_CoSIAdata that is species specific to the mapped_species slot
    CoSIAdata_load <- function(map_species) {
      eh <- ExperimentHub::ExperimentHub()
      merged_CoSIAdata <- data.frame(matrix(ncol = 7, nrow = 0))
      colnames(merged_CoSIAdata) <- c("Anatomical_entity_name", "Ensembl_ID",
                              "Sample_size", "VST", "Experiment_ID",
                              "Anatomical_entity_ID","Species")
      if (any(map_species == "Mus_musculus")) {
        mm_EH_File<-eh[["EH7859"]]
        merged_CoSIAdata <- rbind(merged_CoSIAdata, mm_EH_File)
        merged_CoSIAdata <- as.data.frame(merged_CoSIAdata)
        }
      else if (any(map_species == "Rattus_norvegicus")) {
        rn_EH_File<-eh[["EH7860"]]
        merged_CoSIAdata <- rbind(merged_CoSIAdata, rn_EH_File)
        merged_CoSIAdata <- as.data.frame(merged_CoSIAdata)
      }
      else if (any(map_species == "Danio_rerio")) {
        dr_EH_File<-eh[["EH7861"]]
        merged_CoSIAdata <- rbind(merged_CoSIAdata, dr_EH_File)
        merged_CoSIAdata <- as.data.frame(merged_CoSIAdata)
      }
      else if (any(map_species == "Homo_sapiens")) {
        hs_EH_File<-eh[["EH7858"]]
        merged_CoSIAdata <- rbind(merged_CoSIAdata, hs_EH_File)
        merged_CoSIAdata <- as.data.frame(merged_CoSIAdata)
      }
      else if (any(map_species == "Caenorhabditis_elegans")) {
        ce_EH_File<-eh[["EH7863"]]
        merged_CoSIAdata <- rbind(merged_CoSIAdata, ce_EH_File)
        merged_CoSIAdata <- as.data.frame(merged_CoSIAdata)
      }
      else if (any(map_species == "Drosophila_melanogaster")) {
        dm_EH_File<-eh[["EH7862"]]
        merged_CoSIAdata <- rbind(merged_CoSIAdata, dm_EH_File)
        merged_CoSIAdata <- as.data.frame(merged_CoSIAdata)
      }
      else {
        stop("Error: map_species in CoSIAn Object. Make sure the species in the
        map_species slot are avalible organisms through CoSIA and are in the
        correct format.")
      }
      return(merged_CoSIAdata)
    }

    merged_CoSIAdata <- lapply(map_species, CoSIAdata_load)
    filter_species <- as.data.frame(do.call(rbind, merged_CoSIAdata))

    CV_function <- function(x, na.rm = FALSE) {
        stopifnot(is.numeric(x))
        if (na.rm)
            x <- x[!is.na(x)]
        if (any(x < 0))
            stop("Your data must be greater than zero!")
        stats::sd(x, na.rm = FALSE)/stats::median(x)
    }

    # Code adapted from BioQC
    DS_function <- function(type, mat) {
        Shannon_Entropy <- function(x) ifelse(x == 0, 0, x * log2(x))
        Diversity <- function(x) {
            rel <- x/sum(x, na.rm = TRUE)
            -sum(Shannon_Entropy(rel))
        }
        if (type == "Specificity") {
            mat.rel <- apply(mat, 2L, function(x) x/sum(x, na.rm = TRUE))
            speci <- apply(mat.rel, 1L, function(x) {
                gm <- mean(x, na.rm = TRUE)
                rel <- x/gm
                mean(Shannon_Entropy(rel), na.rm = TRUE)
            })
            speci <- speci/log2(ncol(mat))
            return(speci)
        }
        if (type == "Diversity") {
            res <- apply(mat, 2L, Diversity)
            res <- res/log2(nrow(mat))
            return(res)
        } else {
            stop("Incorrect Type")
        }
    }



    CV_Tissue <- function(map_species, map_tissues) {
        filter_tissue <- dplyr::filter(filter_species, Anatomical_entity_name %in% map_tissues)
        id <- as.vector(t(id_dataframe))
        filter_gene <- dplyr::filter(filter_tissue, Ensembl_ID %in% id)
        filter_gex <- tidyr::separate_rows(filter_gene, VST)
        filter_gex$VST <- as.numeric(filter_gex$VST)
        CV_Tissue <- filter_gex %>%
            group_by(Ensembl_ID, Anatomical_entity_name, Species) %>%
            summarise(CV_Tissue = CV_function(VST, na.rm = FALSE))
        cv_tissue <- tidyr::pivot_wider(CV_Tissue, names_from = Species, values_from = CV_Tissue)
        colnames(cv_tissue)[which(names(cv_tissue) == "Homo_sapiens")] <- "h_sapiens_CV_tissue"
        colnames(cv_tissue)[which(names(cv_tissue) == "Mus_musculus")] <- "m_musculus_CV_tissue"
        colnames(cv_tissue)[which(names(cv_tissue) == "Rattus_norvegicus")] <- "r_norvegicus_CV_tissue"
        colnames(cv_tissue)[which(names(cv_tissue) == "Danio_rerio")] <- "d_rerio_CV_Tissue"
        colnames(cv_tissue)[which(names(cv_tissue) == "Drosophila_melanogaster")] <- "d_melanogaster_CV_tissue"
        colnames(cv_tissue)[which(names(cv_tissue) == "Caenorhabditis_elegans")] <- "c_elegans_CV_tissue"
        for (i in seq_len(length(colnames(id_dataframe)))) {
            id <- colnames(id_dataframe)[i]
            id_dataframe <- id_dataframe %>%
                merge(., cv_tissue, by.x = id, by.y = "Ensembl_ID")
        }
        id_dataframe <- id_dataframe[, colSums(is.na(id_dataframe)) < nrow(id_dataframe)]
        id_dataframe <- id_dataframe %>%
            dplyr::select(order(colnames(id_dataframe), decreasing = TRUE))
        id_dataframe <- id_dataframe %>%
            dplyr::filter(dplyr::if_all(tidyselect::starts_with("Anatomical_entity_name"), ~Anatomical_entity_name.x == .x))
        duplicated_columns <- duplicated(as.list(id_dataframe))
        id_dataframe <- id_dataframe[!duplicated_columns]
        id_dataframe <- id_dataframe %>%
            dplyr::rename_with(~stringr::str_remove(., c(".x")))
        CV_Tissue <- id_dataframe %>%
            dplyr::rename_with(~stringr::str_remove(., c(".y")))
        colnames(CV_Tissue)[which(names(CV_Tissue) == "Anatomical_enti_name.y")] <- "Anatomical_entity_name"
        CV_Tissue <- unique(CV_Tissue)
        return(CV_Tissue)
    }

    CV_Species <- function(map_species, map_tissues) {
        id <- as.vector(t(id_dataframe))
        filter_gene <- dplyr::filter(filter_species, Ensembl_ID %in% id)
        filter_gex <- tidyr::separate_rows(filter_gene, VST)
        filter_gex$VST <- as.numeric(filter_gex$VST)
        CV_Species <- filter_gex %>%
            group_by(Ensembl_ID, Species) %>%
            summarise(CV_Species = CV_function(VST, na.rm = FALSE))
        cv_species <- tidyr::pivot_wider(CV_Species, names_from = Species, values_from = CV_Species)
        colnames(cv_species)[which(names(cv_species) == "Homo_sapiens")] <- "h_sapiens_CV_species"
        colnames(cv_species)[which(names(cv_species) == "Mus_musculus")] <- "m_musculus_CV_species"
        colnames(cv_species)[which(names(cv_species) == "Rattus_norvegicus")] <- "r_norvegicus_CV_species"
        colnames(cv_species)[which(names(cv_species) == "Danio_rerio")] <- "d_rerio_CV_species"
        colnames(cv_species)[which(names(cv_species) == "Drosophila_melanogaster")] <- "d_melanogaster_CV_species"
        colnames(cv_species)[which(names(cv_species) == "Caenorhabditis_elegans")] <- "c_elegans_CV_species"
        for (i in seq_len(length(colnames(id_dataframe)))) {
            id <- colnames(id_dataframe)[i]
            id_dataframe <- id_dataframe %>%
                merge(., cv_species, by.x = id, by.y = "Ensembl_ID")
        }
        id_dataframe <- id_dataframe[, colSums(is.na(id_dataframe)) < nrow(id_dataframe)]
        id_dataframe <- id_dataframe %>%
            dplyr::select(order(colnames(id_dataframe), decreasing = TRUE))
        id_dataframe <- id_dataframe %>%
            dplyr::rename_with(~stringr::str_remove(., c(".x")))
        CV_Species <- id_dataframe %>%
            dplyr::rename_with(~stringr::str_remove(., c(".y")))
        CV_Species <- unique(CV_Species)
        return(CV_Species)
    }

    # DS_Gene : output genes restricted by mapped tissues and gene set
    DS_Gene <- function(map_species, map_tissues) {
        filter_tissue <- dplyr::filter(filter_species, Anatomical_entity_name %in% map_tissues)
        id <- as.vector(t(id_dataframe))
        filter_gene <- dplyr::filter(filter_tissue, Ensembl_ID %in% id)
        filter_gene$Scaled_Median_VST <- as.numeric(filter_gene$Scaled_Median_VST)
        filter_gex <- dplyr::select(filter_gene, Anatomical_entity_name, Scaled_Median_VST, Ensembl_ID)
        filter_gex_D <- filter_gex %>%
            tidyr::pivot_wider(names_from = Ensembl_ID, values_from = Scaled_Median_VST)
        filter_gex_D <- filter_gex_D %>%
            remove_rownames %>%
            tibble::column_to_rownames(var = "Anatomical_entity_name")
        filter_gex_D <- data.matrix(filter_gex_D, )
        ENTROPY_DIVERSITY_G <- data.frame(DS_function("Diversity", filter_gex_D))  # across genes
        colnames(ENTROPY_DIVERSITY_G)[which(names(ENTROPY_DIVERSITY_G) == "DS_function..Diversity...filter_gex_D.")] <- "Diversity"

        filter_gex <- data.frame(filter_gex)
        filter_gex_S <- filter_gex %>%
            tidyr::pivot_wider(names_from = Anatomical_entity_name, values_from = Scaled_Median_VST)
        filter_gex_S <- filter_gex_S %>%
            remove_rownames %>%
            column_to_rownames(var = "Ensembl_ID")
        filter_gex_S <- data.matrix(filter_gex_S, )
        ENTROPY_SPECIFITY_G <- data.frame(DS_function("Specificity", filter_gex_S))  # across tissues
        colnames(ENTROPY_SPECIFITY_G)[which(names(ENTROPY_SPECIFITY_G) == "DS_function..Specificity...filter_gex_S.")] <- "Specificity"

        DS <- merge(ENTROPY_SPECIFITY_G, ENTROPY_DIVERSITY_G, by = "row.names")
        colnames(DS)[which(names(DS) == "Row.names")] <- "Ensembl_ID"
        DS$Ensembl_ID <- as.character(DS$Ensembl_ID)
        Species <- dplyr::select(filter_gene, Species, Ensembl_ID)
        DS <- merge(DS, Species, by = "Ensembl_ID")
        DS <- data.frame(unique(DS))
        rownames(DS) <- NULL
        return(DS)
    }

    # DS_Tissues: output is tissues restricted to mapped tissues and gene set
    DS_Tissue <- function(map_species, map_tissues) {
        DS <- data.frame(matrix(ncol = 4, nrow = 0))
        colnames(DS)[which(names(DS) == "map_tissues")] <- "Anatomical_entity_name"
        for (x in seq_len(length(map_species))) {
            filter_species <- dplyr::filter(filter_species, Species == map_species[x])
            filter_tissue <- dplyr::filter(filter_species, Anatomical_entity_name %in% map_tissues)
            id <- as.vector(t(id_dataframe))
            filter_gene <- dplyr::filter(filter_tissue, Ensembl_ID %in% id)
            filter_gene$Scaled_Median_VST <- as.numeric(filter_gene$Scaled_Median_VST)
            filter_gex <- dplyr::select(filter_gene, Anatomical_entity_name, Scaled_Median_VST, Ensembl_ID)

            filter_gex_D <- filter_gex %>%
                tidyr::pivot_wider(names_from = Anatomical_entity_name, values_from = Scaled_Median_VST)
            filter_gex_D <- filter_gex_D %>%
                remove_rownames %>%
                tibble::column_to_rownames(var = "Ensembl_ID")
            filter_gex_D <- data.matrix(filter_gex_D, )
            ENTROPY_DIVERSITY_T <- data.frame(DS_function("Diversity", filter_gex_D))  # across genes
            colnames(ENTROPY_DIVERSITY_T)[which(names(ENTROPY_DIVERSITY_T) == "DS_function..Diversity...filter_gex_D.")] <- "Diversity"

            filter_gex <- data.frame(filter_gex)
            filter_gex_S <- filter_gex %>%
                tidyr::pivot_wider(names_from = Ensembl_ID, values_from = Scaled_Median_VST)
            filter_gex_S <- filter_gex_S %>%
                remove_rownames %>%
                column_to_rownames(var = "Anatomical_entity_name")
            filter_gex_S <- data.matrix(filter_gex_S, )
            ENTROPY_SPECIFITY_T <- data.frame(DS_function("Specificity", filter_gex_S))  # across tissues
            colnames(ENTROPY_SPECIFITY_T)[which(names(ENTROPY_SPECIFITY_T) == "DS_function..Specificity...filter_gex_S.")] <- "Specificity"

            SDS <- merge(ENTROPY_SPECIFITY_T, ENTROPY_DIVERSITY_T, by = "row.names")
            colnames(SDS)[which(names(SDS) == "Row.names")] <- "Anatomical_entity_name"
            SDS$Anatomical_entity_name <- as.character(SDS$Anatomical_entity_name)
            Species <- dplyr::select(filter_gene, Species, Anatomical_entity_name)
            SDS <- merge(SDS, Species, by = "Anatomical_entity_name")
            DS <- rbind(DS, SDS)
        }
        DS <- data.frame(unique(DS))
        rownames(DS) <- NULL
        return(DS)

    }
    # DS_Gene_all: outputs genes only restricted to selected genes across all tissues
    DS_Gene_all <- function(map_species) {
        DS <- data.frame(matrix(ncol = 4, nrow = 0))
        colnames(DS)[1] <- "Ensembl_ID"
        for (x in seq_len(length(map_species))) {
            filter_species <- dplyr::filter(filter_species, Species == map_species[x])
            id <- as.vector(t(id_dataframe))
            filter_gene <- dplyr::filter(filter_species, Ensembl_ID %in% id)
            filter_gene$Scaled_Median_VST <- as.numeric(filter_gene$Scaled_Median_VST)
            filter_gex <- dplyr::select(filter_gene, Anatomical_entity_name, Scaled_Median_VST, Ensembl_ID)

            filter_gex_D <- filter_gex %>%
                tidyr::pivot_wider(names_from = Ensembl_ID, values_from = Scaled_Median_VST)
            filter_gex_D <- filter_gex_D %>%
                remove_rownames %>%
                tibble::column_to_rownames(var = "Anatomical_entity_name")
            filter_gex_D <- data.matrix(filter_gex_D, )
            ENTROPY_DIVERSITY_G <- data.frame(DS_function("Diversity", filter_gex_D))  # across genes
            colnames(ENTROPY_DIVERSITY_G)[which(names(ENTROPY_DIVERSITY_G) == "DS_function..Diversity...filter_gex_D.")] <- "Diversity"

            filter_gex <- data.frame(filter_gex)
            filter_gex_S <- filter_gex %>%
                tidyr::pivot_wider(names_from = Anatomical_entity_name, values_from = Scaled_Median_VST)
            filter_gex_S <- filter_gex_S %>%
                remove_rownames %>%
                column_to_rownames(var = "Ensembl_ID")
            filter_gex_S <- data.matrix(filter_gex_S, )
            ENTROPY_SPECIFITY_G <- data.frame(DS_function("Specificity", filter_gex_S))  # across tissues
            colnames(ENTROPY_SPECIFITY_G)[which(names(ENTROPY_SPECIFITY_G) == "DS_function..Specificity...filter_gex_S.")] <- "Specificity"

            SDS <- merge(ENTROPY_SPECIFITY_G, ENTROPY_DIVERSITY_G, by = "row.names")
            colnames(SDS)[which(names(SDS) == "Row.names")] <- "Ensembl_ID"
            SDS$Ensembl_ID <- as.character(SDS$Ensembl_ID)
            Species <- dplyr::select(filter_gene, Species, Ensembl_ID)
            SDS <- merge(SDS, Species, by = "Ensembl_ID")
            DS <- rbind(DS, SDS)
        }
        DS <- data.frame(unique(DS))
        rownames(DS) <- NULL
        return(DS)
    }

    # DS_Tissues_All: output is tissues restricted to mapped tissues across all genes
    DS_Tissue_all <- function(map_species, map_tissues) {
        DS <- data.frame(matrix(ncol = 4, nrow = 0))
        colnames(DS)[which(names(DS) == "map_tissues")] <- "Anatomical_entity_name"
        for (x in seq_len(length(map_species))) {
            filter_species <- dplyr::filter(filter_species, Species == map_species[x])
            filter_tissue <- dplyr::filter(filter_species, Anatomical_entity_name %in% map_tissues)
            filter_tissue$Scaled_Median_VST <- as.numeric(filter_tissue$Scaled_Median_VST)
            filter_gex <- dplyr::select(filter_tissue, Anatomical_entity_name, Scaled_Median_VST, Ensembl_ID)

            filter_gex_D <- filter_gex %>%
                tidyr::pivot_wider(names_from = Anatomical_entity_name, values_from = Scaled_Median_VST)
            filter_gex_D <- filter_gex_D %>%
                remove_rownames %>%
                tibble::column_to_rownames(var = "Ensembl_ID")
            filter_gex_D <- data.matrix(filter_gex_D, )
            ENTROPY_DIVERSITY_T <- data.frame(DS_function("Diversity", filter_gex_D))  # across genes
            colnames(ENTROPY_DIVERSITY_T)[which(names(ENTROPY_DIVERSITY_T) == "DS_function..Diversity...filter_gex_D.")] <- "Diversity"

            filter_gex <- data.frame(filter_gex)
            filter_gex_S <- filter_gex %>%
                tidyr::pivot_wider(names_from = Ensembl_ID, values_from = Scaled_Median_VST)
            filter_gex_S <- filter_gex_S %>%
                remove_rownames %>%
                column_to_rownames(var = "Anatomical_entity_name")
            filter_gex_S <- data.matrix(filter_gex_S, )
            ENTROPY_SPECIFITY_T <- data.frame(DS_function("Specificity", filter_gex_S))  # across tissues
            colnames(ENTROPY_SPECIFITY_T)[which(names(ENTROPY_SPECIFITY_T) == "DS_function..Specificity...filter_gex_S.")] <- "Specificity"

            SDS <- merge(ENTROPY_SPECIFITY_T, ENTROPY_DIVERSITY_T, by = "row.names")
            colnames(SDS)[which(names(SDS) == "Row.names")] <- "Anatomical_entity_name"
            SDS$Anatomical_entity_name <- as.character(SDS$Anatomical_entity_name)
            Species <- dplyr::select(filter_tissue, Species, Anatomical_entity_name)
            SDS <- merge(SDS, Species, by = "Anatomical_entity_name")
            DS <- rbind(DS, SDS)
        }
        DS <- data.frame(unique(DS))
        rownames(DS) <- NULL
        return(DS)
    }


    if (metric_type == "CV_Tissue") {
        CV_Tissue <- CV_Tissue(map_species, map_tissues)
        object@metric <- CV_Tissue
    } else if (metric_type == "CV_Species") {
        CV_Species <- CV_Species(map_species, map_tissues)
        object@metric <- CV_Species
    } else if (metric_type == "DS_Gene") {
        DS_Gene <- DS_Gene(map_species, map_tissues)
        object@metric <- DS_Gene
    } else if (metric_type == "DS_Tissue") {
        DS_Tissue <- DS_Tissue(map_species, map_tissues)
        object@metric <- DS_Tissue
    } else if (metric_type == "DS_Gene_all") {
        DS_Gene_all <- DS_Gene_all(map_species, map_tissues)
        object@metric <- DS_Gene_all
    } else if (metric_type == "DS_Tissue_all") {
        DS_Tissue_all <- DS_Tissue_all(map_species, map_tissues)
        object@metric <- DS_Tissue_all
    } else (stop("Error: invalid metric type"))

    return(object)
})
