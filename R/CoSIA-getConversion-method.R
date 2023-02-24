#' getConversion Generic
#'
#' @param object CoSIAn object with all user accessible slots filled
#'
#' @return initializes a generic function for getConversion as preparation for defining the getConversion Method
#' @export
#' @examples
#' Kidney_Genes <- CoSIAn(
#'     gene_set = c(
#'         "ENSG00000008710", "ENSG00000118762",
#'         "ENSG00000152217"
#'     ), i_species = "h_sapiens", input_id = "Ensembl_id",
#'     o_species = c(
#'         "d_melanogaster", "m_musculus", "h_sapiens", "d_rerio", "c_elegans",
#'         "r_norvegicus"
#'     ), output_ids = c("Ensembl_id", "Symbol"), mapping_tool =
#'         "annotationDBI", ortholog_database = "HomoloGene", map_tissues = "heart",
#'     map_species = c("m_musculus"), metric_type = "DS_Gene"
#' )
#' Kidney_gene_conversion <- CoSIA::getConversion(Kidney_Genes)
setGeneric("getConversion", function(object) standardGeneric("getConversion"))

#' getConversion Method
#'
#' @param object CoSIAn object with all user accessible slots filled
#'
#' @return CoSIAn object with converted_id slot filled
#' @export
#' @importFrom magrittr %>%
#' @importFrom stats na.omit
#' @importFrom tidyselect contains
#'
#' @examples
#' Kidney_Genes <- CoSIAn(
#'     gene_set = c(
#'         "ENSG00000008710", "ENSG00000118762",
#'         "ENSG00000152217"
#'     ), i_species = "h_sapiens", input_id = "Ensembl_id",
#'     o_species = c(
#'         "d_melanogaster", "m_musculus", "h_sapiens", "d_rerio", "c_elegans",
#'         "r_norvegicus"
#'     ), output_ids = c("Ensembl_id", "Symbol"), mapping_tool =
#'         "annotationDBI", ortholog_database = "HomoloGene", map_tissues = "heart",
#'     map_species = c("m_musculus"), metric_type = "DS_Gene"
#' )
#' Kidney_gene_conversion <- CoSIA::getConversion(Kidney_Genes)
setMethod("getConversion", signature(object = "CoSIAn"), function(object) {
    # user's input of the function
    # Set each part of the object that this method uses into their own variable that will be used inside the code
    input_species <- object@i_species
    input_id <- object@input_id
    input <- object@gene_set
    input <- unique(input)
    if (input_id == "Ensembl_id") {
        input <- remove_version_numbers(input, input_species)
        # puts the genes through the remove version number function to remove the version number (everything after the dot)
    }
    input <- as.character(input)
    output_ids <- object@output_ids
    output_species <- object@o_species
    tool <- object@mapping_tool
    ortholog_database <- object@ortholog_database
    species_data <- data.frame(input)
    c_name <- tolower(paste(input_species, input_id, sep = "_"))
    colnames(species_data)[which(names(species_data) == "input")] <- c_name
    # changes name to more formal names
    for (index in seq(length(output_species))) {
        Filter_I_Species <- switch(input_species,
            h_sapiens = {
                hs_data <- h_sapiens(input_id, input, output_ids, output_species[index], tool, ortholog_database)
                hs_data <- data.frame(lapply(hs_data, as.character))
                result <- try(
                    {
                        dplyr::full_join(species_data, hs_data)
                    },
                    silent = TRUE
                )
                if (is(result, "try-error")) {
                    warning("No orthologs were found for ", paste(output_species[index], sep = " "), " across all the genes provided by the user.")
                } else {
                    species_data <- dplyr::full_join(species_data, hs_data)
                }
            },
            m_musculus = {
                mm_data <- m_musculus(input_id, input, output_ids, output_species[index], tool, ortholog_database)
                mm_data <- data.frame(lapply(mm_data, as.character))
                result <- try(
                    {
                        dplyr::full_join(species_data, mm_data)
                    },
                    silent = TRUE
                )
                if (is(result, "try-error")) {
                    warning("No orthologs were found for ", paste(output_species[index], sep = " "), " across all the genes provided by the user.")
                } else {
                    species_data <- dplyr::full_join(species_data, mm_data)
                }
            },
            r_norvegicus = {
                rn_data <- r_norvegicus(input_id, input, output_ids, output_species[index], tool, ortholog_database)
                rn_data <- data.frame(lapply(rn_data, as.character))
                result <- try(
                    {
                        dplyr::full_join(species_data, rn_data)
                    },
                    silent = TRUE
                )
                if (is(result, "try-error")) {
                    warning("No orthologs were found for ", paste(output_species[index], sep = " "), " across all the genes provided by the user.")
                } else {
                    species_data <- dplyr::full_join(species_data, rn_data)
                }
            },
            d_rerio = {
                dr_data <- d_rerio(input_id, input, output_ids, output_species[index], tool, ortholog_database)
                dr_data <- data.frame(lapply(dr_data, as.character))
                result <- try(
                    {
                        dplyr::full_join(species_data, dr_data)
                    },
                    silent = TRUE
                )
                if (is(result, "try-error")) {
                    warning("No orthologs were found for ", paste(output_species[index], sep = " "), " across all the genes provided by the user.")
                } else {
                    species_data <- dplyr::full_join(species_data, dr_data)
                }
            },
            c_elegans = {
                ce_data <- c_elegans(input_id, input, output_ids, output_species[index], tool, ortholog_database)
                ce_data <- data.frame(lapply(ce_data, as.character))
                result <- try(
                    {
                        dplyr::full_join(species_data, ce_data)
                    },
                    silent = TRUE
                )
                if (is(result, "try-error")) {
                    warning("No orthologs were found for ", paste(output_species[index], sep = " "), " across all the genes provided by the user.")
                } else {
                    species_data <- dplyr::full_join(species_data, ce_data)
                }
            },
            d_melanogaster = {
                dm_data <- d_melanogaster(input_id, input, output_ids, output_species[index], tool, ortholog_database)
                dm_data <- data.frame(lapply(dm_data, as.character))
                result <- try(
                    {
                        dplyr::full_join(species_data, dm_data)
                    },
                    silent = TRUE
                )
                if (is(result, "try-error")) {
                    warning("No orthologs were found for ", paste(output_species[index], sep = " "), " across all the genes provided by the user.")
                } else {
                    species_data <- dplyr::full_join(species_data, dm_data)
                }
            },
            stop("Error: Invalid i_species in CoSIAn Object. Make sure the species in the i_species slot is an avalible model
                                    organism and is in the correct format.")
        )
    }
    object@converted_id <- data.frame(species_data)
    # converted_id<-slot(object, converted_id)
    return(object)
})

### getConversion Functions

## Species Functions

# Caenorhabditis elegans

c_elegans <- function(input_id, input_dataset, output_ids, output_species, tool, ortholog_database) { # input data funneling in from the cross species conversion function
    if (tool == "biomaRt") { # if the user has choose the tool biomart this it the path the codes follows
        Filter_BO_CE <- switch(output_species,
            h_sapiens = {
                output_data <- biomaRt(input_id, input_dataset, output_ids, input_species = "c_elegans", output_species, 9606, "celegans_gene_ensembl", "hsapiens_gene_ensembl", ortholog_database)
            },
            m_musculus = {
                output_data <- biomaRt(input_id, input_dataset, output_ids, input_species = "c_elegans", output_species, 10090, "celegans_gene_ensembl", "mmusculus_gene_ensembl", ortholog_database)
            },
            r_norvegicus = {
                output_data <- biomaRt(input_id, input_dataset, output_ids, input_species = "c_elegans", output_species, 10116, "celegans_gene_ensembl", "rnorvegicus_gene_ensembl", ortholog_database)
            },
            d_rerio = {
                output_data <- biomaRt(input_id, input_dataset, output_ids, input_species = "c_elegans", output_species, 7955, "celegans_gene_ensembl", "drerio_gene_ensembl", ortholog_database)
            },
            c_elegans = {
                output_data <- biomaRt(input_id, input_dataset, output_ids, input_species = "c_elegans", output_species, 6239, "celegans_gene_ensembl", "celegans_gene_ensembl", ortholog_database)
            },
            d_melanogaster = {
                output_data <- biomaRt(input_id, input_dataset, output_ids, input_species = "c_elegans", output_species, 7227, "celegans_gene_ensembl", "dmelanogaster_gene_ensembl", ortholog_database)
            },
            stop("Error: Invalid o_species in CoSIAn Object. Make sure the species in the o_species slot is an avalible model
                                    organism and is in the correct format.")
        )
    } else if (tool == "annotationDBI") { # code follows this path if the user chooses annotationDBI as their tool of choice
        Filter_AO_CE <- switch(output_species,
            h_sapiens = {
                output_data <- annotationDBI(input_id, input_dataset, output_ids, input_species = "c_elegans", output_species, 9606, org.Ce.eg.db::org.Ce.eg.db, org.Hs.eg.db::org.Hs.eg.db, ortholog_database)
            },
            m_musculus = {
                output_data <- annotationDBI(input_id, input_dataset, output_ids, input_species = "c_elegans", output_species, 10090, org.Ce.eg.db::org.Ce.eg.db, org.Mm.eg.db::org.Mm.eg.db, ortholog_database)
            },
            r_norvegicus = {
                output_data <- annotationDBI(input_id, input_dataset, output_ids, input_species = "c_elegans", output_species, 10116, org.Ce.eg.db::org.Ce.eg.db, org.Rn.eg.db::org.Rn.eg.db, ortholog_database)
            },
            d_rerio = {
                output_data <- annotationDBI(input_id, input_dataset, output_ids, input_species = "c_elegans", output_species, 7955, org.Ce.eg.db::org.Ce.eg.db, org.Dr.eg.db::org.Dr.eg.db, ortholog_database)
            },
            c_elegans = {
                output_data <- annotationDBI(input_id, input_dataset, output_ids, input_species = "c_elegans", output_species, 6239, org.Ce.eg.db::org.Ce.eg.db, org.Ce.eg.db::org.Ce.eg.db, ortholog_database)
            },
            d_melanogaster = {
                output_data <- annotationDBI(input_id, input_dataset, output_ids, input_species = "c_elegans", output_species, 7227, org.Ce.eg.db::org.Ce.eg.db, org.Dm.eg.db::org.Dm.eg.db, ortholog_database)
            },
            stop("Error: Invalid o_species in CoSIAn Object. Make sure the species in the o_species slot is an avalible model
                                    organism and is in the correct format.")
        )
    } else { # code follows this path if the tool that was inputted into the system does not match.
        stop("Error: Invalid tool in CoSIAn Object. Make sure the tool is either annotationDBI or biomaRt.")
    }
}

# Drosophila Melanogaster - Fruit Fly

d_melanogaster <- function(input_id, input_dataset, output_ids, output_species, tool, ortholog_database) { # input data funneling in from the cross species conversion function
    if (tool == "biomaRt") { # if the user has choose the tool biomart this it the path the codes follows
        Filter_BO_DM <- switch(output_species,
            h_sapiens = {
                output_data <- biomaRt(input_id, input_dataset, output_ids, input_species = "d_melanogaster", output_species, 9606, "dmelanogaster_gene_ensembl", "hsapiens_gene_ensembl", ortholog_database)
            },
            m_musculus = {
                output_data <- biomaRt(input_id, input_dataset, output_ids, input_species = "d_melanogaster", output_species, 10090, "dmelanogaster_gene_ensembl", "mmusculus_gene_ensembl", ortholog_database)
            },
            r_norvegicus = {
                output_data <- biomaRt(input_id, input_dataset, output_ids, input_species = "d_melanogaster", output_species, 10116, "dmelanogaster_gene_ensembl", "rnorvegicus_gene_ensembl", ortholog_database)
            },
            d_rerio = {
                output_data <- biomaRt(input_id, input_dataset, output_ids, input_species = "d_melanogaster", output_species, 7955, "dmelanogaster_gene_ensembl", "drerio_gene_ensembl", ortholog_database)
            },
            c_elegans = {
                output_data <- biomaRt(input_id, input_dataset, output_ids, input_species = "d_melanogaster", output_species, 6239, "dmelanogaster_gene_ensembl", "celegans_gene_ensembl", ortholog_database)
            },
            d_melanogaster = {
                output_data <- biomaRt(input_id, input_dataset, output_ids, input_species = "d_melanogaster", output_species, 7227, "dmelanogaster_gene_ensembl", "dmelanogaster_gene_ensembl", ortholog_database)
            },
            stop("Error: Invalid o_species in CoSIAn Object. Make sure the species in the o_species slot is an avalible model
                                    organism and is in the correct format.")
        )
    } else if (tool == "annotationDBI") { # code follows this path if the user chooses annotationDBI as their tool of choice
        Filter_AO_DM <- switch(output_species,
            h_sapiens = {
                output_data <- annotationDBI(input_id, input_dataset, output_ids, input_species = "d_melanogaster", output_species, 9606, org.Dm.eg.db::org.Dm.eg.db, org.Hs.eg.db::org.Hs.eg.db, ortholog_database)
            },
            m_musculus = {
                output_data <- annotationDBI(input_id, input_dataset, output_ids, input_species = "d_melanogaster", output_species, 10090, org.Dm.eg.db::org.Dm.eg.db, org.Mm.eg.db::org.Mm.eg.db, ortholog_database)
            },
            r_norvegicus = {
                output_data <- annotationDBI(input_id, input_dataset, output_ids, input_species = "d_melanogaster", output_species, 10116, org.Dm.eg.db::org.Dm.eg.db, org.Rn.eg.db::org.Rn.eg.db, ortholog_database)
            },
            d_rerio = {
                output_data <- annotationDBI(input_id, input_dataset, output_ids, input_species = "d_melanogaster", output_species, 7955, org.Dm.eg.db::org.Dm.eg.db, org.Dr.eg.db::org.Dr.eg.db, ortholog_database)
            },
            c_elegans = {
                output_data <- annotationDBI(input_id, input_dataset, output_ids, input_species = "d_melanogaster", output_species, 6239, org.Dm.eg.db::org.Dm.eg.db, org.Ce.eg.db::org.Ce.eg.db, ortholog_database)
            },
            d_melanogaster = {
                output_data <- annotationDBI(input_id, input_dataset, output_ids, input_species = "d_melanogaster", output_species, 7227, org.Dm.eg.db::org.Dm.eg.db, org.Dm.eg.db::org.Dm.eg.db, ortholog_database)
            },
            stop("Error: Invalid o_species in CoSIAn Object. Make sure the species in the o_species slot is an avalible model
                                    organism and is in the correct format.")
        )
    } else { # code follows this path if the tool that was inputted into the system does not match.
        stop("Error: Invalid tool in CoSIAn Object. Make sure the tool is either annotationDBI or biomaRt.")
    }
}

# d_rerio - Zebrafish

d_rerio <- function(input_id, input_dataset, output_ids, output_species, tool, ortholog_database) { # input data funneling in from the cross species conversion function
    if (tool == "biomaRt") { # if the user has choose the tool biomart this it the path the codes follows
        Filter_BO_DR <- switch(output_species,
            h_sapiens = {
                output_data <- biomaRt(input_id, input_dataset, output_ids, input_species = "d_rerio", output_species, 9606, "drerio_gene_ensembl", "hsapiens_gene_ensembl", ortholog_database)
            },
            c_elegans = {
                output_data <- biomaRt(input_id, input_dataset, output_ids, input_species = "d_rerio", output_species, 6239, "drerio_gene_ensembl", "celegans_gene_ensembl", ortholog_database)
            },
            d_melanogaster = {
                output_data <- biomaRt(input_id, input_dataset, output_ids, input_species = "d_rerio", output_species, 7227, "drerio_gene_ensembl", "dmelanogaster_gene_ensembl", ortholog_database)
            },
            m_musculus = {
                output_data <- biomaRt(input_id, input_dataset, output_ids, input_species = "d_rerio", output_species, 10090, "drerio_gene_ensembl", "mmusculus_gene_ensembl", ortholog_database)
            },
            d_rerio = {
                output_data <- biomaRt(input_id, input_dataset, output_ids, input_species = "d_rerio", output_species, 7955, "drerio_gene_ensembl", "drerio_gene_ensembl", ortholog_database)
            },
            r_norvegicus = {
                output_data <- biomaRt(input_id, input_dataset, output_ids, input_species = "d_rerio", output_species, 10116, "drerio_gene_ensembl", "rnorvegicus_gene_ensembl", ortholog_database)
            },
            stop("Error: Invalid o_species in CoSIAn Object. Make sure the species in the o_species slot is an avalible model
                                    organism and is in the correct format.")
        )
    } else if (tool == "annotationDBI") { # code follows this path if the user chooses annotationDBI as their tool of choice
        Filter_AO_DR <- switch(output_species,
            h_sapiens = {
                output_data <- annotationDBI(input_id, input_dataset, output_ids,
                    input_species = "d_rerio", output_species, 9606, org.Dr.eg.db::org.Dr.eg.db,
                    org.Hs.eg.db::org.Hs.eg.db, ortholog_database
                )
            },
            m_musculus = {
                output_data <- annotationDBI(input_id, input_dataset, output_ids,
                    input_species = "d_rerio", output_species, 10090, org.Dr.eg.db::org.Dr.eg.db,
                    org.Mm.eg.db::org.Mm.eg.db, ortholog_database
                )
            },
            d_melanogaster = {
                output_data <- annotationDBI(input_id, input_dataset, output_ids,
                    input_species = "d_rerio", output_species, 7227, org.Dr.eg.db::org.Dr.eg.db,
                    org.Dm.eg.db::org.Dm.eg.db, ortholog_database
                )
            },
            d_rerio = {
                output_data <- annotationDBI(input_id, input_dataset, output_ids,
                    input_species = "d_rerio", output_species, 7955, org.Dr.eg.db::org.Dr.eg.db,
                    org.Dr.eg.db::org.Dr.eg.db, ortholog_database
                )
            },
            c_elegans = {
                output_data <- annotationDBI(input_id, input_dataset, output_ids,
                    input_species = "d_rerio", output_species, 6239, org.Dr.eg.db::org.Dr.eg.db,
                    org.Ce.eg.db::org.Ce.eg.db, ortholog_database
                )
            },
            r_norvegicus = {
                output_data <- annotationDBI(input_id, input_dataset, output_ids,
                    input_species = "d_rerio", output_species, 10116, org.Dr.eg.db::org.Dr.eg.db,
                    org.Rn.eg.db::org.Rn.eg.db, ortholog_database
                )
            },
            stop("Error: Invalid o_species in CoSIAn Object. Make sure the species in the o_species slot is an avalible model
                                    organism and is in the correct format.")
        )
    } else { # code follows this path if the tool that was inputted into the system does not match.
        stop("Error: Invalid tool in CoSIAn Object. Make sure the tool is either annotationDBI or biomaRt.")
    }
}

# Homo sapiens - Humans
h_sapiens <- function(input_id, input_dataset, output_ids, output_species, tool, ortholog_database) { # input data funneling in from the cross species conversion function
    if (tool == "biomaRt") { # if the user has choose the tool biomart this it the path the codes follows
        Filter_BO_HS <- switch(output_species,
            c_elegans = {
                output_data <- biomaRt(input_id, input_dataset, output_ids,
                    input_species = "h_sapiens", output_species, 6239, "hsapiens_gene_ensembl",
                    "celegans_gene_ensembl", ortholog_database
                )
            },
            d_melanogaster = {
                output_data <- biomaRt(input_id, input_dataset, output_ids,
                    input_species = "h_sapiens", output_species, 7227, "hsapiens_gene_ensembl",
                    "dmelanogaster_gene_ensembl", ortholog_database
                )
            },
            m_musculus = {
                output_data <- biomaRt(input_id, input_dataset, output_ids,
                    input_species = "h_sapiens", output_species, 10090, "hsapiens_gene_ensembl",
                    "mmusculus_gene_ensembl", ortholog_database
                )
            },
            d_rerio = {
                output_data <- biomaRt(input_id, input_dataset, output_ids,
                    input_species = "h_sapiens", output_species, 7955, "hsapiens_gene_ensembl",
                    "drerio_gene_ensembl", ortholog_database
                )
            },
            h_sapiens = {
                output_data <- biomaRt(input_id, input_dataset, output_ids,
                    input_species = "h_sapiens", output_species, 9606, "hsapiens_gene_ensembl",
                    "hsapiens_gene_ensembl", ortholog_database
                )
            },
            r_norvegicus = {
                output_data <- biomaRt(input_id, input_dataset, output_ids,
                    input_species = "h_sapiens", output_species, 10116, "hsapiens_gene_ensembl",
                    "rnorvegicus_gene_ensembl", ortholog_database
                )
            },
            stop("Error: Invalid o_species in CoSIAn Object. Make sure the species in the o_species slot is an avalible model
                                    organism and is in the correct format.")
        )
        return(output_data)
    } else if (tool == "annotationDBI") { # code follows this path if the user chooses annotationDBI as their tool of choice
        Filter_AO_HS <- switch(output_species,
            h_sapiens = {
                output_data <- annotationDBI(input_id, input_dataset, output_ids,
                    input_species = "h_sapiens", output_species, 9606, org.Hs.eg.db::org.Hs.eg.db,
                    org.Hs.eg.db::org.Hs.eg.db, ortholog_database
                )
            },
            m_musculus = {
                output_data <- annotationDBI(input_id, input_dataset, output_ids,
                    input_species = "h_sapiens", output_species, 10090, org.Hs.eg.db::org.Hs.eg.db,
                    org.Mm.eg.db::org.Mm.eg.db, ortholog_database
                )
            },
            d_melanogaster = {
                output_data <- annotationDBI(input_id, input_dataset, output_ids,
                    input_species = "h_sapiens", output_species, 7227, org.Hs.eg.db::org.Hs.eg.db,
                    org.Dm.eg.db::org.Dm.eg.db, ortholog_database
                )
            },
            d_rerio = {
                output_data <- annotationDBI(input_id, input_dataset, output_ids,
                    input_species = "h_sapiens", output_species, 7955, org.Hs.eg.db::org.Hs.eg.db,
                    org.Dr.eg.db::org.Dr.eg.db, ortholog_database
                )
            },
            c_elegans = {
                output_data <- annotationDBI(input_id, input_dataset, output_ids,
                    input_species = "h_sapiens", output_species, 6239, org.Hs.eg.db::org.Hs.eg.db,
                    org.Ce.eg.db::org.Ce.eg.db, ortholog_database
                )
            },
            r_norvegicus = {
                output_data <- annotationDBI(input_id, input_dataset, output_ids,
                    input_species = "h_sapiens", output_species, 10116, org.Hs.eg.db::org.Hs.eg.db,
                    org.Rn.eg.db::org.Rn.eg.db, ortholog_database
                )
            },
            stop("Error: Invalid o_species in CoSIAn Object. Make sure the species in the o_species slot is an avalible model
                                    organism and is in the correct format.")
        )
        return(output_data)
    } else { # code follows this path if the tool that was inputted into the system does not match.
        stop("Error: Invalid tool in CoSIAn Object. Make sure the tool is either annotationDBI or biomaRt.")
    }
}

# Mus musculus - Mouse

m_musculus <- function(input_id, input_dataset, output_ids, output_species, tool, ortholog_database) { # input data funneling in from the cross species conversion function
    if (tool == "biomaRt") { # if the user has choose the tool biomart this it the path the codes follows
        Filter_BO_MM <- switch(output_species,
            c_elegans = {
                output_data <- biomaRt(input_id, input_dataset, output_ids,
                    input_species = "m_musculus", output_species, 6239, "mmusculus_gene_ensembl",
                    "celegans_gene_ensembl", ortholog_database
                )
            },
            d_melanogaster = {
                output_data <- biomaRt(input_id, input_dataset, output_ids,
                    input_species = "m_musculus", output_species, 7227, "mmusculus_gene_ensembl",
                    "dmelanogaster_gene_ensembl", ortholog_database
                )
            },
            m_musculus = {
                output_data <- biomaRt(input_id, input_dataset, output_ids,
                    input_species = "m_musculus", output_species, 10090, "mmusculus_gene_ensembl",
                    "mmusculus_gene_ensembl", ortholog_database
                )
            },
            d_rerio = {
                output_data <- biomaRt(input_id, input_dataset, output_ids,
                    input_species = "m_musculus", output_species, 7955, "mmusculus_gene_ensembl",
                    "drerio_gene_ensembl", ortholog_database
                )
            },
            h_sapiens = {
                output_data <- biomaRt(input_id, input_dataset, output_ids,
                    input_species = "m_musculus", output_species, 9606, "mmusculus_gene_ensembl",
                    "hsapiens_gene_ensembl", ortholog_database
                )
            },
            r_norvegicus = {
                output_data <- biomaRt(input_id, input_dataset, output_ids,
                    input_species = "m_musculus", output_species, 10116, "mmusculus_gene_ensembl",
                    "rnorvegicus_gene_ensembl", ortholog_database
                )
            },
            stop("Error: Invalid o_species in CoSIAn Object. Make sure the species in the o_species slot is an avalible model
                                    organism and is in the correct format.")
        )
    } else if (tool == "annotationDBI") { # code follows this path if the user chooses annotationDBI as their tool of choice
        Filter_AO_MM <- switch(output_species,
            h_sapiens = {
                output_data <- annotationDBI(input_id, input_dataset, output_ids,
                    input_species = "m_musculus", output_species, 9606, org.Mm.eg.db::org.Mm.eg.db,
                    org.Hs.eg.db::org.Hs.eg.db, ortholog_database
                )
            },
            m_musculus = {
                output_data <- annotationDBI(input_id, input_dataset, output_ids,
                    input_species = "m_musculus", output_species, 10090, org.Mm.eg.db::org.Mm.eg.db,
                    org.Mm.eg.db::org.Mm.eg.db, ortholog_database
                )
            },
            d_melanogaster = {
                output_data <- annotationDBI(input_id, input_dataset, output_ids,
                    input_species = "m_musculus", output_species, 7227, org.Mm.eg.db::org.Mm.eg.db,
                    org.Dm.eg.db::org.Dm.eg.db, ortholog_database
                )
            },
            d_rerio = {
                output_data <- annotationDBI(input_id, input_dataset, output_ids,
                    input_species = "m_musculus", output_species, 7955, org.Mm.eg.db::org.Mm.eg.db,
                    org.Dr.eg.db::org.Dr.eg.db, ortholog_database
                )
            },
            c_elegans = {
                output_data <- annotationDBI(input_id, input_dataset, output_ids,
                    input_species = "m_musculus", output_species, 6239, org.Mm.eg.db::org.Mm.eg.db,
                    org.Ce.eg.db::org.Ce.eg.db, ortholog_database
                )
            },
            r_norvegicus = {
                output_data <- annotationDBI(input_id, input_dataset, output_ids,
                    input_species = "m_musculus", output_species, 10116, org.Mm.eg.db::org.Mm.eg.db,
                    org.Rn.eg.db::org.Rn.eg.db, ortholog_database
                )
            },
            stop("Error: Invalid o_species in CoSIAn Object. Make sure the species in the o_species slot is an avalible model
                                    organism and is in the correct format.")
        )
    } else { # code follows this path if the tool that was inputted into the system does not match.
        stop("Error: Invalid tool in CoSIAn Object. Make sure the tool is either annotationDBI or biomaRt.")
    }
}

# Rattus norvegicus- Rat

r_norvegicus <- function(input_id, input_dataset, output_ids, output_species, tool, ortholog_database) { # input data funneling in from the cross species conversion function
    if (tool == "biomaRt") { # if the user has choose the tool biomart this it the path the codes follows
        Filter_BO_RN <- switch(output_species,
            c_elegans = {
                output_data <- biomaRt(input_id, input_dataset, output_ids,
                    input_species = "r_norvegicus", output_species, 6239, "rnorvegicus_gene_ensembl",
                    "celegans_gene_ensembl", ortholog_database
                )
            },
            d_melanogaster = {
                output_data <- biomaRt(input_id, input_dataset, output_ids,
                    input_species = "r_norvegicus", output_species, 7227, "rnorvegicus_gene_ensembl",
                    "dmelanogaster_gene_ensembl", ortholog_database
                )
            },
            m_musculus = {
                output_data <- biomaRt(input_id, input_dataset, output_ids,
                    input_species = "r_norvegicus", output_species, 10090, "rnorvegicus_gene_ensembl",
                    "mmusculus_gene_ensembl", ortholog_database
                )
            },
            d_rerio = {
                output_data <- biomaRt(input_id, input_dataset, output_ids,
                    input_species = "r_norvegicus", output_species, 7955, "rnorvegicus_gene_ensembl",
                    "drerio_gene_ensembl", ortholog_database
                )
            },
            h_sapiens = {
                output_data <- biomaRt(input_id, input_dataset, output_ids,
                    input_species = "r_norvegicus", output_species, 9606, "rnorvegicus_gene_ensembl",
                    "hsapiens_gene_ensembl", ortholog_database
                )
            },
            r_norvegicus = {
                output_data <- biomaRt(input_id, input_dataset, output_ids,
                    input_species = "r_norvegicus", output_species, 10116, "rnorvegicus_gene_ensembl",
                    "rnorvegicus_gene_ensembl", ortholog_database
                )
            },
            stop("Error: Invalid o_species in CoSIAn Object. Make sure the species in the o_species slot is an avalible model
                                    organism and is in the correct format.")
        )
    } else if (tool == "annotationDBI") { # code follows this path if the user chooses annotationDBI as their tool of choice
        Filter_AO_RN <- switch(output_species,
            h_sapiens = {
                output_data <- annotationDBI(input_id, input_dataset, output_ids,
                    input_species = "r_norvegicus", output_species, 9606, org.Rn.eg.db::org.Rn.eg.db,
                    org.Hs.eg.db::org.Hs.eg.db, ortholog_database
                )
            },
            m_musculus = {
                output_data <- annotationDBI(input_id, input_dataset, output_ids,
                    input_species = "r_norvegicus", output_species, 10090, org.Rn.eg.db::org.Rn.eg.db,
                    org.Mm.eg.db::org.Mm.eg.db, ortholog_database
                )
            },
            d_melanogaster = {
                output_data <- annotationDBI(input_id, input_dataset, output_ids,
                    input_species = "r_norvegicus", output_species, 7227, org.Rn.eg.db::org.Rn.eg.db,
                    org.Dm.eg.db::org.Dm.eg.db, ortholog_database
                )
            },
            d_rerio = {
                output_data <- annotationDBI(input_id, input_dataset, output_ids,
                    input_species = "r_norvegicus", output_species, 7955, org.Rn.eg.db::org.Rn.eg.db,
                    org.Dr.eg.db::org.Dr.eg.db, ortholog_database
                )
            },
            c_elegans = {
                output_data <- annotationDBI(input_id, input_dataset, output_ids,
                    input_species = "r_norvegicus", output_species, 6239, org.Rn.eg.db::org.Rn.eg.db,
                    org.Ce.eg.db::org.Ce.eg.db, ortholog_database
                )
            },
            r_norvegicus = {
                output_data <- annotationDBI(input_id, input_dataset, output_ids,
                    input_species = "r_norvegicus", output_species, 10116, org.Rn.eg.db::org.Rn.eg.db,
                    org.Rn.eg.db::org.Rn.eg.db, ortholog_database
                )
            },
            stop("Error: Invalid o_species in CoSIAn Object. Make sure the species in the o_species slot is an avalible model
                                    organism and is in the correct format.")
        )
    } else { # code follows this path if the tool that was inputted into the system does not match.
        stop("Error: Invalid tool in CoSIAn Object. Make sure the tool is either annotationDBI or biomaRt.")
    }
}

# Tool functions

# annotationDBI

annotationDBI <- function(input_id, input_dataset, output_ids, input_species, output_species, species_number, input_org, output_org, ortholog_database) {
    # This function uses annotationDBI to convert between gene identifier
    # Switch function below is used to transition the id names into the proper format for using annotationdbi
    input_dataset <- as.character(input_dataset)
    ID_SWITCH <- Vectorize(vectorize.args = "ids", FUN = function(ids) {
        switch(as.character(ids),
            Entrez_id = "ENTREZID",
            Ensembl_id = "ENSEMBL",
            Symbol = "SYMBOL",
            stop("Error: Invalid input_id or output_ids. The acceptable ids for annotationDBI are Entrez_id, Ensembl_id, Ensembl_id_version, Symbol, and Gene_name.")
        )
    })
    input_id <- ID_SWITCH(ids = input_id)
    output_ids <- ID_SWITCH(ids = output_ids)
    input_id <- as.character(input_id)
    output_ids <- as.character(output_ids)
    if (output_species == input_species) {
        annotationDbi_data <- AnnotationDbi::select(input_org, keys = input_dataset, columns = output_ids, keytype = input_id) # run the code through annotationDBI
        colnames(annotationDbi_data)[which(names(annotationDbi_data) == "ENSEMBL")] <- paste(output_species, "ensembl_id", sep = "_")
        colnames(annotationDbi_data)[which(names(annotationDbi_data) == "ENTREZID")] <- paste(output_species, "entrez_id", sep = "_")
        colnames(annotationDbi_data)[which(names(annotationDbi_data) == "ENSEMBLIDVERSION")] <- paste(output_species, "ensembl_id_version", sep = "_")
        colnames(annotationDbi_data)[which(names(annotationDbi_data) == "GENENAME")] <- paste(output_species, "gene_name", sep = "_")
        colnames(annotationDbi_data)[which(names(annotationDbi_data) == "SYMBOL")] <- paste(output_species, "symbol", sep = "_")
        return(annotationDbi_data)
    } else {
        output_data <- AnnotationDbi::select(input_org, keys = input_dataset, columns = c("ENTREZID"), keytype = input_id) # runs the gene through annotationdbi to get EntrezID conversion in order to conduct ortholog mapping
        output_data <- na.omit(output_data) # omits the NA values
        id <- homolog(output_data$ENTREZID, species_number, ortholog_database) # send the entrezid ids through the homolog conversion
        names(output_data)[names(output_data) == "ENTREZID"] <- "species_one" # make the EntrezID that were used for the conversion into species one
        merged_data <- merge.data.frame(data.frame(output_data), data.frame(id), by = "species_one") # merge the two dataframes with the original ids and the ortholog ids together
        colnames(merged_data)[which(names(merged_data) == "species_one")] <- paste(input_species, "entrez_id", sep = "_") # change column name so that they are representative of the input species
        if (input_id != "ENTREZID") {
            entrez_id_col_name_input <- paste(input_species, "entrez_id", sep = "_")
            merged_data <- merged_data %>% dplyr::select(-tidyselect::all_of(entrez_id_col_name_input))
        }
        colnames(merged_data)[which(names(merged_data) == "ENTREZID")] <- paste(input_species, "entrez_id", sep = "_") # changes name to more formal names
        colnames(merged_data)[which(names(merged_data) == "ENSEMBL")] <- paste(input_species, "ensembl_id", sep = "_") # changes name to more formal names
        colnames(merged_data)[which(names(merged_data) == "SYMBOL")] <- paste(input_species, "symbol", sep = "_") # changes name to more formal names
        colnames(merged_data)[which(names(merged_data) == "ENSEMBLIDVERSION")] <- paste(input_species, "ensembl_id_version", sep = "_") # changes name to more formal names
        colnames(merged_data)[which(names(merged_data) == "GENENAME")] <- paste(input_species, "gene_name", sep = "_") # changes name to more formal names
        ortho <- AnnotationDbi::select(output_org,
            keys = as.character(as.matrix(merged_data$species_two)), columns = output_ids,
            keytype = "ENTREZID"
        ) # now that you have the entrezids of the converted species we can convert to different gene identifiers
        non_na <- which(is.na(ortho$ENTREZID) == FALSE) # Determine the indices for the non-NA genes
        ortho <- ortho[non_na, ] # Return only the genes with annotations using indices
        non_dups <- which(duplicated(ortho$ENTREZID) == FALSE) # Determine the indices for the non-duplicated genes
        ortho <- ortho[non_dups, ] # Return only the non-duplicated genes using indices
        ortho <- data.frame(ortho) # make the ortholog conversion into a dataframe
        names(ortho)[names(ortho) == "ENTREZID"] <- "species_two" # rename the species into species two in order to merge the orthologs and the original data
        merged_data_ortholog <- merge.data.frame(data.frame(ortho), data.frame(merged_data), by = "species_two") # merge
        colnames(merged_data_ortholog)[names(merged_data_ortholog) == "species_two"] <- paste(output_species, "entrez_id", sep = "_") # changes name to more formal names
        if ("ENTREZID" %in% output_ids == FALSE) {
            entrez_id_col_output <- paste(output_species, "entrez_id", sep = "_")
            merged_data_ortholog <- merged_data_ortholog %>% dplyr::select(-tidyselect::all_of(entrez_id_col_output))
        }
        colnames(merged_data_ortholog)[which(names(merged_data_ortholog) == "ENTREZID")] <- paste(output_species, "entrez_id", sep = "_") # changes name to more formal names
        colnames(merged_data_ortholog)[which(names(merged_data_ortholog) == "ENSEMBL")] <- paste(output_species, "ensembl_id", sep = "_") # changes name to more formal names
        colnames(merged_data_ortholog)[which(names(merged_data_ortholog) == "SYMBOL")] <- paste(output_species, "symbol", sep = "_") # changes name to more formal names
        colnames(merged_data_ortholog)[which(names(merged_data_ortholog) == "ENSEMBLIDVERSION")] <- paste(output_species, "ensembl_id_version", sep = "_") # changes name to more formal names
        colnames(merged_data_ortholog)[which(names(merged_data_ortholog) == "GENENAME")] <- paste(output_species, "gene_name", sep = "_") # changes name to more formal names
        merged_data_ortholog <- merged_data_ortholog[!duplicated(as.list(merged_data_ortholog))] # remove duplicates
        return(merged_data_ortholog) # works
    }
}

# biomaRt

biomaRt <- function(input_id, input_dataset, output_ids, input_species, output_species, species_number, species_dataset, output_species_dataset,
                    ortholog_database) {
    species_dataset <- as.character(species_dataset)
    # Make sure ids are in the class character
    input_id <- as.character(input_id)
    output_ids <- as.character(output_ids)
    # create a switch function that renames the user ids into the designated format for bioconductor
    ID_SWITCH <- Vectorize(vectorize.args = "ids", FUN = function(ids) {
        switch(as.character(ids),
            Entrez_id = "entrezgene_id",
            Ensembl_id = "ensembl_gene_id",
            Symbol = "external_gene_name",
            stop("Error: Invalid input_id or output_ids. The acceptable ids for biomaRt is Entrez_id, Ensembl_id, and Symbol")
        )
    })
    # set the id names to their new formats
    input_id <- ID_SWITCH(ids = input_id)
    output_ids <- ID_SWITCH(ids = output_ids)
    if (input_species == output_species) {
        # goes through this path if the input and output species are the same
        mart <- biomaRt::useMart("ensembl", dataset = species_dataset) # pulls the biomaRt object for the species species that has been choosen
        attributes <- c(output_ids, input_id) # sets the attributes that the user wants makes sure to set both the input and output values
        attributes <- as.character(attributes)
        filters <- input_id # sets the filters as the input vales
        filters <- as.character(filters)
        input_dataset <- as.character(input_dataset)
        output_data <- biomaRt::getBM(attributes = attributes, filters = filters, values = input_dataset, mart = mart, uniqueRows = TRUE, bmHeader = FALSE) # run conversion through biomaRt
        colnames(output_data)[which(names(output_data) == "ensembl_gene_id")] <- paste(input_species, "ensembl_id", sep = "_")
        colnames(output_data)[which(names(output_data) == "entrezgene_id")] <- paste(input_species, "entrez_id", sep = "_")
        colnames(output_data)[which(names(output_data) == "external_gene_name")] <- paste(input_species, "symbol", sep = "_")
        output_data <- output_data %>% dplyr::select(-contains("."))
        return(output_data) # return the biomaRt output
    } else { # goes through this path if the input and output species are different
        mart <- biomaRt::useMart("ensembl", dataset = species_dataset) # sets up biomart for species input
        attributes <- c("entrezgene_id", input_id) # sets input ids and entrezids as attributes (output values)
        filters <- input_id # set the input ids as the filter (input values)
        output_data <- biomaRt::getBM(attributes = attributes, filters = filters, values = input_dataset, mart = mart, uniqueRows = TRUE, bmHeader = FALSE) # run the convesion through biomaRt
        output_data <- na.omit(output_data) # omit the NA values
        id <- homolog(output_data$entrezgene_id, species_number, ortholog_database) # now that we have the entrezids for the input species we can run it through the homolog function and then get the entrezids for the other species
        id <- data.frame(id) # this makes sure that our homologous gene list is in a dataframe
        names(output_data)[names(output_data) == "entrezgene_id"] <- "species_one" # change the name so we can merge the original and the homologous dataframes together
        # change the name of the input id
        colnames(output_data)[which(names(output_data) == "ensembl_gene_id")] <- paste(input_species, "ensembl_id", sep = "_")
        colnames(output_data)[which(names(output_data) == "ensembl_gene_id_version")] <- paste(input_species, "ensembl_id_version", sep = "_")
        colnames(output_data)[which(names(output_data) == "external_gene_name")] <- paste(input_species, "symbol", sep = "_")
        merged_data <- merge.data.frame(data.frame(output_data), data.frame(id), by = "species_one") # merge the two dataframes
        colnames(merged_data)[which(names(merged_data) == "species_one")] <- paste(input_species, "entrez_id", sep = "_") # rename to a more formal name
        marts <- biomaRt::useMart("ensembl", dataset = output_species_dataset) # set the biomart species to the new species
        ortho <- biomaRt::getBM(
            attributes = c(output_ids, "entrezgene_id"), filters = "entrezgene_id", values = id$species_two, mart = marts,
            uniqueRows = TRUE, bmHeader = FALSE
        ) # run the homologous entrezids to the other gene identifiers
        non_na <- which(is.na(ortho$entrezgene_id) == FALSE) # pick non na values from ortholog conversion
        ortho <- ortho[non_na, ] # subset only nonna values
        non_dups <- which(duplicated(ortho$entrezgene_id) == FALSE) # remove duplicates
        ortho <- ortho[non_dups, ] # subset only duplicates
        ortho <- data.frame(ortho) # set the ortho variable as a dataframe
        names(ortho)[names(ortho) == "entrezgene_id"] <- "species_two" # rename the entrezid as species two
        merged_species_data <- merge.data.frame(data.frame(ortho), data.frame(merged_data), by = "species_two") # merge the ortholog dataframe and the original data
        merged_species_data <- merged_species_data[, !names(merged_species_data) %in% c("entrezgene_id.1")]
        names(merged_species_data)[names(merged_species_data) == "species_two"] <- paste(output_species, "entrez_id", sep = "_") # clean up names
        colnames(merged_species_data)[which(names(merged_species_data) == "entrezgene_id")] <- paste(output_species, "entrez_id", sep = "_") # clean up names
        colnames(merged_species_data)[which(names(merged_species_data) == "ensembl_gene_id")] <- paste(output_species, "ensembl_id", sep = "_") # clean up names
        colnames(merged_species_data)[which(names(merged_species_data) == "external_gene_name")] <- paste(output_species, "symbol", sep = "_") # clean up names
        merged_species_data <- merged_species_data[!duplicated(as.list(merged_species_data))] # remove duplicates
        if (input_id != "entrezgene_id") {
            entrez_id_col_name_input <- paste(input_species, "entrez_id", sep = "_")
            merged_species_data <- merged_species_data %>% dplyr::select(-tidyselect::all_of(entrez_id_col_name_input))
        }
        if ("entrezgene_id" %in% output_ids == FALSE) {
            entrez_id_col_output <- paste(output_species, "entrez_id", sep = "_")
            merged_species_data <- merged_species_data %>% dplyr::select(-tidyselect::all_of(entrez_id_col_output))
        }
        return(merged_species_data) # return the final table
    }
}


# homolog files

homolog <- function(entrez_data, species_number, ortholog_database) {
    if (ortholog_database == "HomoloGene") {
        myGenes <- as.character(entrez_data)
        species_one <- entrez_data
        homologene <- homologene::homologeneData2
        data <- annotationTools::getHOMOLOG(myGenes, species_number, homologene, noIDsymbol = NA, clusterCol = 1, speciesCol = 3, idCol = 4)
        data <- as.list(data)
        species_two <- as.character(data)
        homologs <- as.data.frame(species_two)
        homologs <- merge(data.frame(homologs, row.names = NULL), data.frame(species_one, row.names = NULL), by = 0, all = TRUE)[-1]
        manipulated_homologs <- data.frame(homologs)
        x <- nrow(manipulated_homologs)
        for (i in seq_len(x)) {
            ## print(i)
            char <- manipulated_homologs[i, 1]
            ## print(char)
            if ((char == "NA") == TRUE) {
                manipulated_homologs[i, 1] <- NA
            }
        }
        species_two_entrez_ID <- na.omit(manipulated_homologs)
        return(species_two_entrez_ID)
    }
    if (ortholog_database == "NCBIOrtho") {
        myGenes <- as.vector(entrez_data)
        species_one <- entrez_data
        NCBIOrtho <- readr::read_table("https://ftp.ncbi.nih.gov/gene/DATA/gene_orthologs.gz", col_types = "ccccc", show_col_types = TRUE)
        NCBIOrtho <- as.data.frame(NCBIOrtho)
        data <- annotationTools::getHOMOLOG(myGenes, species_number, NCBIOrtho, tableType = "gene_orthologs")
        data <- as.list(data)
        species_two <- as.character(data)
        homologs <- as.data.frame(species_two)
        homologs <- merge(data.frame(homologs, row.names = NULL), data.frame(species_one, row.names = NULL), by = 0, all = TRUE)[-1]
        manipulated_homologs <- data.frame(homologs)
        x <- nrow(manipulated_homologs)
        for (i in seq_len(x)) {
            ## print(i)
            char <- manipulated_homologs[i, 1]
            ## print(char)
            if ((char == "NA") == TRUE) {
                manipulated_homologs[i, 1] <- NA
            }
        }
        species_two_entrez_ID <- na.omit(manipulated_homologs)
        return(species_two_entrez_ID)
    }
}

# remove_version_number (only for inputs with Ensembl_id_version)

remove_version_numbers <- function(input_dataset, species) {
    if (species == "m_musculus") { # ENSMUSG00000032855.7
        input_dataset <- data.frame(input_dataset)
        colnames(input_dataset)[1] <- "ENSEMBL"
        x <- nrow(input_dataset)
        for (i in seq_len(x)) {
            char <- input_dataset[i, 1]
            char1 <- substr(char, start = 1, stop = 18)
            input_dataset[i, 1] <- char1
        }
        colnames(input_dataset)[which(names(input_dataset) == "ENSEMBL")] <- "input"
        return(input_dataset$input)
    }
    if (species == "h_sapiens") { # ENSG00000008710.20
        input_dataset <- data.frame(input_dataset)
        colnames(input_dataset)[1] <- "ENSEMBL"
        x <- nrow(input_dataset)
        for (i in seq_len(x)) {
            char <- input_dataset[i, 1]
            char1 <- substr(char, start = 1, stop = 15)
            input_dataset[i, 1] <- char1
        }
        colnames(input_dataset)[which(names(input_dataset) == "ENSEMBL")] <- "input"
        return(input_dataset$input)
    }

    if (species == "d_rerio") { # ENSDARG00000105344.2 (format)
        input_dataset <- data.frame(input_dataset)
        colnames(input_dataset)[1] <- "ENSEMBL"
        x <- nrow(input_dataset)
        for (i in seq_len(x)) {
            char <- input_dataset[i, 1]
            char1 <- substr(char, start = 1, stop = 18)
            input_dataset[i, 1] <- char1
        }
        colnames(input_dataset)[which(names(input_dataset) == "ENSEMBL")] <- "input"
        return(input_dataset$input)
    }

    if (species == "r_norvegicus") { # ENSRNOG00000010771.8 (format)
        input_dataset <- data.frame(input_dataset)
        colnames(input_dataset)[1] <- "ENSEMBL"
        x <- nrow(input_dataset)
        for (i in seq_len(x)) {
            char <- input_dataset[i, 1]
            char1 <- substr(char, start = 1, stop = 18)
            input_dataset[i, 1] <- char1
        }
        colnames(input_dataset)[which(names(input_dataset) == "ENSEMBL")] <- "input"
        return(input_dataset$input)
    }

    if (species == "c_elegans") {
        input_dataset <- data.frame(input_dataset)
        return(input_dataset$input)
    }
    if (species == "d_melanogaster") {
        input_dataset <- data.frame(input_dataset)
        return(input_dataset$input)
    } else {
        stop("Error. Invalid species. Make sure it matches the proper format.")
    }
}
