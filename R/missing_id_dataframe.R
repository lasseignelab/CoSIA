#' Missing ID Column
#'
#' @param input_data The originally character vector of column that contains the gene identifiers
#' @param output_data The annotated character vector or column that contains the same gene identifiers that were originally inputted
#'
#' @return a data frame that contains TRUE/FALSE for whether the gene identifier was mapped
#' @export
#'
#' @examples
#'
missing_id_column <- function(input_data, output_data) {
    original_dataframe <- input_data
    annotated_dataframe <- output_data
    if (length(original_dataframe) >= length(annotated_dataframe)) {
        tf <- original_dataframe %in% annotated_dataframe
        tf <- data.frame(tf)
        Mapped_IDs <- data.frame(original_dataframe)
        Mapped_IDs["IDs Mapped After Conversion"] <- tf
        names(Mapped_IDs)[names(Mapped_IDs) == "original_dataframe"] <- "Gene_Identifier"
        duplicated <- duplicated(Mapped_IDs$Gene_Identifier)
        Mapped_IDs["IDs_Mapped_Multiple_Times_After_Conversion"] <- duplicated
        return(Mapped_IDs)
    }
    if (length(original_dataframe) < length(annotated_dataframe)) {
        tf <- annotated_dataframe %in% original_dataframe
        tf <- data.frame(tf)
        Mapped_IDs <- data.frame(annotated_dataframe)
        Mapped_IDs["IDs Mapped After Conversion"] <- tf
        names(Mapped_IDs)[names(Mapped_IDs) == "annotated_dataframe"] <- "Gene_Identifier"
        duplicated <- duplicated(Mapped_IDs$Gene_Identifier)
        Mapped_IDs["IDs_Mapped_Multiple_Times_After_Conversion"] <- duplicated
        return(Mapped_IDs)
    }
}

