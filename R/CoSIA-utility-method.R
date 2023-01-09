#CoSIA Utility Functions

#Find the Available/Common Tissues in CoSIA for a given combination of species

#' getTissues
#'
#' @param species name of a species or multiple species that you want to get available tissue list for
#'
#' @return list of tissues that are common/available among the species or multiple species inputted
#' @export
#'
#' @examples
#' load("~/Desktop/EH_Data.RData")
#' tissue<-getTissues(c("m_musculus"))

getTissues<- function(species) {
  List_of_Tissues<-Experimental_Hub_File %>% dplyr::group_by(Anatomical_entity_name) %>% dplyr::summarise(Anatomical_entity_ID = unique(Anatomical_entity_ID), Species=unique(Species))
  Species_SWITCH <- Vectorize(vectorize.args = "species", FUN = function(species) {
    switch(as.character(species),
           h_sapiens ="Homo_sapiens",
           m_musculus = "Mus_musculus",
           r_norvegicus = "Rattus_norvegicus",
           d_rerio = "Danio_rerio",
           d_melanogaster = "Drosophila_melanogaster",
           c_elegans = "Caenorhabditis_elegans",
           stop("Error: Invalid species. Make sure you have choosen an avaliable species in CoSIA")
    )})
  species<-Species_SWITCH(species)
  #load the EH data here
  LOT<-dplyr::filter(List_of_Tissues, Species %in% species)
  LOT<-subset(LOT, select = -c(Species))
  L<-LOT %>% dplyr::summarise(Frequency = table(Anatomical_entity_name))
  frequency_value<-length(species)
  common_tissue<-dplyr::filter(L, Frequency == frequency_value)
  common_tissue<-subset(common_tissue, select = -c(Frequency))
  colnames(common_tissue)[which(names(common_tissue) == "Anatomical_entity_name")] <- 'Common_Anatomical_Entity_Name'
  return(common_tissue)
}