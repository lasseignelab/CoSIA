# #CoSIA Utility Methods
# 
# #Find the Available/Common Tissues in CoSIA
# 
# species<-c("h_sapiens","m_musculus")
# 
# setGeneric("getTissues", function(species) standardGeneric("getTissues"))
# #CoSIAn getGEx
# setMethod("getTissues", function(species) {
#   Species_SWITCH <- Vectorize(vectorize.args = "species", FUN = function(species) {
#     switch(as.character(species), 
#            h_sapiens ="Homo_sapiens",
#            m_musculus = "Mus_musculus",
#            r_norvegicus = "Rattus_norvegicus",
#            d_rerio = "Danio_rerio",
#            d_melanogaster = "Drosophila_melanogaster",
#            c_elegans = "Caenorhabditis_elegans",
#            stop("Error: Invalid species. Make sure you have choosen an avaliable species in CoSIA")
#     )})
#   species<-Species_SWITCH(species)
#   #load the EH data here
#   filter_species<-Experimental_Hub_File %>% dplyr::select(Anatomical_entity_name,Species)
#   filter_species<- unique(filter_species)
#   filter_species <- dplyr::filter(filter_species,Species %in% species)
#   common_tissue<- filter_species[filter_species$Anatomical_entity_name == Reduce(intersect, split(filter_species$Species, filter_species$Anatomical_entity_name)),]
# 
# 
#   
#   
# })