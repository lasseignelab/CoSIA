#Add the new CoSIA Class, Constructor, and Validity
setClass("CoSIAn", 
         slots = c(
           gene_set = "character",
           i_species = "character",
           i_id = "character",
           o_species = "character",
           o_ids = "character",
           tool = "character",
           ortholog_database = "character",
           converted_id = "dataframe",
           map_tissues = "character",
           map_species = "character",
           gex = "dataframe",
           
         )
