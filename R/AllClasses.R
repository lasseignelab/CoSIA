

setClass("CoSIA", representation("VIRTUAL"))# virtual class
setClass("CosiaAnnotate",contains="CoSIA", 
         slots = c(
           input = "character",
           input_species = "character",
           input_id = "character",
           output_ids = "character",
           output_species = "character",
           tool = "character",
           ortholog_database = "character"
         ),
         prototype = list(
           input = character(0),
           input_species = character(0),
           input_id = character(0),
           output_ids = character(0),
           output_species = character(0),
           tool =character(0),
           ortholog_database = character(0)
         )
)
setClass("CosiaExpressSpecies", contains="CoSIA", 
         slots=
           c(list_of_ensembl_ids = "character",
             list_of_respective_species = "character",
             single_tissue = "character",
             pathToData = "character", 
             dataframe = "data.frame"
           ),
         prototype = list(
           list_of_ensembl_ids = character(0),
           list_of_respective_species = character(0),
           single_tissue = character(0),
           pathToData = character(0),
           dataframe = data.frame(0)
         )
)
setClass("CosiaExpressTissue", contains="CoSIA", 
         slots=
           c(single_gene = "character",
             gene_species = "character",
             tissues = "character",
             pathToData = "character",
             dataframe = "data.frame"
           ),
         prototype = list(
           single_gene = character(0),
           gene_species = character(0),
           tissues = character(0),
           pathToData = character(0),
           dataframe = data.frame(0)
         )
)