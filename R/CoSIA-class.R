#' CoSIAn Class
#' @name CoSIAn
#' @rdname CoSIA-class
#' @aliases CoSIA-class
#' @slot gene_set character. A singular gene or a list of genes.
#' @slot i_species character. The species corresponding to gene_set.
#' @slot i_id character. The type of id corresponding to gene_set.
#' @slot o_species character. The singular or list of species that the gene set is being converted to.
#' @slot o_ids character. The singular or list of ide types that the gene set is being converted to.
#' @slot mapping_tool character. The mapping tool, BiomaRt or annotationDBI, being used to map the conversion between IDs.
#' @slot ortholog_database character. The ortholog database, HomoloGene or NCBIOrtho, being used to map the conversion between species.
#' @slot converted_id data frame. Output of getConversion.
#' @slot map_tissues character. A list of tissues being investigated.
#' @slot map_species character. A list of species being investigated. 
#' @slot gex data frame. Output of gene expression data.
#' @slot metric_type character. A list of possible metric the user wants to calculate.
#' @slot metric data frame. Output of gene expression metrics data.
#' @exportClass CoSIAn


setClass("CoSIAn",
         slots = c(
           gene_set = "character",
           i_species = "character",
           i_id = "character",
           o_species = "character",
           o_ids = "character",
           mapping_tool = "character",
           ortholog_database = "character",
           converted_id = "data.frame",
           map_tissues = "character",
           map_species = "character",
           gex = "data.frame",
           metric_type = "character",
           metric = "data.frame"
         ),
         prototype = list(
           gene_set = character(0),
           i_species = character(0),
           i_id = character(0),
           o_species = character(0),
           o_ids = character(0),
           mapping_tool = "annotationDBI", #AnnotationDBI is the default
           ortholog_database = "HomoloGene", #HomoloGene is the default
           converted_id = data.frame(0),
           map_tissues = character(0),
           map_species = character(0),
           gex = data.frame(0),
           metric_type = character(0),
           metric = data.frame(0)
         )
)

## Constructor (COPY OVER THE DESCRIPTION AFTER APPROVAL FROM CLASS)
#' @description The \code{CoSIAn} constructor creates a \code{CoSIAn} object from character vector(s).
#' @rdname CoSIA-class
#' @param gene_set 
#' @param i_species 
#' @param i_id 
#' @param o_species 
#' @param o_ids 
#' @param mapping_tool 
#' @param ortholog_database 
#' @param map_tissues 
#' @param map_species
#' @param metric_type 
#' @return An S4 \code{CoSIAn} object with character vector(s) as slots.
#' @export
#' @examples

CoSIAn <- function(gene_set, i_species, i_id, o_species, o_ids, mapping_tool="annotationDBI", ortholog_database= "HomoloGene", map_tissues, map_species,metric_type) {
  gene_set<- as.character(gene_set)
  i_species<-as.character(i_species)
  i_id<-as.character(i_id)
  o_species<- as.character(o_species)
  o_ids<- as.character(o_ids)
  map_tissues<- as.character(map_tissues)
  map_species<- as.character(map_species)
  metric_type<- as.character(metric_type)
  converted_id<- data.frame(0)
  gex<- data.frame(0)
  metric <- data.frame(0)
  new("CoSIAn", gene_set=gene_set, i_species=i_species, i_id=i_id, o_species=o_species, o_ids=o_ids, mapping_tool=mapping_tool, 
      ortholog_database=ortholog_database,converted_id = converted_id, map_tissues=map_tissues, map_species=map_species, gex = gex, metric_type=metric_type,metric = metric)
}

## Validity of the CoSIAn Class
setValidity("CoSIAn", function(object) {
  if(length(object@gene_set)< 1){"@gene_set needs to be filled"}
  if(length(object@i_species) != 1){"@i_species needs to be a length of 1"}
  if(length(object@i_id) != 1){"@i_id needs to be a length of 1"}
  if (length(object@i_species) != length(i_id)) { "@i_species and @i_id must be the same length" }
  if(length(object@o_species) < 1){ "@o_species needs to be filled" }
  if(length(object@o_ids) < 1){ "@o_ids needs to be filled" }
  if (length(object@o_species) != length(object@o_ids)) {"@o_species and @o_ids must be the same length" }
  if(length(object@mapping_tool) != 1){"@mapping_tool needs to be a length of 1"}
  if(length(object@ortholog_database) != 1){ "@ortholog_database needs to be a length of 1" }
  if(length(object@map_tissues)< 1){"@map_tissues needs to be filled"}
  if(length(object@map_species)< 1){"@map_species needs to be filled"}
  if(length(object@metric_type)< 1){"@metric_type needs to be filled"}
  else {TRUE}
})

