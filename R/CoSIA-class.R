setClass("CoSIA", representation("VIRTUAL"))  # virtual class

#' CoSIAn Class
#' @name CoSIAn
#' @rdname CoSIA-class
#' @slot gene_set character. A singular gene or a list of genes.
#' @slot i_species character. The species corresponding to gene_set.
#' @slot input_id character. The type of id corresponding to gene_set.
#' @slot o_species character. The singular or list of species that the gene set is being converted to.
#' @slot output_ids character. The singular or list of id types that the gene set is being converted to.
#' @slot mapping_tool character. The mapping tool, BiomaRt or annotationDBI, being used to map the conversion between IDs.
#' @slot ortholog_database character. The ortholog database, HomoloGene or NCBIOrtho, being used to map the conversion between species.
#' @slot converted_id data frame. Output of getConversion.
#' @slot map_tissues character. A list of tissues being investigated.
#' @slot map_species character. A list of species being investigated. 
#' @slot gex data frame. Output of gene expression data.
#' @slot metric_type character. A list of possible metric the user wants to calculate.
#' @slot metric data frame. Output of coefficient of variation of gene expression data.
#' @exportClass CoSIAn
setClass("CoSIAn",contains = "CoSIA",
         slots = c(
           gene_set = "character",
           i_species = "character",
           input_id = "character",
           o_species = "character",
           output_ids = "character",
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
           gene_set = NA_character_,
           i_species = NA_character_,
           input_id = NA_character_,
           o_species = NA_character_,
           output_ids = NA_character_,
           mapping_tool = "annotationDBI", #AnnotationDBI is the default
           ortholog_database = "HomoloGene", #HomoloGene is the default
           converted_id = data.frame(0),
           map_tissues = NA_character_,
           map_species = NA_character_,
           gex = data.frame(0),
           metric_type = NA_character_,
           metric = data.frame(0)
         )
)

## Constructor (COPY OVER THE DESCRIPTION AFTER APPROVAL FROM CLASS)
#' @description The \code{CoSIAn} constructor creates a \code{CoSIAn} object from character vector(s).
#' @rdname CoSIA-class
#' @param gene_set A singular gene or a list of genes.
#' @param i_species The species corresponding to gene_set.
#' @param input_id The type of id corresponding to gene_set.
#' @param o_species The singular or list of species that the gene set is being converted to.
#' @param output_ids The singular or list of id types that the gene set is being converted to.
#' @param mapping_tool The mapping tool, BiomaRt or annotationDBI, being used to map the conversion between IDs.
#' @param ortholog_database The ortholog database, HomoloGene or NCBIOrtho, being used to map the conversion between species.
#' @param map_tissues A list of tissues being investigated
#' @param map_species A list of species being investigated
#' @param metric_type A list of possible metric the user wants to calculate.
#' @return An S4 \code{CoSIAn} object with character vector(s) as slots.
#' @export
#' @examples
#' CoSIAn(gene_set = "PKD1",
#' i_species = "h_sapiens",
#' input_id = "Symbol", 
#' o_species = c("h_sapiens","m_musculus"), 
#' output_ids = "Ensembl_id", 
#' mapping_tool="annotationDBI", 
#' ortholog_database= "HomoloGene", 
#' map_tissues = "heart", 
#' map_species = "m_musculus", 
#' metric_type = "DS_Gene"
#' )



#constructor for user who is using all of the functions
CoSIAn <- function(gene_set, i_species, input_id, o_species, output_ids, mapping_tool="annotationDBI", ortholog_database= "HomoloGene", map_tissues, map_species, metric_type) {
  gene_set<- as.character(gene_set)
  i_species<-as.character(i_species)
  input_id<-as.character(input_id)
  o_species<- as.character(o_species)
  output_ids<- as.character(output_ids)
  map_tissues<- as.character(map_tissues)
  map_species<- as.character(map_species)
  metric_type<- as.character(metric_type)
  converted_id<- data.frame(0)
  gex<- data.frame(0)
  metric <- data.frame(0)
  methods::new("CoSIAn", gene_set=gene_set, i_species=i_species, input_id=input_id, o_species=o_species, output_ids=output_ids, mapping_tool=mapping_tool, 
      ortholog_database=ortholog_database,converted_id = converted_id, map_tissues=map_tissues, map_species=map_species, gex = gex, metric_type=metric_type,metric = metric)
}



## Validity of the CoSIAn Class
setValidity("CoSIAn", function(object) {
  if(length(object@gene_set)< 1){"@gene_set needs to be filled"}
  if(length(object@i_species) != 1){"@i_species needs to be a length of 1"}
  if(length(object@input_id) != 1){"@input_id needs to be a length of 1"}
  if (length(object@i_species) != length(input_id)) { "@i_species and @input_id must be the same length" }
  if(length(object@o_species) < 1){ "@o_species needs to be filled" }
  if(length(object@output_ids) < 1){ "@output_ids needs to be filled" }
  #if (length(object@o_species) != length(object@output_ids)) {"@o_species and @output_ids must be the same length" }
  if(length(object@mapping_tool) != 1){"@mapping_tool needs to be a length of 1"}
  if(length(object@ortholog_database) != 1){ "@ortholog_database needs to be a length of 1" }
  if(length(object@map_tissues)< 1){"@map_tissues needs to be filled"}
  if(length(object@map_species)< 1){"@map_species needs to be filled"}
  if(length(object@metric_type)< 1){"@metric_type needs to be filled"}
  else {TRUE}
})

