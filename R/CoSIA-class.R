setClass("CoSIA", representation("VIRTUAL"))  # virtual class

#' CoSIAn Class
#' @name CoSIAn
#' @aliases CoSIAn
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
#' @return CoSIAn object

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
#' Kidney_Genes<-CoSIAn(gene_set = c("ENSG00000008710","ENSG00000118762","ENSG00000152217"),
#' i_species = "h_sapiens",input_id = "Ensembl_id",o_species = c("d_melanogaster","m_musculus",
#' "h_sapiens", "d_rerio","c_elegans","r_norvegicus"),output_ids = c("Ensembl_id","Symbol"), 
#' mapping_tool = "annotationDBI",ortholog_database = "HomoloGene",map_tissues = "heart", 
#' map_species = c("m_musculus"),metric_type = "DS_Gene")


## Constructor for CoSIAn object
CoSIAn <- function(gene_set, i_species, input_id, o_species, output_ids, 
                   mapping_tool = "annotationDBI", ortholog_database = "HomoloGene", map_tissues, 
                   map_species, metric_type){
  if(length(gene_set) < 1){stop("Please provide a valid gene set(gene_set)")}
  if(missing(i_species)){stop("Please provide a valid input species (i_species)")}
  if(length(i_species) > 1){stop("You can provide only one input species (i_species)")}
  if(missing(input_id)){stop("Please provide a valid input gene identifier (input_id). Possible options are Ensembl_id, Entrez_id, Symbol. You can only provide one")}
  if(length(input_id) > 1){stop("You can provide only one input id (input_id)")}
  if(length(i_species) != length(input_id)) {stop("Input species (i_species) and input gene identifier (input_id) must be the same length")}
  if(missing(o_species)){stop("Please provide a valid output species (o_species)")}
  if(missing(output_ids)){stop("Please provide a valid output gene identifier (output_ids). Possible options are Ensembl_id, Enrez_id, Symbol. You can provide multiple ids")}
  if(length(mapping_tool) > 1){stop("You can provide only one database for mapping ids (mapping_tool)")}
  if(length(ortholog_database) > 1){stop("You can provide only one database for mapping orthologs (ortholog_database)")}
  if(missing(map_tissues)){stop("Please provide a valid tissue list (map_tissues)")}
  if(missing(map_species)){stop("Please provide a valid output species (map_species)")}
  if(missing(metric_type)){stop("Please provide a valid metric to be calculated (metric_type)")}
  methods::new("CoSIAn", gene_set = gene_set, i_species = i_species, input_id = input_id, o_species = o_species, output_ids = output_ids, 
      mapping_tool = mapping_tool, ortholog_database = ortholog_database, map_tissues = map_tissues, map_species = map_species, 
      metric_type = metric_type, converted_id = data.frame(0), gex = data.frame(0), metric = data.frame(0))
}
