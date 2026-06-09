# validation/gene_expression.R
#
# Creates the gene expression plot that was not working after updating CoSIA.
#
# Usage (from repo root):
#   Run the script in RStudio.
#
# Output: Gene expression plot in the RStudio Viewer window

devtools::load_all(".")

loaded_version <- as.character(packageVersion("CoSIA"))
loaded_path    <- find.package("CoSIA")

cat("CoSIA version  :", loaded_version, "\n")
cat("CoSIA path     :", loaded_path, "\n\n")

GENE_SET   <- c("ENSG00000008710", "ENSG00000118762")
I_SPECIES  <- "h_sapiens"
INPUT_ID   <- "Ensembl_id"
OUTPUT_IDS <- c("Ensembl_id")
TOOL       <- "annotationDBI"
ORTHO_DB   <- "HomoloGene"
TISSUE     <- "brain"
MULTI_SPECIES <- c("h_sapiens", "m_musculus")

make_obj <- function(map_species, map_tissues, metric_type,
  o_species = map_species) {
  obj <- CoSIAn(
    gene_set          = GENE_SET,
    i_species         = I_SPECIES,
    input_id          = INPUT_ID,
    o_species         = o_species,
    output_ids        = OUTPUT_IDS,
    mapping_tool      = TOOL,
    ortholog_database = ORTHO_DB,
    map_tissues       = map_tissues,
    map_species       = map_species,
    metric_type       = metric_type
  )
  getConversion(obj)
}

run_metric <- function(label, ...) {
  cat("---", label, "---\n")
  tryCatch({
    obj    <- make_obj(...)
    obj    <- getGEx(obj)
    cat("  SUCCESS \n")
    obj
  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n\n")
    NULL
  })
}
  
gex_species_multi = run_metric(
  "DS_Gene multi-species",
  map_species  = MULTI_SPECIES,
  map_tissues  = TISSUE,
  metric_type  = "DS_Gene",
  o_species    = MULTI_SPECIES
)

print(CoSIA::plotSpeciesGEx(gex_species_multi, "brain", "ENSG00000008710"))