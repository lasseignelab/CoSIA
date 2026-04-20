library(CoSIA)
wd <- commandArgs(trailingOnly = TRUE)
print(paste0("Project directory: ", wd))

cat("\nSetting cache...\n")
cachedir <- paste0(wd, "/cache")
setExperimentHubOption("ASK", FALSE)
setExperimentHubOption("CACHE", cachedir)

# Force a clean load of the hub
Sys.setenv(EXPERIMENT_HUB_CACHE = cachedir)
eh <- ExperimentHub(localHub = FALSE)
cat("\nLoading in...\n")
print(paste0(wd, "/inst/extdata/proccessed/monogenic_kidney_genes.rda"))
load(paste0(wd, "/inst/extdata/proccessed/monogenic_kidney_genes.rda"))


# downsampling data for figure
set.seed(42)
monogenic_kid_sample <- dplyr::sample_n(monogenic_kidney_genes, 20)

CoSIA::getTissues("d_rerio")

CoSIA::getTissues(c("m_musculus", "r_norvegicus"))

### Initializing a CoSIAn object

CoSIAn_Obj <- CoSIA::CoSIAn(
    gene_set = unique(monogenic_kid_sample$Gene),
    i_species = "h_sapiens",
    o_species = c(
        "h_sapiens",
        "r_norvegicus"
    ),
    input_id = "Symbol",
    output_ids = "Ensembl_id",
    map_species = c(
        "h_sapiens",
        "r_norvegicus"
    ),
    map_tissues = c(
        "adult mammalian kidney",
        "heart"
    ),
    mapping_tool = "annotationDBI",
    ortholog_database = "HomoloGene",
    metric_type = "CV_Species"
)

str(CoSIAn_Obj)

## Use Cases with Monogenic Kidney Disease-Associated Genes
### Use Case #1: Converting Gene Symbols to Ensembl IDs (`getConversion`)
CoSIAn_Obj_convert <- CoSIA::getConversion(CoSIAn_Obj)

head(CoSIA::viewCoSIAn(CoSIAn_Obj_convert, "converted_id"))

### Use Case #2: Obtaining and visualizing curated non-diseased kidney and heart gene expression data for human, mouse, rat from Bgee
CoSIAn_Obj_gex <- CoSIA::getGEx(CoSIAn_Obj_convert)

head(CoSIA::viewCoSIAn(CoSIAn_Obj_gex, "gex"), 5)

CoSIAn_Obj_gexplot <- CoSIA::plotSpeciesGEx(CoSIAn_Obj_gex, "adult mammalian kidney", "ENSG00000136463")


### Use Case #3: Gene expression variability across species for kidney tissue by calculating and visualizing median-based Coefficient of Variation (CV)

CoSIAn_Obj_CV <- CoSIA::getGExMetrics(CoSIAn_Obj_convert)
saveRDS(CoSIAn_Obj_CV, paste0(wd, "/vignettes/CoSIAn_Obj_CV.rds"))

CoSIAn_Obj_CVplot <- CoSIA::plotCVGEx(CoSIAn_Obj_CV)

### Use Case #4: Gene expression diversity and specificity across tissues and species for monogenic kidney-disease associated genes
CoSIAn_Obj_DS <- CoSIA::CoSIAn(
    gene_set = unique(monogenic_kid_sample$Gene),
    i_species = "h_sapiens",
    o_species = c("h_sapiens", "r_norvegicus"),
    input_id = "Symbol",
    output_ids = "Ensembl_id",
    map_species = c("h_sapiens", "r_norvegicus"),
    map_tissues = c("adult mammalian kidney", "heart"),
    mapping_tool = "annotationDBI",
    ortholog_database = "HomoloGene",
    metric_type = "DS_Tissue"
)

CoSIAn_Obj_DS <- CoSIA::getConversion(CoSIAn_Obj_DS)

CoSIAn_Obj_DS <- CoSIA::getGExMetrics(CoSIAn_Obj_DS)
saveRDS(CoSIAn_Obj_DS, paste0(wd, "/vignettes/CoSIAn_Obj_DS.rds"))

CoSIAn_Obj_DSplot <- CoSIA::plotDSGEx(CoSIAn_Obj_DS)

### Save plots
plots <- c(CoSIAn_Obj_gexplot, CoSIAn_Obj_CVplot, CoSIAn_Obj_DSplot)

pdf(paste0(wd, "/vignettes/plots.pdf"))
for plot in seq_along(plots){
    print(plot)
}
dev.off()

sessionInfo()

fptm <- proc.time()
cat("\nComplete! Total runtime: ", (fptm[3] / 60), " minutes")

