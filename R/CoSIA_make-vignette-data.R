## CoSIA Vignette Data 
## Amanda D. Clark
## Created 01/10/23
## Script is used to create data objects from the monogenic 
## disease-associated kidney gene set from Natera. 


# Load raw data
monogenic_kidney_genes <- readr::read_csv(file = "inst/extdata/raw/385_NateraKidney.csv", col_names = T)

# Output RDS
saveRDS(object = monogenic_kidney_genes, file = "inst/extdata/proccessed/monogenic_kidney_genes.rds")

# Output RDA
save(monogenic_kidney_genes, file = "inst/extdata/proccessed/monogenic_kidney_genes.rda")
