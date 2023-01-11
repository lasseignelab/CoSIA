## CoSIA Vignette Data 
## Amanda D. Clark
## Created 01/10/23
## Script is used to create data objects from the monogenic 
## disease-associated kidney gene set from Natera. 


# Set up
library(readr)

# Load raw data
monogenic_kidney_genes <- read_csv(file = "data/raw/385_NateraKidney.csv", col_names = T)

# Output RDS
saveRDS(object = monogenic_kidney_genes, file = "data/proccessed/monogenic_kidney_genes.rds")

# Output RDA
save(monogenic_kidney_genes, file = "data/proccessed/monogenic_kidney_genes.rda")
