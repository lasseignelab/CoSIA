## CoSIA Vignette Data
## Author: Amanda D. Clark
## Date Created: 01/10/23
## Purpose: Script is used to create data objects from the monogenic
## disease-associated kidney gene set from Natera.


## Natera's Monogenic Kidney Disease Panel (https://www.natera.com/resource-library/renasight/385-genes-associated-with-monogenic-disorders-linked-to-kidney-disease)
## can be downloaded using the link above. The panel is only available as a pdf
## Genes from the pdf were manually typed into the csv 
## file available in "/inst/extdata/raw"

# Load raw data
monogenic_kidney_genes <-
    readr::read_csv(
        file = "inst/extdata/raw/385_NateraKidney.csv",
        col_names = TRUE
    )

# Output RDS
saveRDS(
    object = monogenic_kidney_genes,
    file = "inst/extdata/proccessed/monogenic_kidney_genes.rds"
)

# Output RDA
save(monogenic_kidney_genes,
    file = "inst/extdata/proccessed/monogenic_kidney_genes.rda"
)
