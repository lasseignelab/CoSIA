# CoSIA: **C**r**o**ss **S**pecies **I**nvestigation and **A**nalysis

Cross-Species Investigation and Analysis (CoSIA) is a package that provides 
  researchers with an alternative methodology for comparing across species and 
  tissues using normal wild-type RNA-Seq Gene Expression data from Bgee. Using 
  RNA-Seq Gene Expression data, CoSIA provides multiple visualization tools to 
  explore the transcriptome diversity and variation across genes, tissues, and 
  species. CoSIA uses the Coefficient of Variation and Shannon Entropy and 
  Specificity to calculate transcriptome diversity and variation. CoSIA also 
  provides additional conversion tools and utilities to provide a streamlined 
  methodology for cross-species comparison.

CoSIA is an R package that provides researchers with the tools to measure and visualize gene-expression metrics in order to compare across five commonly used biomedical research species (Mus musculus, Rattus norvegicus, Danio rerio, Drosophila melanogaster, and Caenorhabditis elegans) in addition to Homo sapiens and their tissues.

Specifically, CoSIA uses curated non-diseased wild-type RNA-sequencing expression data, from Bgee, to visualize a gene's expression across tissues and species. CoSIA also streamlines conversions between gene identifiers among the same species and different species.

<img src="inst/images/CoSIA_Workflow.png" alt="Figure 1. CoSIA_Workflow" width="703"/>

CoSIA is split into 3 modules that provide various resources in order for researchers to conduct cross species analysis using gene expression metrics.

The first module and second module have a shared method, `getConversion`, that is used to conduct conversion between different gene identifiers in the same species as well as in different species. The third module has a method, `getTissueExpression`, that is used to visualize the raw gene expression values of a gene across multiple tissues in one model organism. The fourth module has a method, `getSpeciesExpression`, that is used to visualize the raw gene expression values of a gene in one tissues across multiple model organisms.The fifth module has a method, `getTranscriptomeDiversity`, that is used to visualize median-based Coefficient of Variation and Shannon Entropy to look at the variation and diversity of gene expression across tissues and model organisms.

## Getting Started

### Installing

In R:

``` r
library(devtools)
install_github("lasseignelab/CoSIA", ref= "main", auth_token = "<PAT>")
```

## How to use the CoSIA package

### Load the package

``` r
library(CoSIA)
```

## Authors

-   Anisha Haldar
-   Vishal H. Oza
-   Nathaniel DeVoss
-   Amanda Clark
-   Brittany N. Lasseigne

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
