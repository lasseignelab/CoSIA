# CoSIA: **C**r**o**ss **S**pecies **I**nvestigation and **A**nalysis

CoSIA is an R package that provides researchers with the tools to measure and visualize gene-expression metrics in order to compare across model organisms and their tissues. Specifically, CoSIA uses curated non-diseased wild-type RNA-sequencing expression data, from Bgee, to visualize a gene's expression across tissues and model organisms. CoSIA also streamlines conversions between gene identifiers among the same species and different species.

<img src="images/CoSIA_Workflow.png" alt="Figure 1. CoSIA_Workflow" width="705"/>

CoSIA is split into 5 modules that provide various resources in order for researchers to conduct cross species analysis using gene expression metrics.

The first module and second module have a shared method, `getConversion`, that is used to conduct conversion between different gene identifiers in the same species as well as in different species. The third module has a method, `getTissueExpression`, that is used to visualize the raw gene expression values of a gene across multiple tissues in one model organism. The fourth module has a method, `getSpeciesExpression`, that is used to visualize the raw gene expression values of a gene in one tissues across multiple model organisms.The fifth module has a method, `getTranscriptomeDiversity`, that is used to visualize median-based Coefficient of Variation and Shannon Entropy to look at the variation and diversity of gene expression across tissues and model organisms.

## Getting Started

### Installing

In R:

``` r
install.packages("devtools")
devtools::install_github("lasseignelab/CoSIA", ref= "main", auth_token = "")
#add personal authorization token for the auth_token argument in order to properly load the CoSIA package
```

## How to use the CoSIA package

### Load the package

``` r
library(CoSIA)
```

### Test Case: 10 Monogenic Kidney Disease Genes

In the following example, we will aim to portray how CoSIA can be used using 10 genes from [Natera's Monogenic Kidney Gene Panel](https://www.natera.com/resource-library/renasight/385-genes-associated-with-monogenic-disorders-linked-to-kidney-disease). First we will create a data frame with 10 Human Gene Symbols and Conditions from this panel.

``` r
Monogenic_Kidney_Genes <- data.frame (
  Human_Gene_Symbols = 
)
```

## Authors

-   Anisha Haldar
-   Vishal H. Oza
-   Nathaniel DeVoss
-   Amanda Clark
-   Brittany N. Lasseigne

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
