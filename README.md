[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7674870.svg)](https://doi.org/10.5281/zenodo.7674870)

# CoSIA: **C**r**o**ss **S**pecies **I**nvestigation and **A**nalysis

**C**r**o**ss **S**pecies **I**nvestigation and **A**nalysis (`CoSIA`) is a package that provides researchers with an alternative methodology for comparing across species and tissues using normal wild-type RNA-Seq Gene Expression data from Bgee. Using RNA-Seq Gene Expression data, CoSIA provides multiple visualization tools to explore the transcriptome diversity and variation across genes, tissues, and species. CoSIA uses Coefficient of Variation and Shannon Entropy and Specificity to calculate transcriptome diversity and variation. CoSIA also provides additional conversion tools and utilities to provide a streamlined methodology for cross-species comparison across the tissues and genes of five commonly used biomedical research species (*Mus musculus*, *Rattus norvegicus*, *Danio rerio*, *Drosophila melanogaster*, and *Caenorhabditis elegans*) in addition to *Homo sapiens*.

<img src="inst/images/CoSIA_Workflow.png" alt="Figure 1. CoSIA_Workflow" width="703"/>

CoSIA is split into 3 methods that provide various resources in order for researchers to conduct cross species analysis using gene expression metrics.

The first method `getConversion` is used to conduct gene identifier and ortholog conversions.

The second method `getGEx` is used to obtain the raw gene expression values (Variance Stabilizing Transformation of Read Counts) of a gene or set of gene among tissues or species of choice. The data obtained in `getGEx` can then be visualized for a single gene *across multiple tissues in single model organism* through `plotTissueGEx` or *across multiple species in a single tissue* through `plotSpeciesGEx`.

The third method `getGExMetrics` provides various methodologies in calculating median-based Coefficient of Variation (variability) and Shannon Entropy (diversity & specificity). Once calculated, `plotCVGEx` & `plotDSGEx` can be used to visualize the variation and diversity & specificity (DS) of gene expression across genes, tissues, and species.

## Getting Started

### Installing

**In preparation for the Bioconductor 3.17 release, we have developed CoSIA within a Bioconductor docker. Follow the instructions on [Bioconductor's Docker Help Page](https://www.bioconductor.org/help/docker/) to install and run CoSIA within the `bioconductor_docker:devel` container.**

In R:

``` r
library(devtools)
install_github("lasseignelab/CoSIA", ref= "main", auth_token = "<PAT>")
library(CoSIA)
```

### Using CoSIA

(add link to Vignette for users to be able to get more information on 
running CoSIA)

## Authors

-   Anisha Haldar
-   Vishal H. Oza
-   Nathaniel DeVoss
-   Amanda D. Clark
-   Brittany N. Lasseigne

## Lasseigne Lab

[What is Happening in the Lasseigne Lab?](https://www.lasseigne.org/)

<img src="https://www.lasseigne.org/img/main/lablogo.png" width="75" height="75">

## Funding

This work was supported in part by the UAB Lasseigne Lab Start-Up funds 
(BNL, AH, ND, ADC and VHO), the UAB Pilot Center for Precision Animal Modeling 
(C-PAM) (1U54OD030167) (BNL and VHO), UAB Pilot Center for Precision Animal 
Modeling (C-PAM) - Diversity Supplement (3U54OD030167-03S1) (ADC), and 
Mentored Experiences in Research, Instruction, and Teaching (MERIT) 
Program (K12 GM088010) (ADC).

## Acknowledgements

The authors would like to thank the members of the Lassigne lab for their 
support and feedback, in particular, Elizabeth J. Wilk, Jordan Whitlock, 
and Timothy C. Howton.

## License

This project is licensed under the MIT License - 
see the [LICENSE.md](LICENSE.md) file for details
