# CoSIA: **C**r**o**ss **S**pecies **I**nvestigation and **A**nalysis

CoSIA is an R package that provides researchers with the tools to measure and visualize gene-expression metrics in order to compare across model organisms and their tissues. Specifically, CoSIA uses curated non-diseased wild-type RNA-sequencing expression data, from Bgee, to visualize a gene’s expression across tissues and model organisms. CoSIA also streamlines conversions between gene identifiers among the same species and different species.

<img src="images/CoSIA Workflow.png" alt="Figure 1. CoSIA Workflow" width="705"/>

CoSIA is split into 5 modules that provide various resources in order for researchers to conduct cross species analysis using gene expression metrics.

The first module and second module have a shared method, `getConversion`, that is used to conduction conversion between different gene identifiers in the same species as well as in different species. The third module has a method, `getTissueExpression`, that is used to visualize the raw gene expression values of a gene across multiple tissues in one model organism. The fourth module has a method, `getSpeciesExpression`, that is used to visualize the raw gene expression values of a gene in one tissues across multiple model organisms.The fifth module has a method, `getTranscriptomeDiversity`, that is used to visualize median-based Coefficient of Variation and Shannon Entropy to look at the variation and diversity of gene expression across tissues and model organisms.

## Getting Started

### Installing

In R:

``` R
install_github("lasseignelab/CoSIA", ref= "main", auth_token = "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")
```
### How to use CoSIA


## Running the tests

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

    Give an example

### And coding style tests

Explain what these tests test and why

    Give an example

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

-   [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
-   [Maven](https://maven.apache.org/) - Dependency Management
-   [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags).

## Authors

-   **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

-   Hat tip to anyone whose code was used
-   Inspiration
-   etc
