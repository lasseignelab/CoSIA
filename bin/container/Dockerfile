FROM bioconductor/bioconductor:devel-r-4.6.0
# includes: BIOCONDUCTOR_DOCKER_VERSION=3.23.22

##----------------------------------------------------------------------------##
## SYSTEM PACKAGES
##----------------------------------------------------------------------------##

RUN apt-get update --yes && \
  apt-get upgrade --yes

RUN apt-get install --yes \
  python3-pip \
  python3-dev \
  python3-pip \
  python3-tk


RUN apt-get install --yes \
  libcurl4-openssl-dev \
  libssl-dev \
  libcairo2-dev libxt-dev \
  libxml2-dev \
  libudunits2-dev \
  libhdf5-dev \
  libv8-dev \
  libgdal-dev \
  xorg libx11-dev libglu1-mesa-dev \
  samtools libboost-all-dev \
  libgsl-dev

RUN apt-get install -y \
	wget \
	unzip \
	git


##----------------------------------------------------------------------------##
##  R PACKAGES
##----------------------------------------------------------------------------##

#### Leave the top line the same for each new project. Uncomment lines to add ####
RUN R -e 'install.packages(c("tidyverse", "remotes",  "BiocManager", "devtools", "knitr", "lintr", "styler", "here", "markdown", "rmarkdown"))' \
#  && R -e 'install.packages(c("< insert CRAN repos here >"))' \
#  && R -e 'devtools::install_github("< insert GitHub repos here >")' \
#  && R -e 'remotes::install_github("< insert GitHub repos here >")' \
  && R -e 'BiocManager::install(c("BiocFileCache", "BiocGenerics", "AnnotationHub", "ExperimentHub", "CoSIA", "CoSIAdata"))'


##----------------------------------------------------------------------------##
## PYTHON PACKAGES
##----------------------------------------------------------------------------##


# RUN pip3 install <Insert Python packages here>


##----------------------------------------------------------------------------##
## 3RD PARTY PROGRAMS
##----------------------------------------------------------------------------##

#### Uncomment these lines if you need a conda env in your docker ####

# Install miniconda
# ENV CONDA_DIR /opt/conda
# RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
#     /bin/bash ~/miniconda.sh -b -p /opt/conda
# ENV PATH=$CONDA_DIR/bin:$PATH


##----------------------------------------------------------------------------##
## CONDA ENV
##----------------------------------------------------------------------------##

#### Uncomment these lines if you need a conda env in your docker ####

#RUN conda config --set channel_priority false && \
#    conda env create -f opt/BLAZE-1.1.0/conda_env/environment.yml

#SHELL ["conda", "run", "-n", "blaze", "/bin/bash", "-c"]

