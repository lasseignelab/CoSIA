#!/bin/bash

#SBATCH --ntasks=1  ## number of tasks (analyses) to run, i.e. # of nodes to request resource from
#SBATCH --cpus-per-task=12  ## the number of threads allocated to each task
#SBATCH --mem=64G  # memory per CPU core
#SBATCH --partition=express ## the partition to run in (short == 12h max run time)

cap_container -c singularity "lizzyr/bioc_cosia:1.10.0"
mkdir -p "${CAP_PROJECT_PATH}"/cache

singularity exec --cleanenv \
                 --containall \
                 -B "${CAP_PROJECT_PATH}" \
                 "${CAP_CONTAINER_PATH}/bioc_cosia_1.10.0.sif" \
                 Rscript --vanilla "${CAP_PROJECT_PATH}"/vignettes/vanilla_Intro.R "${CAP_PROJECT_PATH}"
