#!/bin/bash

##################
## Marker Genes ##
##################

cd $PBS_O_WORKDIR
cd ..

module load singularity/3.5.2

singularity exec -eCB `pwd` -H `pwd` scrnaseq_software_seurat_velocytor_0.3.sif \
Rscript scripts/seurat_marker_genes.R
