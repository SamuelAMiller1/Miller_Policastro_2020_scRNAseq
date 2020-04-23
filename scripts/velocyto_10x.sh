#!/bin/bash

##########################################
## RNA-velocity on 10X Chromium Samples ##
##########################################

cd $PBS_O_WORKDIR
cd ..

module load singularity

singularity exec -eCB `pwd` -H `pwd` scrnaseq_software_velocyto_0.17.17.sif \
velocyto run10x \
  -@ $NCORES \
  --samtools-memory 5000 \
  -m genome/repeat_mask.gtf \
  aligned/$SAMPLE \
  genome/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf
