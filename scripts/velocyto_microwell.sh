#!/bin/bash

###########################################
## RNA-velocity on Microwell-seq Samples ##
###########################################

cd $PBS_O_WORKDIR
cd ..

module load singularity

if [ ! -d aligned/${SAMPLE}_velocyto ]; then
  mkdir -p aligned/${SAMPLE}_velocyto
fi

singularity exec -eCB `pwd` -H `pwd` scrnaseq_software_velocyto_0.17.17.sif \
velocyto run \
  -o aligned/${SAMPLE}_velocyto \
  -@ $NCORES \
  --samtools-memory 5000 \
  -m genome/repeat_mask.gtf \
  aligned/$SAMPLE/aligned_${SAMPLE}_possorted.bam \
  genome/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf
