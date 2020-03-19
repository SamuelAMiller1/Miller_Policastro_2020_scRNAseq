#!/bin/bash

## Variables

workdir=/N/dc2/scratch/rpolicas/sam_scRNAseq
container=scrnaseq_software_cellranger_3.1.0.sif

###############################
## CellRanger Count of Reads ##
###############################

cd $workdir

## preparation Steps
## ----------

## Activate singularity on HPC.

module load singularity/3.5.2

## Download singularity container.

if [ ! -f ${workdir}/${container} ]; then
   singularity pull --arch amd64 library://rpolicastro/default/scrnaseq_software:cellranger_3.1.0
fi

## Download 10X hg38 genome.

if [ ! -d ${workdir}/genome ]; then
   mkdir -p ${workdir}/genome
fi

if [ ! -d ${workdir}/genome/refdata-cellranger-GRCh38-3.0.0 ]; then
   cd ${workdir}/genome
   wget http://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-GRCh38-3.0.0.tar.gz
   tar -xzvf refdata-cellranger-GRCh38-3.0.0.tar.gz
   cd ${workdir}
fi

## Run CellRanger
## ----------

## Count H508 EV

singularity exec -eCB $workdir -H $workdir $container \
cellranger count \
   --id H508_EV \
   --fastqs ${workdir}/sequences/H508/BL_Chrm_001_Ohagan_1 \
   --sample EV \
   --transcriptome ${workdir}/genome/refdata-cellranger-GRCh38-3.0.0 \
   --expect-cells 10000

## Count H508 LSD1 KD

singularity exec -eCB $workdir -H $workdir $container \
cellranger count \
   --id H508_LSD1_KD \
   --fastqs ${workdir}/sequences/H508/BL_Chrm_001_Ohagan_2 \
   --sample KO \
   --transcriptome ${workdir}/genome/refdata-cellranger-GRCh38-3.0.0 \
   --expect-cells 10000

## Count HT29 EV

singularity exec -eCB $workdir -H $workdir $container \
cellranger count \
   --id HT29_EV \
   --fastqs ${workdir}/sequences/HT29/Chrm_066_OHagan_BL_1 \
   --sample sm_ssrna1 \
   --transcriptome ${workdir}/genome/refdata-cellranger-GRCh38-3.0.0 \
   --expect-cells 10000

## Count HT29 LSD1 KD

singularity exec -eCB $workdir -H $workdir $container \
cellranger count \
   --id HT29_LSD1_KD \
   --fastqs ${workdir}/sequences/HT29/Chrm_066_OHagan_BL_2 \
   --sample sm_ssrna2 \
   --transcriptome ${workdir}/genome/refdata-cellranger-GRCh38-3.0.0 \
   --expect-cells 10000

## Count Colon 1

singularity exec -eCB $workdir -H $workdir $container \
cellranger count \
   --id COLON_1 \
   --fastqs ${workdir}/sequences/colon_1 \
   --sample colon1 \
   --transcriptome ${workdir}/genome/refdata-cellranger-GRCh38-3.0.0 \
   --expect-cells 10000
