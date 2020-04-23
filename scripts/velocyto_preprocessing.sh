#!/bin/bash

microwell=(ascending_colon, sigmoid_colon, transverse_colon)
chromium=(H508_EV H508_LSD1_KD HT29_EV HT29_LSD1_KD)
ncores=8

################################
## RNA-velocity Preprocessing ##
################################

module load singularity

## Prepare Sigularity Container
## ----------

if [ ! -f scrnaseq_software_velocyto_0.17.17.sif ]; then
  singularity pull --arch amd64 library://rpolicastro/default/scrnaseq_software:velocyto_0.17.17
fi

## 10X Chromium Samples
## ----------

for sample in ${chromium[@]}; do
  qsub -v SAMPLE=${sample},NCORES=${ncores} -l nodes=1:ppn=8,vmem=120gb,walltime=24:00:00 velocyto_10x.sh
done

## Microwell-seq Samples
## ----------

for sample in ${microwell[@]}; do
  qsub -v SAMPLE=${sample},NCORES=${ncores} -l nodes=1:ppn=8,vmem=120gb,walltime=24:00:00 velocyto_microwell.sh
done
