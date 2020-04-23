#!/bin/bash

## Variables

container=scrnaseq_software_drop_seq_2.3.0.sif
ncores=8

########################################
## Microwell-seq Human Cell Landscape ##
########################################

cd $PBS_O_WORKDIR
cd ..

module load singularity

## Preparation Steps
## ----------

## Download singularity container.

if [ ! -f ${container} ]; then
   singularity pull --arch amd64 library://rpolicastro/default/scrnaseq_software:drop_seq_2.3.0
fi

## Download 10X hg38 genome.

if [ ! -d genome ]; then
   mkdir -p genome
fi

if [ ! -d genome/refdata-cellranger-GRCh38-3.0.0 ]; then
   cd genome
   wget http://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-GRCh38-3.0.0.tar.gz
   tar -xzvf refdata-cellranger-GRCh38-3.0.0.tar.gz
   cd ..
fi

## Analysis of Microwell-seq Data
## ----------

## Create sequence dictionary from genome assembly.

if [ ! -d microwell_meta_data ]; then
 mkdir -p microwell_meta_data
fi

singularity exec -eCB "$(pwd)" -H "$(pwd)" scrnaseq_software_drop_seq_2.3.0.sif \
picard CreateSequenceDictionary \
  R=genome/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa \
  O=microwell_meta_data/genome.dict \
  SPECIES=Homo_sapiens

## Convert genome annotation to RefFlat.

singularity exec -eCB "$(pwd)" -H "$(pwd)" scrnaseq_software_drop_seq_2.3.0.sif \
ConvertToRefFlat \
  ANNOTATIONS_FILE=genome/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf \
  SEQUENCE_DICTIONARY=microwell_meta_data/genome.dict \
  OUTPUT=microwell_meta_data/genome.refFlat

## Convert fastq files to unaligned BAM.

samples=(ascending_colon transverse_colon sigmoid_colon)

for sample in ${samples[@]}; do
  if [ ! -d aligned/$sample ]; then
    mkdir -p aligned/$sample
  fi
done

mkdir -p tempdir

# Ascending colon.
singularity exec -eCB "$(pwd)" -H "$(pwd)" scrnaseq_software_drop_seq_2.3.0.sif \
picard FastqToSam \
  F1=sequences/ascending_colon/SRR9843410_1.fastq \
  F2=sequences/ascending_colon/SRR9843410_2.fastq \
  O=aligned/ascending_colon/unaligned_ascending_colon.bam \
  SAMPLE_NAME=ascending_colon \
  TMP_DIR=tempdir

# Sigmoid colon.
singularity exec -eCB "$(pwd)" -H "$(pwd)" scrnaseq_software_drop_seq_2.3.0.sif \
picard FastqToSam \
  F1=sequences/sigmoid_colon/SRR9887686_1.fastq \
  F2=sequences/sigmoid_colon/SRR9887686_2.fastq \
  O=aligned/sigmoid_colon/unaligned_sigmoid_colon.bam \
  SAMPLE_NAME=sigmoid_colon \
  TMP_DIR=tempdir

# Transverse colon.
singularity exec -eCB "$(pwd)" -H "$(pwd)" scrnaseq_software_drop_seq_2.3.0.sif \
picard FastqToSam \
  F1=sequences/transverse_colon/SRR9887678_1.fastq \
  F2=sequences/transverse_colon/SRR9887678_2.fastq \
  O=aligned/transverse_colon/unaligned_transverse_colon.bam \
  SAMPLE_NAME=transverse_colon \
  TMP_DIR=tempdir

## Extract cell barcode.

for sample in ${samples[@]}; do
  singularity exec -eCB "$(pwd)" -H "$(pwd)" scrnaseq_software_drop_seq_2.3.0.sif \
  TagBamWithReadSequenceExtended \
    INPUT=aligned/${sample}/unaligned_${sample}.bam \
    OUTPUT=aligned/${sample}/unaligned_${sample}_cellbc.bam \
    SUMMARY=aligned/${sample}/unaligned_${sample}_cellbc.bam_summary.txt \
    BASE_RANGE=1-6:22-27:43-48 \
    BASE_QUALITY=10 \
    BARCODED_READ=1 \
    DISCARD_READ=False \
    TAG_NAME=XC \
    NUM_BASES_BELOW_QUALITY=1
done

## Extract UMI.

for sample in ${samples[@]}; do
  singularity exec -eCB "$(pwd)" -H "$(pwd)" scrnaseq_software_drop_seq_2.3.0.sif \
  TagBamWithReadSequenceExtended \
    INPUT=aligned/${sample}/unaligned_${sample}_cellbc.bam \
    OUTPUT=aligned/${sample}/unaligned_${sample}_cellbc_umi.bam \
    SUMMARY=aligned/${sample}/unaligned_${sample}_cellbc_umi.bam_summary.txt \
    BASE_RANGE=49-54 \
    BASE_QUALITY=10 \
    BARCODED_READ=1 \
    DISCARD_READ=True \
    TAG_NAME=XM \
    NUM_BASES_BELOW_QUALITY=1
done

## Remove reads with low cell or UMI barcode quality.

for sample in ${samples[@]}; do
  singularity exec -eCB "$(pwd)" -H "$(pwd)" scrnaseq_software_drop_seq_2.3.0.sif \
  FilterBam \
    TAG_REJECT=XQ \
    INPUT=aligned/${sample}/unaligned_${sample}_cellbc_umi.bam \
    OUTPUT=aligned/${sample}/unaligned_${sample}_filtered.bam
done

## Trim barcodes from reads.

for sample in ${samples[@]}; do
  singularity exec -eCB "$(pwd)" -H "$(pwd)" scrnaseq_software_drop_seq_2.3.0.sif \
  TrimStartingSequence \
    INPUT=aligned/${sample}/unaligned_${sample}_filtered.bam \
    OUTPUT=aligned/${sample}/unaligned_${sample}_adapter_trimmed.bam \
    OUTPUT_SUMMARY=aligned/${sample}/unaligned_${sample}_adapter_trimmed.bam_summary.txt \
    SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG \
    MISMATCHES=0 \
    NUM_BASES=5
done

## Trim polyA tails from reads.

for sample in ${samples[@]}; do
  singularity exec -eCB "$(pwd)" -H "$(pwd)" scrnaseq_software_drop_seq_2.3.0.sif \
  PolyATrimmer \
    INPUT=aligned/${sample}/unaligned_${sample}_adapter_trimmed.bam \
    OUTPUT=aligned/${sample}/unaligned_${sample}_adapter_polya_trimmed.bam \
    OUTPUT_SUMMARY=aligned/${sample}/unaligned_${sample}_adapter_polya_trimmed.bam_summary.txt \
    MISMATCHES=0 \
    NUM_BASES=6 \
    USE_NEW_TRIMMER=true
done

## Convert SAM files back to FASTQ for alignment.

for sample in ${samples[@]}; do
  singularity exec -eCB "$(pwd)" -H "$(pwd)" scrnaseq_software_drop_seq_2.3.0.sif \
  picard SamToFastq \
    INPUT=aligned/${sample}/unaligned_${sample}_adapter_polya_trimmed.bam \
    FASTQ=aligned/${sample}/unaligned_${sample}_filtered.fastq
done

## Create STAR genome index.

mkdir -p microwell_meta_data/star_index

singularity exec -eCB "$(pwd)" -H "$(pwd)" scrnaseq_software_drop_seq_2.3.0.sif \
STAR \
  --runThreadN $ncores \
  --runMode genomeGenerate \
  --genomeDir microwell_meta_data/star_index \
  --genomeFastaFiles genome/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa \
  --sjdbGTFfile genome/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf

## Align reads with STAR.

for sample in ${samples[@]}; do
  singularity exec -eCB "$(pwd)" -H "$(pwd)" scrnaseq_software_drop_seq_2.3.0.sif \
  STAR \
    --runThreadN $ncores \
    --genomeDir microwell_meta_data/star_index \
    --readFilesIn aligned/${sample}/unaligned_${sample}_filtered.fastq \
    --outFileNamePrefix aligned/${sample}/${sample}
done

## Sort SAM files by queryname and turn to BAM.

for sample in ${samples[@]}; do
  singularity exec -eCB "$(pwd)" -H "$(pwd)" scrnaseq_software_drop_seq_2.3.0.sif \
  picard SortSam \
    I=aligned/${sample}/${sample}Aligned.out.sam \
    O=aligned/${sample}/aligned_${sample}_sorted.bam \
    SO=queryname \
    TMP_DIR=tempdir
done

## Merge back the barcode information.

for sample in ${samples[@]}; do
  singularity exec -eCB "$(pwd)" -H "$(pwd)" scrnaseq_software_drop_seq_2.3.0.sif \
  picard MergeBamAlignment \
    REFERENCE_SEQUENCE=genome/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa \
    UNMAPPED_BAM=aligned/${sample}/unaligned_${sample}_adapter_polya_trimmed.bam \
    ALIGNED_BAM=aligned/${sample}/aligned_${sample}_sorted.bam \
    OUTPUT=aligned/${sample}/aligned_${sample}_merged.bam \
    INCLUDE_SECONDARY_ALIGNMENTS=false \
    PAIRED_RUN=false \
    TMP_DIR=tempdir
done

## Detect bead substitution errors.

#for sample in ${samples[@]}; do
#  singularity exec -eCB "$(pwd)" -H "$(pwd)" scrnaseq_software_drop_seq_2.3.0.sif \
#  DetectBeadSubstitutionErrors \
#    I=aligned/${sample}/aligned_${sample}_merged.bam \
#    O=aligned/${sample}/aligned_${sample}_bead_substitution.bam \
#    OUTPUT_REPORT=aligned/${sample}/aligned_${sample}_bead_substitution.bam_summary.txt \
#    NUM_THREADS=$ncores
#    CELL_BARCODE_TAG="XC:Z" \
#    MOLECULAR_BARCODE_TAG="XM:Z"
#done

## Detect bead synthesis errors.

#for sample in ${samples[@]}; do
#  singularity exec -eCB "$(pwd)" -H "$(pwd)" scrnaseq_software_drop_seq_2.3.0.sif \
#  DetectBeadSynthesisErrors \
#    I=aligned/${sample}/aligned_${sample}_bead_substitution.bam \
#    O=aligned/${sample}/aligned_${sample}_clean.bam \
#    REPORT=aligned/${sample}/aligned_${sample}_clean.bam_report.txt \
#    OUTPUT_STATS=aligned/${sample}/aligned_${sample}_clean.bam_stats.txt \
#    SUMMARY=aligned/${sample}/aligned_${sample}_clean.bam_summary.txt \
#    PRIMER_SEQUENCE=AAGCAGTGGTATCAACGCAGAGTAC \
#    NUM_THREADS=$ncores
#done

## Digital gene expression.

for sample in ${samples[@]}; do
  singularity exec -eCB "$(pwd)" -H "$(pwd)" scrnaseq_software_drop_seq_2.3.0.sif \
  DigitalExpression \
    I=aligned/${sample}/aligned_${sample}_merged.bam \
    O=aligned/${sample}/${sample}_count_matrix.tsv \
    CELL_BARCODE_TAG="XC" \
    MOLECULAR_BARCODE_TAG="XM" \
    NUM_CORE_BARCODES=5000 \
    SUMMARY=aligned/${sample}/${sample}_count_matrix.tsv_summary.txt
done
