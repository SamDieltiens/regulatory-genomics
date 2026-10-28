#!/usr/bin/env bash
# Same STAR command as in part 1 of the notebook; waits until each sample is downloaded.
set -uo pipefail
cd "$(dirname "$0")/../data"
GENOME_DIR=genome
LOG=logs/download.log   # written by download_fastq.sh
mkdir -p alignment
for SAMPLE_NAME in $(cut -f 2 fastq_files/samples.tsv)
do
    until grep -qx "done ${SAMPLE_NAME}" $LOG; do sleep 60; done
    [ -s alignment/${SAMPLE_NAME}.Log.final.out ] && continue
    echo -e "Aligning ${SAMPLE_NAME} $(date)"
    STAR \
        --genomeDir ${GENOME_DIR}/GRCh38.gencode_v46.star_d \
        --runThreadN 16 \
        --readFilesIn fastq_files/${SAMPLE_NAME}_1.fastq.gz fastq_files/${SAMPLE_NAME}_2.fastq.gz \
        --readFilesCommand zcat \
        --genomeLoad LoadAndKeep \
        --limitBAMsortRAM 30000000000 \
        --outFileNamePrefix alignment/${SAMPLE_NAME}. \
        --outSAMtype BAM SortedByCoordinate
done
STAR --genomeDir ${GENOME_DIR}/GRCh38.gencode_v46.star_d --genomeLoad Remove --outFileNamePrefix alignment/remove_genome.
