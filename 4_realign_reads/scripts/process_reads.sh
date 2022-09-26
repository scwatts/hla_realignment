#!/usr/bin/env bash
set -euo pipefail

bam_gds="${1}"
reference_fp="${2}"

bam_gds_fn=${bam_gds##*/}
sample=${bam_gds_fn/.bam/}

# Get HLA reads; we ignore read mates not in specified region in this approach
mkdir -p 1_sliced/${sample}/
bam_fp_sliced=1_sliced/${sample}/${sample}.bam

bam_presign=$(ica files get -o json ${bam_gds} | jq -r '.presignedUrl')
bai_presign=$(ica files get -o json ${bam_gds}.bai | jq -r '.presignedUrl')

samtools view \
  -M \
  -L ../3_compile_regions/regions.bed \
  -o ${bam_fp_sliced} \
  "${bam_presign}##idx##${bai_presign}"

# Convert HLA reads to FASTQ
sambamba sort -n ${bam_fp_sliced} -o ${bam_fp_sliced/.bam/_sorted.bam}
samtools fastq -@2 ${bam_fp_sliced/.bam/_sorted.bam} \
    -1 1_sliced/${sample}/${sample}_R1.fastq.gz \
    -2 1_sliced/${sample}/${sample}_R2.fastq.gz \
    -0 1_sliced/${sample}/${sample}_other.fastq.gz \
    -s 1_sliced/${sample}/${sample}_singleton.fastq.gz;

# Align
mkdir -p 2_alignments/${sample}/
bwa mem \
  -t2 \
  -Y \
  ${reference_fp} \
  1_sliced/${sample}/${sample}_R1.fastq.gz \
  1_sliced/${sample}/${sample}_R2.fastq.gz | \
  samtools sort -T tmp -o 2_alignments/${sample}/${sample}_realigned.pe.bam

bwa mem \
  -t2 \
  -Y \
  ${reference_fp} \
  1_sliced/${sample}/${sample}_singleton.fastq.gz | \
  samtools sort -T tmp -o 2_alignments/${sample}/${sample}_realigned.se.bam

# Combine alignments
sambamba merge \
  2_alignments/${sample}/${sample}_realigned.bam \
  2_alignments/${sample}/${sample}_realigned.pe.bam \
  2_alignments/${sample}/${sample}_realigned.se.bam
