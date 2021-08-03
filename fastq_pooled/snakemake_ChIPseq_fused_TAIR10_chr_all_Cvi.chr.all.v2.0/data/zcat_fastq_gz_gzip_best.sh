#!/bin/bash

# Usage:
# csmit -m 20G -c 10 "bash ./zcat_fastq_gz_gzip_best.sh ColCviF1_DMC1_V5_Rep1 ColCviF1_DMC1_V5_Rep1_ChIP_R1 ColCviF1_DMC1_V5_Rep1_ChIP_R2"
# csmit -m 20G -c 10 "bash ./zcat_fastq_gz_gzip_best.sh ColColF1_DMC1_V5_Rep1 ColColF1_DMC1_V5_Rep1_ChIP_R1 ColColF1_DMC1_V5_Rep1_ChIP_R2"

inName=$1
outNameR1=$2
outNameR2=$3

if [ ! -f "$outNameR1.fastq.gz" ]; then
  zcat $inName"_Run1_ChIP_R1.fastq.gz" \
       $inName"_Run2_ChIP_R1.fastq.gz" \
       | gzip -c -k --best > $outNameR1.fastq.gz;
else
  echo "skipping $outNameR1"
fi

if [ ! -f "$outNameR2.fastq.gz" ]; then
  zcat $inName"_Run1_ChIP_R2.fastq.gz" \
       $inName"_Run2_ChIP_R2.fastq.gz" \
       | gzip -c -k --best > $outNameR2.fastq.gz;
else
  echo "skipping $outNameR2"
fi
