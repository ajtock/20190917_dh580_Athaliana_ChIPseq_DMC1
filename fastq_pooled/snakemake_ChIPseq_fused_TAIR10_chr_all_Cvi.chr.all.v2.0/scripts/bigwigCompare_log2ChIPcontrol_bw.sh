#!/bin/bash

# Usage:
# csmit -m 50G -c 48 "bash ./bigwigCompare_log2ChIPcontrol_bw.sh WT_CENH3_Rep1_ChIP_SRR4430537 WT_REC8_Myc_Rep1_input 10000 10kb 48"

ChIPName=$1
controlName=$2
binSize=$3
binName=$4
threads=$5

[ -d ../mapped/both/log2ChIPcontrol ] || mkdir -p ../mapped/both/log2ChIPcontrol

source activate ChIPseq_mapping

bigwigCompare -b1 "../mapped/both/bw/"${ChIPName}"_MappedOn_t2t-col.20210610_lowXM_both_sort_norm.bw" \
              -b2 "/home/ajt200/analysis/REC8_pooled/snakemake_ChIPseq_t2t-col.20210610/mapped/both/bw/"${controlName}"_MappedOn_t2t-col.20210610_lowXM_both_sort_norm.bw" \
              -o "../mapped/both/log2ChIPcontrol/log2_"${ChIPName}"_"${controlName}"_MappedOn_t2t-col.20210610_lowXM_both_sort_norm_binSize"${binName}".bw" \
              -of "bigwig" \
              --pseudocount 1 \
              --operation log2 \
              --binSize ${binSize} \
              -p ${threads}

conda deactivate
