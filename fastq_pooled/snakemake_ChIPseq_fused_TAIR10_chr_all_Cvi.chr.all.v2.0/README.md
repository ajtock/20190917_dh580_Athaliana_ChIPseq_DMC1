# Automated workflow for paired-end ChIP-seq and MNase-seq data processing and alignment

This is a Snakemake workflow for automated processing and alignment of paired-end read data derived from chromatin immunoprecipitation or micrococcal nuclease digestion followed by high-throughput sequencing (ChIP-seq or MNase-seq).

### Requirements

- Installation of [Snakemake](https://snakemake.readthedocs.io/en/stable/) and optionally [conda](https://conda.io/docs/)
- Demultiplexed paired-end reads in gzipped FASTQ format located in the `data/` directory. These should be named according to the following naming convention: `{sample}_R1.fastq.gz` and `{sample}_R2.fastq.gz`
- A samtools-indexed reference genome in FASTA format and a chromosome sizes file (e.g., `wheat_v1.0.fa`, `wheat_v1.0.fa.fai`, and `wheat_v1.0.fa.sizes`, the latter two of which generated with `samtools faidx wheat_v1.0.fa; cut -f1,2 wheat_v1.0.fa.fai > wheat_v1.0.fa.sizes`), each located in `data/index/`. The [IWGSC RefSeq v1.0 Chinese Spring genome assembly](https://urgi.versailles.inra.fr/download/iwgsc/IWGSC_RefSeq_Assemblies/v1.0/) was used in this case
- A reference genome index for bowtie2, located in `data/index/` (see `data/index/bowtie2_index.sh` for an example of how to generate this)
- `Snakefile` in this repository. This contains "rules" that each execute a step in the workflow
- `config.yaml` in this repository. This contains customizable parameters including `reference`, which should be the reference genome file name without the `.fa` extension (e.g., `wheat_v1.0`)
- Optional: `environment.yaml` in this repository, used to create the software environment if conda is used
- If conda is not used, the tools listed in environment.yaml must be specified in the PATH variable

### Creating the conda environment

```
conda env create --file environment.yaml --name ChIPseq_mapping
```

### Usage

In a Unix shell, navigate to the base directory containing `Snakefile`, `config.yaml`, `environment.yaml`, and the `data\` subdirectory, which should have a directory tree structure like this:

```
.
├── alignment_summary.sh
├── config.yaml
├── data
│   ├── DMC1_Rep1_ChIP_R1.fastq.gz
│   ├── DMC1_Rep1_ChIP_R2.fastq.gz
│   └── index
│       ├── bowtie2_index.sh
│       ├── samtools_faidx_chr_sizes.sh
│       ├── wheat_v1.0.1.bt2l
│       ├── wheat_v1.0.2.bt2l
│       ├── wheat_v1.0.3.bt2l
│       ├── wheat_v1.0.4.bt2l
│       ├── wheat_v1.0.fa
│       ├── wheat_v1.0.fa.fai
│       ├── wheat_v1.0.fa.sizes
│       ├── wheat_v1.0.rev.1.bt2l
│       └── wheat_v1.0.rev.2.bt2l
├── environment.yaml
├── README.md
├── scripts
│   └── keepPaired.py
└── Snakefile
```

Then run the following commands in the base directory (`--cores` should match the `THREADS` parameter in `config.yaml`):

```
conda activate ChIPseq_mapping
snakemake -p --cores 48
conda deactivate
```

### Useful Snakemake parameters

- `--cores` specifies the maximum number of threads
- `-n` performs a dry run
- `-p` prints commands
- `--use-conda`
- `--conda-prefix ~/.myconda`
- `--forcerun calc_coverage` forces rerun of a given rule (e.g., `calc_coverage`)

### Outputs

Below is the directory tree structure including files generated once the Snakemake workflow has run to completion.

```
.
├── alignment_summary.sh
├── config.yaml
├── data
│   ├── dedup
│   │   ├── DMC1_Rep1_ChIP_dedup_singletons.fastq.gz
│   │   ├── DMC1_Rep1_ChIP_R1_dedup.fastq.gz
│   │   ├── DMC1_Rep1_ChIP_R1_dedup_repair.fastq.gz
│   │   ├── DMC1_Rep1_ChIP_R2_dedup_repair.fastq.gz
│   │   └── trimmed
│   ├── DMC1_Rep1_ChIP_R1.fastq.gz
│   ├── DMC1_Rep1_ChIP_R2.fastq.gz
│   └── index
│       ├── wheat_v1.0.1.bt2l
│       ├── wheat_v1.0.2.bt2l
│       ├── wheat_v1.0.3.bt2l
│       ├── wheat_v1.0.4.bt2l
│       ├── wheat_v1.0.fa
│       ├── wheat_v1.0.fa.fai
│       ├── wheat_v1.0.fa.sizes
│       ├── wheat_v1.0.rev.1.bt2l
│       └── wheat_v1.0.rev.2.bt2l
├── environment.yaml
├── logs
│   ├── alignment_stats
│   │   ├── DMC1_Rep1_ChIP_alignment_summary.txt
│   │   ├── DMC1_Rep1_ChIP_MappedOn_wheat_v1.0.bam.stats
│   │   ├── DMC1_Rep1_ChIP_MappedOn_wheat_v1.0_lowXM_both_sort.bam.stats
│   │   ├── DMC1_Rep1_ChIP_MappedOn_wheat_v1.0_lowXM_unique_sort.bam.stats
│   │   ├── DMC1_Rep1_ChIP_R1_dedup_repair.fastq.gz.stats
│   │   ├── DMC1_Rep1_ChIP_R1_dedup_repair_trimmed.fastq.gz.stats
│   │   ├── DMC1_Rep1_ChIP_R1.fastq.gz.stats
│   │   ├── DMC1_Rep1_ChIP_R2_dedup_repair.fastq.gz.stats
│   │   ├── DMC1_Rep1_ChIP_R2_dedup_repair_trimmed.fastq.gz.stats
│   │   └── DMC1_Rep1_ChIP_R2.fastq.gz.stats
│   ├── bamCoverage
│   │   ├── DMC1_Rep1_ChIP_MappedOn_wheat_v1.0_lowXM_both_sort_norm_binSize1Mb.log
│   │   ├── DMC1_Rep1_ChIP_MappedOn_wheat_v1.0_lowXM_both_sort_norm.log
│   │   ├── DMC1_Rep1_ChIP_MappedOn_wheat_v1.0_lowXM_unique_sort_norm_binSize1Mb.log
│   │   └── DMC1_Rep1_ChIP_MappedOn_wheat_v1.0_lowXM_unique_sort_norm.log
│   ├── bowtie2
│   │   └── DMC1_Rep1_ChIP_MappedOn_wheat_v1.0_sort.log
│   ├── cutadapt
│   │   └── DMC1_Rep1_ChIP_dedup_repair_trimmed.log
│   ├── dedup
│   │   └── DMC1_Rep1_ChIP_R1_dedup.log
│   ├── fastqc
│   │   ├── raw
│   │   └── trimmed
│   ├── repair
│   │   └── DMC1_Rep1_ChIP_R1_R2_dedup_repair.log
│   └── samtools
│       ├── DMC1_Rep1_ChIP_MappedOn_wheat_v1.0_lowXM_both_sort.log
│       ├── DMC1_Rep1_ChIP_MappedOn_wheat_v1.0_lowXM_unique_sort.log
│       └── stats
├── mapped
│   ├── both
│   │   ├── bg
│   │   │   ├── DMC1_Rep1_ChIP_MappedOn_wheat_v1.0_lowXM_both_sort_norm.bedgraph
│   │   │   └── DMC1_Rep1_ChIP_MappedOn_wheat_v1.0_lowXM_both_sort_norm_binSize1Mb.bedgraph
│   │   ├── bw
│   │   │   └── DMC1_Rep1_ChIP_MappedOn_wheat_v1.0_lowXM_both_sort_norm.bw
│   │   ├── DMC1_Rep1_ChIP_MappedOn_wheat_v1.0_lowXM_both_sort.bam
│   │   └── DMC1_Rep1_ChIP_MappedOn_wheat_v1.0_lowXM_both_sort.bam.csi
│   └── unique
│       ├── bg
│       │   ├── DMC1_Rep1_ChIP_MappedOn_wheat_v1.0_lowXM_both_sort_norm.bedgraph
│       │   └── DMC1_Rep1_ChIP_MappedOn_wheat_v1.0_lowXM_both_sort_norm_binSize1Mb.bedgraph
│       ├── bw
│       │   └── DMC1_Rep1_ChIP_MappedOn_wheat_v1.0_lowXM_both_sort_norm.bw
│       ├── DMC1_Rep1_ChIP_MappedOn_wheat_v1.0_lowXM_both_sort.bam
│       └── DMC1_Rep1_ChIP_MappedOn_wheat_v1.0_lowXM_both_sort.bam.csi
├── scripts
│   └── keepPaired.py
├── README.md
└── Snakefile
```

### Updating the conda environment

```
conda env update --file environment.yaml --name ChIPseq_mapping
```
