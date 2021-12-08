# Coverage, SNP frequency and TE frequency profiles around ChIP-seq peaks

This subdirectoy contains a Snakemake workflow for creating a matrix of windowed coverage values within ChIP-seq peaks and in flanking regions.

It also contains R scripts for creating a matrix of windowed SNP (`1000exomes_SNP_profiles_around_peaks_commandArgs.R`) or transposable element (`TE_superfamily_profiles_around_peaks_commandArgs.R`) frequency values within ChIP-seq peaks and in flanking regions.

### Requirements

- Installation of [Snakemake](https://snakemake.readthedocs.io/en/stable/) and optionally [conda](https://conda.io/docs/)
- `Snakefile` in this repository. This contains "rules" that each execute a step in the workflow
- `config.yaml` in this repository. This contains customizable parameters including `reference`, which should be the reference genome file name without the `.fa` extension (e.g., `wheat_v1.0`)
- Optional: `environment.yaml` in [scripts/read_alignment/snakemake_ChIPseq_MNaseseq/](https://github.com/ajtock/Wheat_DMC1_ASY1_paper/tree/master/scripts/read_alignment/snakemake_ChIPseq_MNaseseq/), used to create the software environment if conda is used
- If conda is not used, [deepTools](https://deeptools.readthedocs.io/en/develop/) must be installed and specified in the PATH variable
- Peak coordinates and, separately, random locus coordinates in BED6 format: column 1 = chromosome ID; column 2 = 0-based start coordinates; column 3 = 1-based end coordinates; column 4 = sequential or otherwise unique numbers (this speeds up computation; see comment from dpryan79 on 13/09/2018 under GitHub issue [computeMatrix has problem with multi processors #760](https://github.com/deeptools/deepTools/issues/760)]); column 5 = fill with NA; column 6 = fill with \*
- A bigWig coverage file (generated using deepTools bamCoverage as part of the [snakemake_ChIPseq_MNaseseq/](https://github.com/ajtock/Wheat_DMC1_ASY1_paper/tree/master/scripts/read_alignment/snakemake_ChIPseq_MNaseseq/) pipeline), to be used for calculating coverage profiles around peaks and random loci (e.g., `DMC1_Rep1_ChIP_MappedOn_wheat_v1.0_lowXM_both_sort_norm.bw`)
- A variant call format (VCF) file containing ~3 million exome sequencing-derived SNP sites ([all.GP08_mm75_het3_publication01142019.vcf](http://wheatgenomics.plantpath.ksu.edu/1000EC/)), from [He et al. (2019) *Nat. Genet.* **51**. DOI: 10.1038/s41588-019-0382-2](https://www.nature.com/articles/s41588-019-0382-2)
- Transposable elements (TEs) from the [IWGSC RefSeq v1.0 annotation](https://urgi.versailles.inra.fr/download/iwgsc/IWGSC_RefSeq_Annotations/v1.0/), with genomic coordinates for elements in each of 14 superfamilies in BED6 format, including strand information (one BED6-format file for each TE superfamily and for each set of randomly positioned loci)

### Creating the conda environment

```
conda env create --file environment.yaml --name ChIPseq_mapping
```

### Usage

In a Unix shell, navigate to the base directory containing `Snakefile` and  `config.yaml`.
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

### Updating the conda environment

```
conda env update --file environment.yaml --name ChIPseq_mapping
```
