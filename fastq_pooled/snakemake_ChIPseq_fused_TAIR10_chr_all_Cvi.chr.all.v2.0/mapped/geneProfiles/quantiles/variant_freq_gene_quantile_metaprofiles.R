#!/applications/R/R-4.0.0/bin/Rscript

# author: Andy Tock
# contact: ajt200@cam.ac.uk
# date: 02.08.2021

# Calculate and plot metaprofiles of ChIP-seq
# (feature windowed means and 95% confidence intervals, CIs)
# for each group of features, defined either by
# decreasing varType in orderRegion or randomly
# Features were divided into quantiles based on varType frequency per kb
# in a given feature region (e.g., promoters)

# Usage:
# /applications/R/R-4.0.0/bin/Rscript variant_freq_gene_quantile_metaprofiles.R SNP_INDEL '/data/public_data/arabidopsis/MPIPZ_Jiao_Schneeberger_2020_NatCommun/Cvi' unique 'Chr1,Chr2,Chr3,Chr4,Chr5' 2200 2000 2kb '2 kb' 10 10bp genes genomewide 6 fused_TAIR10_chr_all_Cvi.chr.all.v2.0 '0.02,0.96'

#varType <- "SNP_INDEL"
#dirName <- "/data/public_data/arabidopsis/MPIPZ_Jiao_Schneeberger_2020_NatCommun/Cvi"
#align <- "unique"
#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split = ","))
#bodyLength <- 2200
#upstream <- 2000
#downstream <- 2000
#flankName <- "2kb"
#flankNamePlot <- "2 kb"
#binSize <- 10
#binName <- "10bp"
#orderRegion <- "genes"
#genomeRegion <- "genomewide"
#quantiles <- 6
#refbase <- "fused_TAIR10_chr_all_Cvi.chr.all.v2.0"
## top left
#legendPos <- as.numeric(unlist(strsplit("0.02,0.96",
#                                        split = ",")))
## top centre
#legendPos <- as.numeric(unlist(strsplit("0.38,0.96",
#                                        split = ",")))
## top right
#legendPos <- as.numeric(unlist(strsplit("0.75,0.96",
#                                        split = ",")))
## bottom left
#legendPos <- as.numeric(unlist(strsplit("0.02,0.30",
#                                        split = ",")))

args <- commandArgs(trailingOnly = T)
varType <- args[1]
dirName <- args[2]
align <- args[3]
chrName <- unlist(strsplit(args[4],
                               split = ","))
bodyLength <- as.numeric(args[5])
upstream <- as.numeric(args[6])
downstream <- as.numeric(args[6])
flankName <- args[7]
flankNamePlot <- args[8]
binSize <- as.numeric(args[9])
binName <- args[10]
orderRegion <- args[11]
genomeRegion <- args[12]
quantiles <- as.numeric(args[13])
refbase <- args[14]
legendPos <- as.numeric(unlist(strsplit(args[15],
                                        split = ",")))

options(stringsAsFactors = F)
library(data.table)
library(parallel)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(grid)
library(gridExtra)
library(extrafont)
library(viridis)

outDir <- paste0(paste0(chrName, collapse = "_"),
                 "/quantiles_", genomeRegion, "_by_", varType,
                 "_freq_in_", orderRegion, "/")
plotDir <- paste0(outDir, "plots/")
system(paste0("[ -d ", outDir, " ] || mkdir -p ", outDir))
system(paste0("[ -d ", plotDir, " ] || mkdir -p ", plotDir))

# Define plot titles
featureAcc1NamePlot <- paste0(varType, " quantiles (Col)")
ranFeatAcc1NamePlot <- "Random quantiles (Col)"
ranLocAcc1NamePlot <- "RanLoc quantiles (Col)"
featureAcc2NamePlot <- paste0(varType, " quantiles (", substr(x = refbase, start = 22, stop = 24), ")")
ranFeatAcc2NamePlot <- paste0("Random quantiles (", substr(x = refbase, start = 22, stop = 24), ")")
ranLocAcc2NamePlot <- paste0("RanLoc quantiles (", substr(x = refbase, start = 22, stop = 24), ")")

# Define quantile colours
quantileColours <- rev(viridis(quantiles))
quantileColours <- rev(plasma(quantiles))
quantileColours[1] <- "gold"

# Define feature start and end labels for plotting
featureStartLab <- "TSS"
featureEndLab <- "TTS"

# Load table of featuresAcc1_ortho_DF grouped into quantiles
featuresAcc1_ortho_DF <- fread(
  paste0(outDir,
         "featuresAcc1_ortho_", quantiles, "quantiles",
         "_", genomeRegion, "_by_", varType, "_freq_in_", orderRegion,
         "_of_Acc1_Chr_genes_in_", refbase, "_",
         paste0(chrName, collapse = "_"), ".tsv"),
  header = TRUE, sep = "\t")
# Load features to confirm feature (row) ordering in "featuresAcc1_ortho_DF" is the same
# as in "features" (which was used for generating the coverage matrices)
featuresAcc1_ortho_DF_bed <- read.table(
  paste0(outDir,
         "featuresAcc1_ortho_", quantiles, "quantiles",
         "_", genomeRegion, "_by_", varType, "_freq_in_", orderRegion,
         "_of_Acc1_Chr_genes_in_", refbase, "_",
         paste0(chrName, collapse = "_"), ".bed"),
  header = FALSE, sep = "\t")
colnames(featuresAcc1_ortho_DF_bed) <- c("chr", "start", "end", "ID", "score", "strand")

stopifnot(identical(as.character(featuresAcc1_ortho_DF$ID),
                    as.character(featuresAcc1_ortho_DF_bed$ID)))
rm(featuresAcc1_ortho_DF_bed); gc()

# Get row indices for each feature quantile
quantileIndicesAcc1 <- lapply(1:quantiles, function(k) {
  which(featuresAcc1_ortho_DF$quantile == paste0("Quantile ", k))
})

## Genomic definitions
#fai <- read.table(paste0(dirName, "/", refbase, "/", refbase, ".fa.fai"), header = F)
#chrs <- fai$V1[which(gsub(pattern = ".+_", replacement = "", x = fai$V1) %in% chrName)]
#chrLens <- fai$V2[which(gsub(pattern = ".+_", replacement = "", x = fai$V1) %in% chrName)]

## Random feature quantiles
# Define function to randomly select n rows from
# a data.frame
selectRandomFeatures <- function(features, n) {
  return(features[sample(x = dim(features)[1],
                         size = n,
                         replace = FALSE),])
}

# Divide features into random sets of equal number,
# with the same number of CEN180s per chromosome as
# above-defined varType-defined feature quantiles
chrs <- sort(unique(featuresAcc1_ortho_DF$seqid))
randomPCIndicesAcc1 <- lapply(1:quantiles, function(k) {
  randomPCIndicesAcc1k <- NULL
  for(i in 1:length(chrs)) {
  # Define seed so that random selections are reproducible
  set.seed(93750174)
    randomPCfeatureskChr <- selectRandomFeatures(features = featuresAcc1_ortho_DF[featuresAcc1_ortho_DF$seqid == chrs[i],],
                                                 n = dim(featuresAcc1_ortho_DF[featuresAcc1_ortho_DF$quantile == paste0("Quantile ", k) &
                                                                    featuresAcc1_ortho_DF$seqid == chrs[i],])[1])
    randomPCIndicesAcc1k <- c(randomPCIndicesAcc1k, as.integer(rownames(randomPCfeatureskChr)))
  }
  randomPCIndicesAcc1k
})
# Confirm per-chromosome feature numbers are the same for quantiles and random groupings
lapply(seq_along(1:quantiles), function(k) {
  sapply(seq_along(chrs), function(x) {
    if(!identical(dim(featuresAcc1_ortho_DF[randomPCIndicesAcc1[[k]],][featuresAcc1_ortho_DF[randomPCIndicesAcc1[[k]],]$chr == chrs[x],]),
                  dim(featuresAcc1_ortho_DF[quantileIndicesAcc1[[k]],][featuresAcc1_ortho_DF[quantileIndicesAcc1[[k]],]$chr == chrs[x],])))    {
      stop("Quantile features and random features do not consist of the same number of features per chromosome")
    }
  })
})


# Load feature matrices for each chromatin dataset, calculate log2(ChIP/control),
# and sort by decreasing log2mat1RegionRowMeans
ChIPNames <- c(
               "Col_DMC1_V5_Rep1_ChIP",
               "Col_DMC1_V5_Rep2_ChIP",
               "ColCviF1_DMC1_V5_Rep1_ChIP",
               "ColLerF1_DMC1_V5_Rep1_ChIP",
               "ColColF1_DMC1_V5_Rep1_ChIP",
               "Col_DMC1_V5_Rep1_leaf",
               "Col_DMC1_V5_Rep1_mock",
               "Col_DMC1_V5_Rep2_mock",
               "Col_REC8_HA_Rep2_ChIP"
              )
ChIPNamesDir <- c(
                  rep(paste0("20190917_dh580_Athaliana_ChIPseq_DMC1/fastq_pooled/snakemake_ChIPseq_", refbase), length(ChIPNames))
                 )
log2ChIPNamesPlot <- c(
                       "Col DMC1 Rep1",
                       "Col DMC1 Rep2",
                       "Col/Cvi DMC1 Rep1",
                       "Col/Ler DMC1 Rep1",
                       "Col/Col DMC1 Rep1",
                       "Col DMC1 Leaf",
                       "Col DMC1 Mock1",
                       "Col DMC1 Mock2",
                       "Col REC8-HA Rep2"
                      )
ChIPNamesPlot <- log2ChIPNamesPlot
log2ChIPColours <- c(
                     rep("black", length(log2ChIPNamesPlot))
                    )
ChIPColours <- log2ChIPColours
ChIPDirs <- sapply(seq_along(ChIPNames), function(x) {
  paste0("/home/ajt200/analysis/",
         ChIPNamesDir[x],
         "/mapped/geneProfiles/matrices/")
})

controlNames <- c(
                  "Col_DMC1_V5_Rep1_art75",
                  "Col_DMC1_V5_Rep2_art150",
                  "ColCviF1_DMC1_V5_Rep1_art150",
                  "ColLerF1_DMC1_V5_Rep1_art150",
                  "ColColF1_DMC1_V5_Rep1_art150",
                  "Col_REC8_Myc_Rep1_input"
                 )
controlNamesDir <- c(
                     rep(paste0("20190917_dh580_Athaliana_ChIPseq_DMC1/fastq_pooled/snakemake_ChIPseq_", refbase), length(controlNames))
                    )
controlNamesPlot <- c(
                      "Col DMC1 Rep1 art75",
                      "Col DMC1 Rep2 art150",
                      "Col/Cvi DMC1 Rep1 art150",
                      "Col/Ler DMC1 Rep1 art150",
                      "Col/Col DMC1 Rep1 art150",
                      "Col REC8-Myc Rep1 input"
                     )
controlColours <- c(
                    rep("black", length(controlNamesPlot))
                   )
controlDirs <- sapply(seq_along(controlNames), function(x) {
  paste0("/home/ajt200/analysis/",
         controlNamesDir[x],
         "/mapped/geneProfiles/matrices/")
})

## ChIP
# featureAcc1
ChIP_featureAcc1Mats <- mclapply(seq_along(ChIPNames), function(x) {
  if(grepl("F1", ChIPNames[x])) {
    as.matrix(read.table(paste0(ChIPDirs[x],
                                ChIPNames[x],
                                "_MappedOn_", refbase, "_lowXM_", align, "_sort_TPM_Acc1_Chr_genes_in_",
                                paste0(chrName, collapse = "_"), "_", genomeRegion,
                                "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
  } else {
    as.matrix(read.table(paste0(ChIPDirs[x],
                                ChIPNames[x],
                                "_MappedOn_", refbase, "_lowXM_", align, "_sort_TPM_Acc1_Chr_genes_in_",
                                paste0(chrName, collapse = "_"), "_", genomeRegion,
                                "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))/2
 }
}, mc.cores = length(ChIPNames))

# ranLocAcc1
ChIP_ranLocAcc1Mats <- mclapply(seq_along(ChIPNames), function(x) {
  if(grepl("F1", ChIPNames[x])) {
    as.matrix(read.table(paste0(ChIPDirs[x],
                                ChIPNames[x],
                                "_MappedOn_", refbase, "_lowXM_", align, "_sort_TPM_Acc1_Chr_genes_in_",
                                paste0(chrName, collapse = "_"), "_", genomeRegion,
                                "_ranLocAcc1_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
  } else {
    as.matrix(read.table(paste0(ChIPDirs[x],
                                ChIPNames[x],
                                "_MappedOn_", refbase, "_lowXM_", align, "_sort_TPM_Acc1_Chr_genes_in_",
                                paste0(chrName, collapse = "_"), "_", genomeRegion,
                                "_ranLocAcc1_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))/2
 }
}, mc.cores = length(ChIPNames))

# featureAcc2
ChIP_featureAcc2Mats <- mclapply(seq_along(ChIPNames), function(x) {
  if(grepl("F1", ChIPNames[x])) {
    as.matrix(read.table(paste0(ChIPDirs[x],
                                ChIPNames[x],
                                "_MappedOn_", refbase, "_lowXM_", align, "_sort_TPM_Acc2_Chr_genes_in_",
                                paste0(chrName, collapse = "_"), "_", genomeRegion,
                                "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
  } else {
    as.matrix(read.table(paste0(ChIPDirs[x],
                                ChIPNames[x],
                                "_MappedOn_", refbase, "_lowXM_", align, "_sort_TPM_Acc2_Chr_genes_in_",
                                paste0(chrName, collapse = "_"), "_", genomeRegion,
                                "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))/2
 }
}, mc.cores = length(ChIPNames))

# ranLocAcc2
ChIP_ranLocAcc2Mats <- mclapply(seq_along(ChIPNames), function(x) {
  if(grepl("F1", ChIPNames[x])) {
    as.matrix(read.table(paste0(ChIPDirs[x],
                                ChIPNames[x],
                                "_MappedOn_", refbase, "_lowXM_", align, "_sort_TPM_Acc2_Chr_genes_in_",
                                paste0(chrName, collapse = "_"), "_", genomeRegion,
                                "_ranLocAcc2_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
  } else {
    as.matrix(read.table(paste0(ChIPDirs[x],
                                ChIPNames[x],
                                "_MappedOn_", refbase, "_lowXM_", align, "_sort_TPM_Acc2_Chr_genes_in_",
                                paste0(chrName, collapse = "_"), "_", genomeRegion,
                                "_ranLocAcc2_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))/2
 }
}, mc.cores = length(ChIPNames))


### control
## featureAcc1
#control_featureAcc1Mats <- mclapply(seq_along(controlNames), function(x) {
#  if(grepl("F1", controlNames[x])) {
#    as.matrix(read.table(paste0(controlDirs[x],
#                                controlNames[x],
#                                "_MappedOn_", refbase, "_lowXM_", align, "_sort_TPM_Acc1_Chr_genes_in_",
#                                paste0(chrName, collapse = "_"), "_", genomeRegion,
#                                "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
#                         header = F, skip = 3))
#  } else {
#    as.matrix(read.table(paste0(controlDirs[x],
#                                controlNames[x],
#                                "_MappedOn_", refbase, "_lowXM_", align, "_sort_TPM_Acc1_Chr_genes_in_",
#                                paste0(chrName, collapse = "_"), "_", genomeRegion,
#                                "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
#                         header = F, skip = 3))/2
# }
#}, mc.cores = length(controlNames))
#
## ranLocAcc1
#control_ranLocAcc1Mats <- mclapply(seq_along(controlNames), function(x) {
#  if(grepl("F1", controlNames[x])) {
#    as.matrix(read.table(paste0(controlDirs[x],
#                                controlNames[x],
#                                "_MappedOn_", refbase, "_lowXM_", align, "_sort_TPM_Acc1_Chr_genes_in_",
#                                paste0(chrName, collapse = "_"), "_", genomeRegion,
#                                "_ranLocAcc1_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
#                         header = F, skip = 3))
#  } else {
#    as.matrix(read.table(paste0(controlDirs[x],
#                                controlNames[x],
#                                "_MappedOn_", refbase, "_lowXM_", align, "_sort_TPM_Acc1_Chr_genes_in_",
#                                paste0(chrName, collapse = "_"), "_", genomeRegion,
#                                "_ranLocAcc1_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
#                         header = F, skip = 3))/2
# }
#}, mc.cores = length(controlNames))
#
## featureAcc2
#control_featureAcc2Mats <- mclapply(seq_along(controlNames), function(x) {
#  if(grepl("F1", controlNames[x])) {
#    as.matrix(read.table(paste0(controlDirs[x],
#                                controlNames[x],
#                                "_MappedOn_", refbase, "_lowXM_", align, "_sort_TPM_Acc2_Chr_genes_in_",
#                                paste0(chrName, collapse = "_"), "_", genomeRegion,
#                                "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
#                         header = F, skip = 3))
#  } else {
#    as.matrix(read.table(paste0(controlDirs[x],
#                                controlNames[x],
#                                "_MappedOn_", refbase, "_lowXM_", align, "_sort_TPM_Acc2_Chr_genes_in_",
#                                paste0(chrName, collapse = "_"), "_", genomeRegion,
#                                "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
#                         header = F, skip = 3))/2
# }
#}, mc.cores = length(controlNames))
#
## ranLocAcc2
#control_ranLocAcc2Mats <- mclapply(seq_along(controlNames), function(x) {
#  if(grepl("F1", controlNames[x])) {
#    as.matrix(read.table(paste0(controlDirs[x],
#                                controlNames[x],
#                                "_MappedOn_", refbase, "_lowXM_", align, "_sort_TPM_Acc2_Chr_genes_in_",
#                                paste0(chrName, collapse = "_"), "_", genomeRegion,
#                                "_ranLocAcc2_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
#                         header = F, skip = 3))
#  } else {
#    as.matrix(read.table(paste0(controlDirs[x],
#                                controlNames[x],
#                                "_MappedOn_", refbase, "_lowXM_", align, "_sort_TPM_Acc2_Chr_genes_in_",
#                                paste0(chrName, collapse = "_"), "_", genomeRegion,
#                                "_ranLocAcc2_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
#                         header = F, skip = 3))/2
# }
#}, mc.cores = length(controlNames))


## Conditionally calculate log2(ChIP/control)
## for each matrix depending on library
## feature
#log2ChIP_featureAcc1Mats <- mclapply(seq_along(ChIP_featureAcc1Mats), function(x) {
#  if ( grepl("MNase", ChIPNames[x]) ) {
#    print(paste0(ChIPNames[x], " library; using ", controlNames[2], " for log2((MNase+1)/(gDNA+1)) calculation"))
#    log2((ChIP_featureAcc1Mats[[x]]+1)/(control_featureAcc1Mats[[2]]+1))
#  } else if ( grepl("SPO11oligos", ChIPNames[x]) ) {
#    print(paste0(ChIPNames[x], " library; using ", controlNames[3], " for log2((SPO11-1-oligos+1)/(gDNA+1)) calculation"))
#    log2((ChIP_featureAcc1Mats[[x]]+1)/(control_featureAcc1Mats[[3]]+1))
#  } else if ( grepl("set1", ChIPNames[x]) ) {
#    print(paste0(ChIPNames[x], " library; using ", controlNames[4], " for log2((ChIP+1)/(input+1)) calculation"))
#    log2((ChIP_featureAcc1Mats[[x]]+1)/(control_featureAcc1Mats[[4]]+1))
#  } else if ( grepl("set2", ChIPNames[x]) ) {
#    print(paste0(ChIPNames[x], " library; using ", controlNames[5], " for log2((ChIP+1)/(input+1)) calculation"))
#    log2((ChIP_featureAcc1Mats[[x]]+1)/(control_featureAcc1Mats[[5]]+1))
#  } else if ( grepl("set3", ChIPNames[x]) ) {
#    print(paste0(ChIPNames[x], " library; using ", controlNames[6], " for log2((ChIP+1)/(input+1)) calculation"))
#    log2((ChIP_featureAcc1Mats[[x]]+1)/(control_featureAcc1Mats[[6]]+1))
#  } else if ( grepl("set5", ChIPNames[x]) ) {
#    print(paste0(ChIPNames[x], " library; using ", controlNames[7], " for log2((ChIP+1)/(input+1)) calculation"))
#    log2((ChIP_featureAcc1Mats[[x]]+1)/(control_featureAcc1Mats[[7]]+1))
#  } else if ( grepl("set7", ChIPNames[x]) ) {
#    print(paste0(ChIPNames[x], " library; using ", controlNames[8], " for log2((ChIP+1)/(input+1)) calculation"))
#    log2((ChIP_featureAcc1Mats[[x]]+1)/(control_featureAcc1Mats[[8]]+1))
#  } else {
#    print(paste0(ChIPNames[x], " library; using ", controlNames[1], " for log2((ChIP+1)/(input+1)) calculation"))
#    log2((ChIP_featureAcc1Mats[[x]]+1)/(control_featureAcc1Mats[[1]]+1))
#  }
#}, mc.cores = length(ChIP_featureAcc1Mats))
#
## Conditionally calculate log2(ChIP/control)
## for each matrix depending on library
## ranLoc
#log2ChIP_ranLocAcc1Mats <- mclapply(seq_along(ChIP_ranLocAcc1Mats), function(x) {
#  if ( grepl("MNase", ChIPNames[x]) ) {
#    print(paste0(ChIPNames[x], " library; using ", controlNames[2], " for log2((MNase+1)/(gDNA+1)) calculation"))
#    log2((ChIP_ranLocAcc1Mats[[x]]+1)/(control_ranLocAcc1Mats[[2]]+1))
#  } else if ( grepl("SPO11oligos", ChIPNames[x]) ) {
#    print(paste0(ChIPNames[x], " library; using ", controlNames[3], " for log2((SPO11-1-oligos+1)/(gDNA+1)) calculation"))
#    log2((ChIP_ranLocAcc1Mats[[x]]+1)/(control_ranLocAcc1Mats[[3]]+1))
#  } else if ( grepl("set1", ChIPNames[x]) ) {
#    print(paste0(ChIPNames[x], " library; using ", controlNames[4], " for log2((ChIP+1)/(input+1)) calculation"))
#    log2((ChIP_ranLocAcc1Mats[[x]]+1)/(control_ranLocAcc1Mats[[4]]+1))
#  } else if ( grepl("set2", ChIPNames[x]) ) {
#    print(paste0(ChIPNames[x], " library; using ", controlNames[5], " for log2((ChIP+1)/(input+1)) calculation"))
#    log2((ChIP_ranLocAcc1Mats[[x]]+1)/(control_ranLocAcc1Mats[[5]]+1))
#  } else if ( grepl("set3", ChIPNames[x]) ) {
#    print(paste0(ChIPNames[x], " library; using ", controlNames[6], " for log2((ChIP+1)/(input+1)) calculation"))
#    log2((ChIP_ranLocAcc1Mats[[x]]+1)/(control_ranLocAcc1Mats[[6]]+1))
#  } else if ( grepl("set5", ChIPNames[x]) ) {
#    print(paste0(ChIPNames[x], " library; using ", controlNames[7], " for log2((ChIP+1)/(input+1)) calculation"))
#    log2((ChIP_ranLocAcc1Mats[[x]]+1)/(control_ranLocAcc1Mats[[7]]+1))
#  } else if ( grepl("set7", ChIPNames[x]) ) {
#    print(paste0(ChIPNames[x], " library; using ", controlNames[8], " for log2((ChIP+1)/(input+1)) calculation"))
#    log2((ChIP_ranLocAcc1Mats[[x]]+1)/(control_ranLocAcc1Mats[[8]]+1))
#  } else {
#    print(paste0(ChIPNames[x], " library; using ", controlNames[1], " for log2((ChIP+1)/(input+1)) calculation"))
#    log2((ChIP_ranLocAcc1Mats[[x]]+1)/(control_ranLocAcc1Mats[[1]]+1))
#  }
#}, mc.cores = length(ChIP_ranLocAcc1Mats))
#
## log2ChIP
## Add column names
#for(x in seq_along(log2ChIP_featureAcc1Mats)) {
#  colnames(log2ChIP_featureAcc1Mats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
#                                           paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
#                                           paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
#  colnames(log2ChIP_ranLocAcc1Mats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
#                                          paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
#                                          paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
#}
#
## Subdivide coverage matrices into above-defined quantiles and random groupings
#log2ChIP_mats_quantiles <- mclapply(seq_along(log2ChIP_featureAcc1Mats), function(x) {
#  list(
#       # feature quantiles
#       lapply(1:quantiles, function(k) {
#         log2ChIP_featureAcc1Mats[[x]][quantileIndicesAcc1[[k]],]
#       }),
#       # feature random groupings
#       lapply(1:quantiles, function(k) {
#         log2ChIP_featureAcc1Mats[[x]][randomPCIndicesAcc1[[k]],]
#       }),
#       # random loci groupings
#       lapply(1:quantiles, function(k) {
#         log2ChIP_ranLocAcc1Mats[[x]][quantileIndicesAcc1[[k]],]
#       })
#      ) 
#}, mc.cores = length(log2ChIP_featureAcc1Mats))
#
## Transpose matrix and convert into dataframe
## in which first column is window name
#wideDFfeature_list_log2ChIP <- mclapply(seq_along(log2ChIP_mats_quantiles), function(x) {
#  lapply(seq_along(log2ChIP_mats_quantiles[[x]]), function(y) {
#    lapply(seq_along(log2ChIP_mats_quantiles[[x]][[y]]), function(k) {
#      data.frame(window = colnames(log2ChIP_mats_quantiles[[x]][[y]][[k]]),
#                 t(log2ChIP_mats_quantiles[[x]][[y]][[k]]))
#    })
#  })
#}, mc.cores = length(log2ChIP_mats_quantiles))
#
## Convert into tidy data.frame (long format)
#tidyDFfeature_list_log2ChIP  <- mclapply(seq_along(wideDFfeature_list_log2ChIP), function(x) {
#  lapply(seq_along(log2ChIP_mats_quantiles[[x]]), function(y) {
#    lapply(seq_along(log2ChIP_mats_quantiles[[x]][[y]]), function(k) {
#      gather(data  = wideDFfeature_list_log2ChIP[[x]][[y]][[k]],
#             key   = feature,
#             value = coverage,
#             -window)
#    })
#  }) 
#}, mc.cores = length(wideDFfeature_list_log2ChIP))
#
## Order levels of factor "window" so that sequential levels
## correspond to sequential windows
#for(x in seq_along(tidyDFfeature_list_log2ChIP)) {
#  for(y in seq_along(log2ChIP_mats_quantiles[[x]])) {
#    for(k in seq_along(log2ChIP_mats_quantiles[[x]][[y]])) {
#      tidyDFfeature_list_log2ChIP[[x]][[y]][[k]]$window <- factor(tidyDFfeature_list_log2ChIP[[x]][[y]][[k]]$window,
#                                                                  levels = as.character(wideDFfeature_list_log2ChIP[[x]][[y]][[k]]$window))
#    }
#  }
#}
#
## Create summary data.frame in which each row corresponds to a window (Column 1),
## Column2 is the number of coverage values (features) per window,
## Column3 is the mean of coverage values per window,
## Column4 is the standard deviation of coverage values per window,
## Column5 is the standard error of the mean of coverage values per window,
## Column6 is the lower bound of the 95% confidence interval, and
## Column7 is the upper bound of the 95% confidence interval
#summaryDFfeature_list_log2ChIP  <- mclapply(seq_along(tidyDFfeature_list_log2ChIP), function(x) {
#  lapply(seq_along(log2ChIP_mats_quantiles[[x]]), function(y) {
#    lapply(seq_along(log2ChIP_mats_quantiles[[x]][[y]]), function(k) {
#      data.frame(window = as.character(wideDFfeature_list_log2ChIP[[x]][[y]][[k]]$window),
#                 n      = tapply(X     = tidyDFfeature_list_log2ChIP[[x]][[y]][[k]]$coverage,
#                                 INDEX = tidyDFfeature_list_log2ChIP[[x]][[y]][[k]]$window,
#                                 FUN   = length),
#                 mean   = tapply(X     = tidyDFfeature_list_log2ChIP[[x]][[y]][[k]]$coverage,
#                                 INDEX = tidyDFfeature_list_log2ChIP[[x]][[y]][[k]]$window,
#                                 FUN   = mean,
#                                 na.rm = TRUE),
#                 sd     = tapply(X     = tidyDFfeature_list_log2ChIP[[x]][[y]][[k]]$coverage,
#                                 INDEX = tidyDFfeature_list_log2ChIP[[x]][[y]][[k]]$window,
#                                 FUN   = sd,
#                                 na.rm = TRUE))
#    })
#  })
#}, mc.cores = length(tidyDFfeature_list_log2ChIP))
#
#for(x in seq_along(summaryDFfeature_list_log2ChIP)) {
#  for(y in seq_along(log2ChIP_mats_quantiles[[x]])) {
#    for(k in seq_along(log2ChIP_mats_quantiles[[x]][[y]])) {
#      summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$window <- factor(summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$window,
#                                                                     levels = as.character(wideDFfeature_list_log2ChIP[[x]][[y]][[k]]$window))
#      summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$winNo <- factor(1:dim(summaryDFfeature_list_log2ChIP[[x]][[y]][[k]])[1])
#      summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$sem <- summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$sd/sqrt(summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$n-1)
#      summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$CI_lower <- summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$mean -
#        qt(0.975, df = summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$n-1)*summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$sem
#      summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$CI_upper <- summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$mean +
#        qt(0.975, df = summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$n-1)*summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$sem
#    }
#  }
#}
#
#quantileNames <- paste0(rep("Quantile ", quantiles), 1:quantiles)
#randomPCNames <- paste0(rep("Random ", quantiles), 1:quantiles)
#for(x in seq_along(summaryDFfeature_list_log2ChIP)) {
#  # feature quantiles
#  names(summaryDFfeature_list_log2ChIP[[x]][[1]]) <- quantileNames
#  # feature random groupings
#  names(summaryDFfeature_list_log2ChIP[[x]][[2]]) <- randomPCNames
#  # random loci groupings
#  names(summaryDFfeature_list_log2ChIP[[x]][[3]]) <- randomPCNames
#}
#
## Convert list of lists of lists of feature quantiles summaryDFfeature_list_log2ChIP into
## a list of lists of single data.frames containing all feature quantiles for plotting
#summaryDFfeature_log2ChIP  <- mclapply(seq_along(summaryDFfeature_list_log2ChIP), function(x) {
#  lapply(seq_along(log2ChIP_mats_quantiles[[x]]), function(y) {
#    bind_rows(summaryDFfeature_list_log2ChIP[[x]][[y]], .id = "quantile")
#  })
#}, mc.cores = length(summaryDFfeature_list_log2ChIP))
#for(x in seq_along(summaryDFfeature_log2ChIP)) {
#  # feature quantiles
#  summaryDFfeature_log2ChIP[[x]][[1]]$quantile <- factor(summaryDFfeature_log2ChIP[[x]][[1]]$quantile,
#                                                         levels = names(summaryDFfeature_list_log2ChIP[[x]][[1]]))
#  # feature random groupings
#  summaryDFfeature_log2ChIP[[x]][[2]]$quantile <- factor(summaryDFfeature_log2ChIP[[x]][[2]]$quantile,
#                                                         levels = names(summaryDFfeature_list_log2ChIP[[x]][[2]]))
#  # random loci groupings
#  summaryDFfeature_log2ChIP[[x]][[3]]$quantile <- factor(summaryDFfeature_log2ChIP[[x]][[3]]$quantile,
#                                                         levels = names(summaryDFfeature_list_log2ChIP[[x]][[3]]))
#}
#
## Define y-axis limits
#ymin_list_log2ChIP <- lapply(seq_along(summaryDFfeature_log2ChIP), function(x) {
#  min(c(summaryDFfeature_log2ChIP[[x]][[1]]$CI_lower,
#        summaryDFfeature_log2ChIP[[x]][[2]]$CI_lower,
#        summaryDFfeature_log2ChIP[[x]][[3]]$CI_lower))
#})
#ymax_list_log2ChIP <- lapply(seq_along(summaryDFfeature_log2ChIP), function(x) {
#  max(c(summaryDFfeature_log2ChIP[[x]][[1]]$CI_upper,
#        summaryDFfeature_log2ChIP[[x]][[2]]$CI_upper,
#        summaryDFfeature_log2ChIP[[x]][[3]]$CI_upper))
#})
#
## Define legend labels
#legendLabs_featureAcc1 <- lapply(seq_along(quantileNames), function(x) {
#  grobTree(textGrob(bquote(.(quantileNames[x])),
#                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
#                    gp = gpar(col = quantileColours[x], fontsize = 18)))
#})
#legendLabs_ranFeatAcc1 <- lapply(seq_along(randomPCNames), function(x) {
#  grobTree(textGrob(bquote(.(randomPCNames[x])),
#                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
#                    gp = gpar(col = quantileColours[x], fontsize = 18)))
#})
#legendLabs_ranLocAcc1 <- lapply(seq_along(randomPCNames), function(x) {
#  grobTree(textGrob(bquote(.(randomPCNames[x])),
#                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
#                    gp = gpar(col = quantileColours[x], fontsize = 18)))
#})
#
## Plot average profiles with 95% CI ribbon
### feature
#ggObj1_combined_log2ChIP <- mclapply(seq_along(log2ChIPNamesPlot), function(x) {
#  summaryDFfeature <- summaryDFfeature_log2ChIP[[x]][[1]]
#  ggplot(data = summaryDFfeature,
#         mapping = aes(x = winNo,
#                       y = mean,
#                       group = quantile)
#        ) +
#  geom_line(data = summaryDFfeature,
#            mapping = aes(colour = quantile),
#            size = 1) +
#  scale_colour_manual(values = quantileColours) +
#  geom_ribbon(data = summaryDFfeature,
#              mapping = aes(ymin = CI_lower,
#                            ymax = CI_upper,
#                            fill = quantile),
#              alpha = 0.4) +
#  scale_fill_manual(values = quantileColours) +
#  scale_y_continuous(limits = c(ymin_list_log2ChIP[[x]], ymax_list_log2ChIP[[x]]),
#                     labels = function(x) sprintf("%6.3f", x)) +
#  scale_x_discrete(breaks = c(1,
#                              (upstream/binSize)+1,
#                              (dim(summaryDFfeature_log2ChIP[[x]][[1]])[1]/quantiles)-(downstream/binSize),
#                              dim(summaryDFfeature_log2ChIP[[x]][[1]])[1]/quantiles),
#                   labels = c(paste0("-", flankName),
#                              featureStartLab,
#                              featureEndLab,
#                              paste0("+", flankName))) +
#  geom_vline(xintercept = c((upstream/binSize)+1,
#                            (dim(summaryDFfeature_log2ChIP[[x]][[1]])[1]/quantiles)-(downstream/binSize)),
#             linetype = "dashed",
#             size = 1) +
#  labs(x = "",
#       y = log2ChIPNamesPlot[x]) +
#  annotation_custom(legendLabs_featureAcc1[[1]]) +
#  annotation_custom(legendLabs_featureAcc1[[2]]) +
#  annotation_custom(legendLabs_featureAcc1[[3]]) +
#  annotation_custom(legendLabs_featureAcc1[[4]]) +
#  theme_bw() +
#  theme(
#        axis.ticks = element_line(size = 1.0, colour = "black"),
#        axis.ticks.length = unit(0.25, "cm"),
#        axis.text.x = element_text(size = 22, colour = "black"),
#        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
#        axis.title = element_text(size = 30, colour = log2ChIPColours[x]),
#        legend.position = "none",
#        panel.grid = element_blank(),
#        panel.border = element_rect(size = 3.5, colour = "black"),
#        panel.background = element_blank(),
#        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
#        plot.title = element_text(hjust = 0.5, size = 30)) +
#  ggtitle(bquote(.(featureAcc1NamePlot) ~ "(" * italic("n") ~ "=" ~
#                 .(prettyNum(summaryDFfeature$n[1],
#                             big.mark = ",", trim = T)) *
#                 ")"))
#}, mc.cores = length(log2ChIPNamesPlot))
#
### ranFeat
#ggObj2_combined_log2ChIP <- mclapply(seq_along(log2ChIPNamesPlot), function(x) {
#  summaryDFfeature <- summaryDFfeature_log2ChIP[[x]][[2]]
#  ggplot(data = summaryDFfeature,
#         mapping = aes(x = winNo,
#                       y = mean,
#                       group = quantile)
#        ) +
#  geom_line(data = summaryDFfeature,
#            mapping = aes(colour = quantile),
#            size = 1) +
#  scale_colour_manual(values = quantileColours) +
#  geom_ribbon(data = summaryDFfeature,
#              mapping = aes(ymin = CI_lower,
#                            ymax = CI_upper,
#                            fill = quantile),
#              alpha = 0.4) +
#  scale_fill_manual(values = quantileColours) +
#  scale_y_continuous(limits = c(ymin_list_log2ChIP[[x]], ymax_list_log2ChIP[[x]]),
#                     labels = function(x) sprintf("%6.3f", x)) +
#  scale_x_discrete(breaks = c(1,
#                              (upstream/binSize)+1,
#                              (dim(summaryDFfeature_log2ChIP[[x]][[2]])[1]/quantiles)-(downstream/binSize),
#                              dim(summaryDFfeature_log2ChIP[[x]][[2]])[1]/quantiles),
#                   labels = c(paste0("-", flankName),
#                              featureStartLab,
#                              featureEndLab,
#                              paste0("+", flankName))) +
#  geom_vline(xintercept = c((upstream/binSize)+1,
#                            (dim(summaryDFfeature_log2ChIP[[x]][[2]])[1]/quantiles)-(downstream/binSize)),
#             linetype = "dashed",
#             size = 1) +
#  labs(x = "",
#       y = log2ChIPNamesPlot[x]) +
#  annotation_custom(legendLabs_ranFeatAcc1[[1]]) +
#  annotation_custom(legendLabs_ranFeatAcc1[[2]]) +
#  annotation_custom(legendLabs_ranFeatAcc1[[3]]) +
#  annotation_custom(legendLabs_ranFeatAcc1[[4]]) +
#  theme_bw() +
#  theme(
#        axis.ticks = element_line(size = 1.0, colour = "black"),
#        axis.ticks.length = unit(0.25, "cm"),
#        axis.text.x = element_text(size = 22, colour = "black"),
#        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
#        axis.title = element_text(size = 30, colour = log2ChIPColours[x]),
#        legend.position = "none",
#        panel.grid = element_blank(),
#        panel.border = element_rect(size = 3.5, colour = "black"),
#        panel.background = element_blank(),
#        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
#        plot.title = element_text(hjust = 0.5, size = 30)) +
#  ggtitle(bquote(.(ranFeatAcc1NamePlot) ~ "(" * italic("n") ~ "=" ~
#                 .(prettyNum(summaryDFfeature$n[1],
#                             big.mark = ",", trim = T)) *
#                 ")"))
#}, mc.cores = length(log2ChIPNamesPlot))
#
### ranLoc
#ggObj3_combined_log2ChIP <- mclapply(seq_along(log2ChIPNamesPlot), function(x) {
#  summaryDFfeature <- summaryDFfeature_log2ChIP[[x]][[3]]
#  ggplot(data = summaryDFfeature,
#         mapping = aes(x = winNo,
#                       y = mean,
#                       group = quantile)
#        ) +
#  geom_line(data = summaryDFfeature,
#            mapping = aes(colour = quantile),
#            size = 1) +
#  scale_colour_manual(values = quantileColours) +
#  geom_ribbon(data = summaryDFfeature,
#              mapping = aes(ymin = CI_lower,
#                            ymax = CI_upper,
#                            fill = quantile),
#              alpha = 0.4) +
#  scale_fill_manual(values = quantileColours) +
#  scale_y_continuous(limits = c(ymin_list_log2ChIP[[x]], ymax_list_log2ChIP[[x]]),
#                     labels = function(x) sprintf("%6.3f", x)) +
#  scale_x_discrete(breaks = c(1,
#                              (upstream/binSize)+1,
#                              (dim(summaryDFfeature_log2ChIP[[x]][[3]])[1]/quantiles)-(downstream/binSize),
#                              dim(summaryDFfeature_log2ChIP[[x]][[3]])[1]/quantiles),
#                   labels = c(paste0("-", flankName),
#                              "Start",
#                              "End",
#                              paste0("+", flankName))) +
#  geom_vline(xintercept = c((upstream/binSize)+1,
#                            (dim(summaryDFfeature_log2ChIP[[x]][[3]])[1]/quantiles)-(downstream/binSize)),
#             linetype = "dashed",
#             size = 1) +
#  labs(x = "",
#       y = log2ChIPNamesPlot[x]) +
#  annotation_custom(legendLabs_ranLocAcc1[[1]]) +
#  annotation_custom(legendLabs_ranLocAcc1[[2]]) +
#  annotation_custom(legendLabs_ranLocAcc1[[3]]) +
#  annotation_custom(legendLabs_ranLocAcc1[[4]]) +
#  theme_bw() +
#  theme(
#        axis.ticks = element_line(size = 1.0, colour = "black"),
#        axis.ticks.length = unit(0.25, "cm"),
#        axis.text.x = element_text(size = 22, colour = "black"),
#        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
#        axis.title = element_text(size = 30, colour = log2ChIPColours[x]),
#        legend.position = "none",
#        panel.grid = element_blank(),
#        panel.border = element_rect(size = 3.5, colour = "black"),
#        panel.background = element_blank(),
#        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
#        plot.title = element_text(hjust = 0.5, size = 30)) +
#  ggtitle(bquote(.(ranLocAcc1NamePlot) ~ "(" * italic("n") ~ "=" ~
#                 .(prettyNum(summaryDFfeature$n[1],
#                             big.mark = ",", trim = T)) *
#                 ")"))
#}, mc.cores = length(log2ChIPNamesPlot))
#
#ggObjGA_combined <- grid.arrange(grobs = c(
#                                           ggObj1_combined_log2ChIP,
#                                           ggObj2_combined_log2ChIP,
#                                           ggObj3_combined_log2ChIP
#                                          ),
#                                 layout_matrix = cbind(
#                                                       1:length(c(log2ChIPNamesPlot)),
#                                                       (length(c(log2ChIPNamesPlot))+1):(length(c(log2ChIPNamesPlot))*2),
#                                                       ((length(c(log2ChIPNamesPlot))*2)+1):(length(c(log2ChIPNamesPlot))*3)
#                                                      ))
#ggsave(paste0(plotDir,
#              "log2ChIPcontrol_", align, "_avgProfiles_around_", quantiles, "quantiles",
#               "_by_", varType, "_in_", orderRegion,
#               "_of_Araport11_genes.pdf"),
#       plot = ggObjGA_combined,
#       height = 6.5*length(c(log2ChIPNamesPlot)), width = 21, limitsize = FALSE)
#
##### Free up memory by removing no longer required objects
#rm(
#   log2ChIP_featureAcc1Mats, log2ChIP_ranLocAcc1Mats,
#   log2ChIP_mats_quantiles,
#   wideDFfeature_list_log2ChIP,
#   tidyDFfeature_list_log2ChIP,
#   summaryDFfeature_list_log2ChIP,
#   summaryDFfeature_log2ChIP
#  ) 
#gc()
######
#
#
## control
## Add column names
#for(x in seq_along(control_featureAcc1Mats)) {
#  colnames(control_featureAcc1Mats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
#                                           paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
#                                           paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
#  colnames(control_ranLocAcc1Mats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
#                                          paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
#                                          paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
#}
#
## Subdivide coverage matrices into above-defined quantiles and random groupings
#control_mats_quantiles <- mclapply(seq_along(control_featureAcc1Mats), function(x) {
#  list(
#       # feature quantiles
#       lapply(1:quantiles, function(k) {
#         control_featureAcc1Mats[[x]][quantileIndicesAcc1[[k]],]
#       }),
#       # feature random groupings
#       lapply(1:quantiles, function(k) {
#         control_featureAcc1Mats[[x]][randomPCIndicesAcc1[[k]],]
#       }),
#       # random loci groupings
#       lapply(1:quantiles, function(k) {
#         control_ranLocAcc1Mats[[x]][quantileIndicesAcc1[[k]],]
#       })
#      ) 
#}, mc.cores = length(control_featureAcc1Mats))
#
## Transpose matrix and convert into dataframe
## in which first column is window name
#wideDFfeature_list_control <- mclapply(seq_along(control_mats_quantiles), function(x) {
#  lapply(seq_along(control_mats_quantiles[[x]]), function(y) {
#    lapply(seq_along(control_mats_quantiles[[x]][[y]]), function(k) {
#      data.frame(window = colnames(control_mats_quantiles[[x]][[y]][[k]]),
#                 t(control_mats_quantiles[[x]][[y]][[k]]))
#    })
#  })
#}, mc.cores = length(control_mats_quantiles))
#
## Convert into tidy data.frame (long format)
#tidyDFfeature_list_control  <- mclapply(seq_along(wideDFfeature_list_control), function(x) {
#  lapply(seq_along(control_mats_quantiles[[x]]), function(y) {
#    lapply(seq_along(control_mats_quantiles[[x]][[y]]), function(k) {
#      gather(data  = wideDFfeature_list_control[[x]][[y]][[k]],
#             key   = feature,
#             value = coverage,
#             -window)
#    })
#  }) 
#}, mc.cores = length(wideDFfeature_list_control))
#
## Order levels of factor "window" so that sequential levels
## correspond to sequential windows
#for(x in seq_along(tidyDFfeature_list_control)) {
#  for(y in seq_along(control_mats_quantiles[[x]])) {
#    for(k in seq_along(control_mats_quantiles[[x]][[y]])) {
#      tidyDFfeature_list_control[[x]][[y]][[k]]$window <- factor(tidyDFfeature_list_control[[x]][[y]][[k]]$window,
#                                                                  levels = as.character(wideDFfeature_list_control[[x]][[y]][[k]]$window))
#    }
#  }
#}
#
## Create summary data.frame in which each row corresponds to a window (Column 1),
## Column2 is the number of coverage values (features) per window,
## Column3 is the mean of coverage values per window,
## Column4 is the standard deviation of coverage values per window,
## Column5 is the standard error of the mean of coverage values per window,
## Column6 is the lower bound of the 95% confidence interval, and
## Column7 is the upper bound of the 95% confidence interval
#summaryDFfeature_list_control  <- mclapply(seq_along(tidyDFfeature_list_control), function(x) {
#  lapply(seq_along(control_mats_quantiles[[x]]), function(y) {
#    lapply(seq_along(control_mats_quantiles[[x]][[y]]), function(k) {
#      data.frame(window = as.character(wideDFfeature_list_control[[x]][[y]][[k]]$window),
#                 n      = tapply(X     = tidyDFfeature_list_control[[x]][[y]][[k]]$coverage,
#                                 INDEX = tidyDFfeature_list_control[[x]][[y]][[k]]$window,
#                                 FUN   = length),
#                 mean   = tapply(X     = tidyDFfeature_list_control[[x]][[y]][[k]]$coverage,
#                                 INDEX = tidyDFfeature_list_control[[x]][[y]][[k]]$window,
#                                 FUN   = mean,
#                                 na.rm = TRUE),
#                 sd     = tapply(X     = tidyDFfeature_list_control[[x]][[y]][[k]]$coverage,
#                                 INDEX = tidyDFfeature_list_control[[x]][[y]][[k]]$window,
#                                 FUN   = sd,
#                                 na.rm = TRUE))
#    })
#  })
#}, mc.cores = length(tidyDFfeature_list_control))
#
#for(x in seq_along(summaryDFfeature_list_control)) {
#  for(y in seq_along(control_mats_quantiles[[x]])) {
#    for(k in seq_along(control_mats_quantiles[[x]][[y]])) {
#      summaryDFfeature_list_control[[x]][[y]][[k]]$window <- factor(summaryDFfeature_list_control[[x]][[y]][[k]]$window,
#                                                                     levels = as.character(wideDFfeature_list_control[[x]][[y]][[k]]$window))
#      summaryDFfeature_list_control[[x]][[y]][[k]]$winNo <- factor(1:dim(summaryDFfeature_list_control[[x]][[y]][[k]])[1])
#      summaryDFfeature_list_control[[x]][[y]][[k]]$sem <- summaryDFfeature_list_control[[x]][[y]][[k]]$sd/sqrt(summaryDFfeature_list_control[[x]][[y]][[k]]$n-1)
#      summaryDFfeature_list_control[[x]][[y]][[k]]$CI_lower <- summaryDFfeature_list_control[[x]][[y]][[k]]$mean -
#        qt(0.975, df = summaryDFfeature_list_control[[x]][[y]][[k]]$n-1)*summaryDFfeature_list_control[[x]][[y]][[k]]$sem
#      summaryDFfeature_list_control[[x]][[y]][[k]]$CI_upper <- summaryDFfeature_list_control[[x]][[y]][[k]]$mean +
#        qt(0.975, df = summaryDFfeature_list_control[[x]][[y]][[k]]$n-1)*summaryDFfeature_list_control[[x]][[y]][[k]]$sem
#    }
#  }
#}
#
#quantileNames <- paste0(rep("Quantile ", quantiles), 1:quantiles)
#randomPCNames <- paste0(rep("Random ", quantiles), 1:quantiles)
#for(x in seq_along(summaryDFfeature_list_control)) {
#  # feature quantiles
#  names(summaryDFfeature_list_control[[x]][[1]]) <- quantileNames
#  # feature random groupings
#  names(summaryDFfeature_list_control[[x]][[2]]) <- randomPCNames
#  # random loci groupings
#  names(summaryDFfeature_list_control[[x]][[3]]) <- randomPCNames
#}
#
## Convert list of lists of lists of feature quantiles summaryDFfeature_list_control into
## a list of lists of single data.frames containing all feature quantiles for plotting
#summaryDFfeature_control  <- mclapply(seq_along(summaryDFfeature_list_control), function(x) {
#  lapply(seq_along(control_mats_quantiles[[x]]), function(y) {
#    bind_rows(summaryDFfeature_list_control[[x]][[y]], .id = "quantile")
#  })
#}, mc.cores = length(summaryDFfeature_list_control))
#for(x in seq_along(summaryDFfeature_control)) {
#  # feature quantiles
#  summaryDFfeature_control[[x]][[1]]$quantile <- factor(summaryDFfeature_control[[x]][[1]]$quantile,
#                                                         levels = names(summaryDFfeature_list_control[[x]][[1]]))
#  # feature random groupings
#  summaryDFfeature_control[[x]][[2]]$quantile <- factor(summaryDFfeature_control[[x]][[2]]$quantile,
#                                                         levels = names(summaryDFfeature_list_control[[x]][[2]]))
#  # random loci groupings
#  summaryDFfeature_control[[x]][[3]]$quantile <- factor(summaryDFfeature_control[[x]][[3]]$quantile,
#                                                         levels = names(summaryDFfeature_list_control[[x]][[3]]))
#}
#
## Define y-axis limits
#ymin_list_control <- lapply(seq_along(summaryDFfeature_control), function(x) {
#  min(c(summaryDFfeature_control[[x]][[1]]$CI_lower,
#        summaryDFfeature_control[[x]][[2]]$CI_lower,
#        summaryDFfeature_control[[x]][[3]]$CI_lower))
#})
#ymax_list_control <- lapply(seq_along(summaryDFfeature_control), function(x) {
#  max(c(summaryDFfeature_control[[x]][[1]]$CI_upper,
#        summaryDFfeature_control[[x]][[2]]$CI_upper,
#        summaryDFfeature_control[[x]][[3]]$CI_upper))
#})
#
## Define legend labels
#legendLabs_featureAcc1 <- lapply(seq_along(quantileNames), function(x) {
#  grobTree(textGrob(bquote(.(quantileNames[x])),
#                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
#                    gp = gpar(col = quantileColours[x], fontsize = 18)))
#})
#legendLabs_ranFeatAcc1 <- lapply(seq_along(randomPCNames), function(x) {
#  grobTree(textGrob(bquote(.(randomPCNames[x])),
#                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
#                    gp = gpar(col = quantileColours[x], fontsize = 18)))
#})
#legendLabs_ranLocAcc1 <- lapply(seq_along(randomPCNames), function(x) {
#  grobTree(textGrob(bquote(.(randomPCNames[x])),
#                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
#                    gp = gpar(col = quantileColours[x], fontsize = 18)))
#})
#
## Plot average profiles with 95% CI ribbon
### feature
#ggObj1_combined_control <- mclapply(seq_along(controlNamesPlot), function(x) {
#  summaryDFfeature <- summaryDFfeature_control[[x]][[1]]
#  ggplot(data = summaryDFfeature,
#         mapping = aes(x = winNo,
#                       y = mean,
#                       group = quantile)
#        ) +
#  geom_line(data = summaryDFfeature,
#            mapping = aes(colour = quantile),
#            size = 1) +
#  scale_colour_manual(values = quantileColours) +
#  geom_ribbon(data = summaryDFfeature,
#              mapping = aes(ymin = CI_lower,
#                            ymax = CI_upper,
#                            fill = quantile),
#              alpha = 0.4) +
#  scale_fill_manual(values = quantileColours) +
#  scale_y_continuous(limits = c(ymin_list_control[[x]], ymax_list_control[[x]]),
#                     labels = function(x) sprintf("%6.3f", x)) +
#  scale_x_discrete(breaks = c(1,
#                              (upstream/binSize)+1,
#                              (dim(summaryDFfeature_control[[x]][[1]])[1]/quantiles)-(downstream/binSize),
#                              dim(summaryDFfeature_control[[x]][[1]])[1]/quantiles),
#                   labels = c(paste0("-", flankName),
#                              featureStartLab,
#                              featureEndLab,
#                              paste0("+", flankName))) +
#  geom_vline(xintercept = c((upstream/binSize)+1,
#                            (dim(summaryDFfeature_control[[x]][[1]])[1]/quantiles)-(downstream/binSize)),
#             linetype = "dashed",
#             size = 1) +
#  labs(x = "",
#       y = controlNamesPlot[x]) +
#  annotation_custom(legendLabs_featureAcc1[[1]]) +
#  annotation_custom(legendLabs_featureAcc1[[2]]) +
#  annotation_custom(legendLabs_featureAcc1[[3]]) +
#  annotation_custom(legendLabs_featureAcc1[[4]]) +
#  theme_bw() +
#  theme(
#        axis.ticks = element_line(size = 1.0, colour = "black"),
#        axis.ticks.length = unit(0.25, "cm"),
#        axis.text.x = element_text(size = 22, colour = "black"),
#        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
#        axis.title = element_text(size = 30, colour = controlColours[x]),
#        legend.position = "none",
#        panel.grid = element_blank(),
#        panel.border = element_rect(size = 3.5, colour = "black"),
#        panel.background = element_blank(),
#        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
#        plot.title = element_text(hjust = 0.5, size = 30)) +
#  ggtitle(bquote(.(featureAcc1NamePlot) ~ "(" * italic("n") ~ "=" ~
#                 .(prettyNum(summaryDFfeature$n[1],
#                             big.mark = ",", trim = T)) *
#                 ")"))
#}, mc.cores = length(controlNamesPlot))
#
### ranFeat
#ggObj2_combined_control <- mclapply(seq_along(controlNamesPlot), function(x) {
#  summaryDFfeature <- summaryDFfeature_control[[x]][[2]]
#  ggplot(data = summaryDFfeature,
#         mapping = aes(x = winNo,
#                       y = mean,
#                       group = quantile)
#        ) +
#  geom_line(data = summaryDFfeature,
#            mapping = aes(colour = quantile),
#            size = 1) +
#  scale_colour_manual(values = quantileColours) +
#  geom_ribbon(data = summaryDFfeature,
#              mapping = aes(ymin = CI_lower,
#                            ymax = CI_upper,
#                            fill = quantile),
#              alpha = 0.4) +
#  scale_fill_manual(values = quantileColours) +
#  scale_y_continuous(limits = c(ymin_list_control[[x]], ymax_list_control[[x]]),
#                     labels = function(x) sprintf("%6.3f", x)) +
#  scale_x_discrete(breaks = c(1,
#                              (upstream/binSize)+1,
#                              (dim(summaryDFfeature_control[[x]][[2]])[1]/quantiles)-(downstream/binSize),
#                              dim(summaryDFfeature_control[[x]][[2]])[1]/quantiles),
#                   labels = c(paste0("-", flankName),
#                              featureStartLab,
#                              featureEndLab,
#                              paste0("+", flankName))) +
#  geom_vline(xintercept = c((upstream/binSize)+1,
#                            (dim(summaryDFfeature_control[[x]][[2]])[1]/quantiles)-(downstream/binSize)),
#             linetype = "dashed",
#             size = 1) +
#  labs(x = "",
#       y = controlNamesPlot[x]) +
#  annotation_custom(legendLabs_ranFeatAcc1[[1]]) +
#  annotation_custom(legendLabs_ranFeatAcc1[[2]]) +
#  annotation_custom(legendLabs_ranFeatAcc1[[3]]) +
#  annotation_custom(legendLabs_ranFeatAcc1[[4]]) +
#  theme_bw() +
#  theme(
#        axis.ticks = element_line(size = 1.0, colour = "black"),
#        axis.ticks.length = unit(0.25, "cm"),
#        axis.text.x = element_text(size = 22, colour = "black"),
#        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
#        axis.title = element_text(size = 30, colour = controlColours[x]),
#        legend.position = "none",
#        panel.grid = element_blank(),
#        panel.border = element_rect(size = 3.5, colour = "black"),
#        panel.background = element_blank(),
#        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
#        plot.title = element_text(hjust = 0.5, size = 30)) +
#  ggtitle(bquote(.(ranFeatAcc1NamePlot) ~ "(" * italic("n") ~ "=" ~
#                 .(prettyNum(summaryDFfeature$n[1],
#                             big.mark = ",", trim = T)) *
#                 ")"))
#}, mc.cores = length(controlNamesPlot))
#
### ranLoc
#ggObj3_combined_control <- mclapply(seq_along(controlNamesPlot), function(x) {
#  summaryDFfeature <- summaryDFfeature_control[[x]][[3]]
#  ggplot(data = summaryDFfeature,
#         mapping = aes(x = winNo,
#                       y = mean,
#                       group = quantile)
#        ) +
#  geom_line(data = summaryDFfeature,
#            mapping = aes(colour = quantile),
#            size = 1) +
#  scale_colour_manual(values = quantileColours) +
#  geom_ribbon(data = summaryDFfeature,
#              mapping = aes(ymin = CI_lower,
#                            ymax = CI_upper,
#                            fill = quantile),
#              alpha = 0.4) +
#  scale_fill_manual(values = quantileColours) +
#  scale_y_continuous(limits = c(ymin_list_control[[x]], ymax_list_control[[x]]),
#                     labels = function(x) sprintf("%6.3f", x)) +
#  scale_x_discrete(breaks = c(1,
#                              (upstream/binSize)+1,
#                              (dim(summaryDFfeature_control[[x]][[3]])[1]/quantiles)-(downstream/binSize),
#                              dim(summaryDFfeature_control[[x]][[3]])[1]/quantiles),
#                   labels = c(paste0("-", flankName),
#                              "Start",
#                              "End",
#                              paste0("+", flankName))) +
#  geom_vline(xintercept = c((upstream/binSize)+1,
#                            (dim(summaryDFfeature_control[[x]][[3]])[1]/quantiles)-(downstream/binSize)),
#             linetype = "dashed",
#             size = 1) +
#  labs(x = "",
#       y = controlNamesPlot[x]) +
#  annotation_custom(legendLabs_ranLocAcc1[[1]]) +
#  annotation_custom(legendLabs_ranLocAcc1[[2]]) +
#  annotation_custom(legendLabs_ranLocAcc1[[3]]) +
#  annotation_custom(legendLabs_ranLocAcc1[[4]]) +
#  theme_bw() +
#  theme(
#        axis.ticks = element_line(size = 1.0, colour = "black"),
#        axis.ticks.length = unit(0.25, "cm"),
#        axis.text.x = element_text(size = 22, colour = "black"),
#        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
#        axis.title = element_text(size = 30, colour = controlColours[x]),
#        legend.position = "none",
#        panel.grid = element_blank(),
#        panel.border = element_rect(size = 3.5, colour = "black"),
#        panel.background = element_blank(),
#        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
#        plot.title = element_text(hjust = 0.5, size = 30)) +
#  ggtitle(bquote(.(ranLocAcc1NamePlot) ~ "(" * italic("n") ~ "=" ~
#                 .(prettyNum(summaryDFfeature$n[1],
#                             big.mark = ",", trim = T)) *
#                 ")"))
#}, mc.cores = length(controlNamesPlot))
#
#ggObjGA_combined <- grid.arrange(grobs = c(
#                                           ggObj1_combined_control,
#                                           ggObj2_combined_control,
#                                           ggObj3_combined_control
#                                          ),
#                                 layout_matrix = cbind(
#                                                       1:length(c(controlNamesPlot)),
#                                                       (length(c(controlNamesPlot))+1):(length(c(controlNamesPlot))*2),
#                                                       ((length(c(controlNamesPlot))*2)+1):(length(c(controlNamesPlot))*3)
#                                                      ))
#ggsave(paste0(plotDir,
#              "control_", align, "_avgProfiles_around_", quantiles, "quantiles",
#               "_by_", varType, "_in_", orderRegion,
#               "_of_Araport11_genes.pdf"),
#       plot = ggObjGA_combined,
#       height = 6.5*length(c(controlNamesPlot)), width = 21, limitsize = FALSE)
#
##### Free up memory by removing no longer required objects
#rm(
#   control_featureAcc1Mats, control_ranLocAcc1Mats,
#   control_mats_quantiles,
#   wideDFfeature_list_control,
#   tidyDFfeature_list_control,
#   summaryDFfeature_list_control,
#   summaryDFfeature_control
#  ) 
#gc()
######


# ChIP
# Add column names
for(x in seq_along(ChIP_featureAcc1Mats)) {
  colnames(ChIP_featureAcc1Mats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                           paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                           paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
  colnames(ChIP_ranLocAcc1Mats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                          paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                          paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
  colnames(ChIP_featureAcc2Mats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                           paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                           paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
  colnames(ChIP_ranLocAcc2Mats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                          paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                          paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
}

# Subdivide coverage matrices into above-defined quantiles and random groupings
ChIP_mats_quantiles <- mclapply(seq_along(ChIP_featureAcc1Mats), function(x) {
  list(
       # featureAcc1 quantiles
       lapply(1:quantiles, function(k) {
         ChIP_featureAcc1Mats[[x]][quantileIndicesAcc1[[k]],]
       }),
       # featureAcc1 random groupings
       lapply(1:quantiles, function(k) {
         ChIP_featureAcc1Mats[[x]][randomPCIndicesAcc1[[k]],]
       }),
       # ranLocAcc1 groupings
       lapply(1:quantiles, function(k) {
         ChIP_ranLocAcc1Mats[[x]][quantileIndicesAcc1[[k]],]
       }),
       # featureAcc2 quantiles
       lapply(1:quantiles, function(k) {
         ChIP_featureAcc2Mats[[x]][quantileIndicesAcc1[[k]],]
       }),
       # featureAcc2 random groupings
       lapply(1:quantiles, function(k) {
         ChIP_featureAcc2Mats[[x]][randomPCIndicesAcc1[[k]],]
       }),
       # ranLocAcc2 groupings
       lapply(1:quantiles, function(k) {
         ChIP_ranLocAcc2Mats[[x]][quantileIndicesAcc1[[k]],]
       })
      ) 
}, mc.cores = length(ChIP_featureAcc1Mats))

# Transpose matrix and convert into dataframe
# in which first column is window name
wideDFfeature_list_ChIP <- mclapply(seq_along(ChIP_mats_quantiles), function(x) {
  lapply(seq_along(ChIP_mats_quantiles[[x]]), function(y) {
    lapply(seq_along(ChIP_mats_quantiles[[x]][[y]]), function(k) {
      data.frame(window = colnames(ChIP_mats_quantiles[[x]][[y]][[k]]),
                 t(ChIP_mats_quantiles[[x]][[y]][[k]]))
    })
  })
}, mc.cores = length(ChIP_mats_quantiles))

# Convert into tidy data.frame (long format)
tidyDFfeature_list_ChIP  <- mclapply(seq_along(wideDFfeature_list_ChIP), function(x) {
  lapply(seq_along(ChIP_mats_quantiles[[x]]), function(y) {
    lapply(seq_along(ChIP_mats_quantiles[[x]][[y]]), function(k) {
      gather(data  = wideDFfeature_list_ChIP[[x]][[y]][[k]],
             key   = feature,
             value = coverage,
             -window)
    })
  }) 
}, mc.cores = length(wideDFfeature_list_ChIP))

# Order levels of factor "window" so that sequential levels
# correspond to sequential windows
for(x in seq_along(tidyDFfeature_list_ChIP)) {
  for(y in seq_along(ChIP_mats_quantiles[[x]])) {
    for(k in seq_along(ChIP_mats_quantiles[[x]][[y]])) {
      tidyDFfeature_list_ChIP[[x]][[y]][[k]]$window <- factor(tidyDFfeature_list_ChIP[[x]][[y]][[k]]$window,
                                                                  levels = as.character(wideDFfeature_list_ChIP[[x]][[y]][[k]]$window))
    }
  }
}

# Create summary data.frame in which each row corresponds to a window (Column 1),
# Column2 is the number of coverage values (features) per window,
# Column3 is the mean of coverage values per window,
# Column4 is the standard deviation of coverage values per window,
# Column5 is the standard error of the mean of coverage values per window,
# Column6 is the lower bound of the 95% confidence interval, and
# Column7 is the upper bound of the 95% confidence interval
summaryDFfeature_list_ChIP  <- mclapply(seq_along(tidyDFfeature_list_ChIP), function(x) {
  lapply(seq_along(ChIP_mats_quantiles[[x]]), function(y) {
    lapply(seq_along(ChIP_mats_quantiles[[x]][[y]]), function(k) {
      data.frame(window = as.character(wideDFfeature_list_ChIP[[x]][[y]][[k]]$window),
                 n      = tapply(X     = tidyDFfeature_list_ChIP[[x]][[y]][[k]]$coverage,
                                 INDEX = tidyDFfeature_list_ChIP[[x]][[y]][[k]]$window,
                                 FUN   = length),
                 mean   = tapply(X     = tidyDFfeature_list_ChIP[[x]][[y]][[k]]$coverage,
                                 INDEX = tidyDFfeature_list_ChIP[[x]][[y]][[k]]$window,
                                 FUN   = mean,
                                 na.rm = TRUE),
                 sd     = tapply(X     = tidyDFfeature_list_ChIP[[x]][[y]][[k]]$coverage,
                                 INDEX = tidyDFfeature_list_ChIP[[x]][[y]][[k]]$window,
                                 FUN   = sd,
                                 na.rm = TRUE))
    })
  })
}, mc.cores = length(tidyDFfeature_list_ChIP))

for(x in seq_along(summaryDFfeature_list_ChIP)) {
  for(y in seq_along(ChIP_mats_quantiles[[x]])) {
    for(k in seq_along(ChIP_mats_quantiles[[x]][[y]])) {
      summaryDFfeature_list_ChIP[[x]][[y]][[k]]$window <- factor(summaryDFfeature_list_ChIP[[x]][[y]][[k]]$window,
                                                                     levels = as.character(wideDFfeature_list_ChIP[[x]][[y]][[k]]$window))
      summaryDFfeature_list_ChIP[[x]][[y]][[k]]$winNo <- factor(1:dim(summaryDFfeature_list_ChIP[[x]][[y]][[k]])[1])
      summaryDFfeature_list_ChIP[[x]][[y]][[k]]$sem <- summaryDFfeature_list_ChIP[[x]][[y]][[k]]$sd/sqrt(summaryDFfeature_list_ChIP[[x]][[y]][[k]]$n-1)
      summaryDFfeature_list_ChIP[[x]][[y]][[k]]$CI_lower <- summaryDFfeature_list_ChIP[[x]][[y]][[k]]$mean -
        qt(0.975, df = summaryDFfeature_list_ChIP[[x]][[y]][[k]]$n-1)*summaryDFfeature_list_ChIP[[x]][[y]][[k]]$sem
      summaryDFfeature_list_ChIP[[x]][[y]][[k]]$CI_upper <- summaryDFfeature_list_ChIP[[x]][[y]][[k]]$mean +
        qt(0.975, df = summaryDFfeature_list_ChIP[[x]][[y]][[k]]$n-1)*summaryDFfeature_list_ChIP[[x]][[y]][[k]]$sem
    }
  }
}

quantileNames <- paste0(rep("Quantile ", quantiles), 1:quantiles)
randomPCNames <- paste0(rep("Random ", quantiles), 1:quantiles)
for(x in seq_along(summaryDFfeature_list_ChIP)) {
  # featureAcc1 quantiles
  names(summaryDFfeature_list_ChIP[[x]][[1]]) <- quantileNames
  # featureAcc1 random groupings
  names(summaryDFfeature_list_ChIP[[x]][[2]]) <- randomPCNames
  # ranLocAcc1 groupings
  names(summaryDFfeature_list_ChIP[[x]][[3]]) <- randomPCNames
  # featureAcc2 quantiles
  names(summaryDFfeature_list_ChIP[[x]][[4]]) <- quantileNames
  # featureAcc2 random groupings
  names(summaryDFfeature_list_ChIP[[x]][[5]]) <- randomPCNames
  # ranLocAcc2 groupings
  names(summaryDFfeature_list_ChIP[[x]][[6]]) <- randomPCNames
}

# Convert list of lists of lists of feature quantiles summaryDFfeature_list_ChIP into
# a list of lists of single data.frames containing all feature quantiles for plotting
summaryDFfeature_ChIP  <- mclapply(seq_along(summaryDFfeature_list_ChIP), function(x) {
  lapply(seq_along(ChIP_mats_quantiles[[x]]), function(y) {
    bind_rows(summaryDFfeature_list_ChIP[[x]][[y]], .id = "quantile")
  })
}, mc.cores = length(summaryDFfeature_list_ChIP))
for(x in seq_along(summaryDFfeature_ChIP)) {
  # featureAcc1 quantiles
  summaryDFfeature_ChIP[[x]][[1]]$quantile <- factor(summaryDFfeature_ChIP[[x]][[1]]$quantile,
                                                     levels = names(summaryDFfeature_list_ChIP[[x]][[1]]))
  # featureAcc1 random groupings
  summaryDFfeature_ChIP[[x]][[2]]$quantile <- factor(summaryDFfeature_ChIP[[x]][[2]]$quantile,
                                                     levels = names(summaryDFfeature_list_ChIP[[x]][[2]]))
  # ranLocAcc1 groupings
  summaryDFfeature_ChIP[[x]][[3]]$quantile <- factor(summaryDFfeature_ChIP[[x]][[3]]$quantile,
                                                     levels = names(summaryDFfeature_list_ChIP[[x]][[3]]))
  # featureAcc2 quantiles
  summaryDFfeature_ChIP[[x]][[4]]$quantile <- factor(summaryDFfeature_ChIP[[x]][[4]]$quantile,
                                                     levels = names(summaryDFfeature_list_ChIP[[x]][[4]]))
  # featureAcc2 random groupings
  summaryDFfeature_ChIP[[x]][[5]]$quantile <- factor(summaryDFfeature_ChIP[[x]][[5]]$quantile,
                                                     levels = names(summaryDFfeature_list_ChIP[[x]][[5]]))
  # ranLocAcc2 groupings
  summaryDFfeature_ChIP[[x]][[6]]$quantile <- factor(summaryDFfeature_ChIP[[x]][[6]]$quantile,
                                                     levels = names(summaryDFfeature_list_ChIP[[x]][[6]]))
}

# Define y-axis limits
ymin_list_ChIP <- lapply(seq_along(summaryDFfeature_ChIP), function(x) {
  min(c(summaryDFfeature_ChIP[[x]][[1]]$CI_lower,
        summaryDFfeature_ChIP[[x]][[2]]$CI_lower,
        summaryDFfeature_ChIP[[x]][[3]]$CI_lower,
        summaryDFfeature_ChIP[[x]][[4]]$CI_lower,
        summaryDFfeature_ChIP[[x]][[5]]$CI_lower,
        summaryDFfeature_ChIP[[x]][[6]]$CI_lower))
})
ymax_list_ChIP <- lapply(seq_along(summaryDFfeature_ChIP), function(x) {
  max(c(summaryDFfeature_ChIP[[x]][[1]]$CI_upper,
        summaryDFfeature_ChIP[[x]][[2]]$CI_upper,
        summaryDFfeature_ChIP[[x]][[3]]$CI_upper,
        summaryDFfeature_ChIP[[x]][[4]]$CI_upper,
        summaryDFfeature_ChIP[[x]][[5]]$CI_upper,
        summaryDFfeature_ChIP[[x]][[6]]$CI_upper))
})

# Repeated ymin and ymax values for as many ChIPNames (and defined as list)
# for consistency with above definitions and convenience
ymin_list_ChIP <- as.list(rep(min(unlist(ymin_list_ChIP)), length(ymin_list_ChIP)))
ymax_list_ChIP <- as.list(rep(max(unlist(ymax_list_ChIP)), length(ymax_list_ChIP)))

# Define legend labels
legendLabs_featureAcc1 <- lapply(seq_along(quantileNames), function(x) {
  grobTree(textGrob(bquote(.(quantileNames[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = quantileColours[x], fontsize = 18)))
})
legendLabs_ranFeatAcc1 <- lapply(seq_along(randomPCNames), function(x) {
  grobTree(textGrob(bquote(.(randomPCNames[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = quantileColours[x], fontsize = 18)))
})
legendLabs_ranLocAcc1 <- lapply(seq_along(randomPCNames), function(x) {
  grobTree(textGrob(bquote(.(randomPCNames[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = quantileColours[x], fontsize = 18)))
})
legendLabs_featureAcc2 <- lapply(seq_along(quantileNames), function(x) {
  grobTree(textGrob(bquote(.(quantileNames[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = quantileColours[x], fontsize = 18)))
})
legendLabs_ranFeatAcc2 <- lapply(seq_along(randomPCNames), function(x) {
  grobTree(textGrob(bquote(.(randomPCNames[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = quantileColours[x], fontsize = 18)))
})
legendLabs_ranLocAcc2 <- lapply(seq_along(randomPCNames), function(x) {
  grobTree(textGrob(bquote(.(randomPCNames[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = quantileColours[x], fontsize = 18)))
})

# Plot average profiles with 95% CI ribbon
## feature
ggObj1_combined_ChIP <- mclapply(seq_along(ChIPNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_ChIP[[x]][[1]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_ChIP[[x]], ymax_list_ChIP[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_ChIP[[x]][[1]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_ChIP[[x]][[1]])[1]/quantiles),
                   labels = c(paste0("-", flankName),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankName))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_ChIP[[x]][[1]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = ChIPNamesPlot[x]) +
  annotation_custom(legendLabs_featureAcc1[[1]]) +
  annotation_custom(legendLabs_featureAcc1[[2]]) +
  annotation_custom(legendLabs_featureAcc1[[3]]) +
  annotation_custom(legendLabs_featureAcc1[[4]]) +
  annotation_custom(legendLabs_featureAcc1[[5]]) +
  annotation_custom(legendLabs_featureAcc1[[6]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = ChIPColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 0.5, size = 30)) +
  ggtitle(bquote(.(featureAcc1NamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(ChIPNamesPlot))

## ranFeat
ggObj2_combined_ChIP <- mclapply(seq_along(ChIPNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_ChIP[[x]][[2]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_ChIP[[x]], ymax_list_ChIP[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_ChIP[[x]][[2]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_ChIP[[x]][[2]])[1]/quantiles),
                   labels = c(paste0("-", flankName),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankName))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_ChIP[[x]][[2]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = ChIPNamesPlot[x]) +
  annotation_custom(legendLabs_ranFeatAcc1[[1]]) +
  annotation_custom(legendLabs_ranFeatAcc1[[2]]) +
  annotation_custom(legendLabs_ranFeatAcc1[[3]]) +
  annotation_custom(legendLabs_ranFeatAcc1[[4]]) +
  annotation_custom(legendLabs_ranFeatAcc1[[5]]) +
  annotation_custom(legendLabs_ranFeatAcc1[[6]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = ChIPColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 0.5, size = 30)) +
  ggtitle(bquote(.(ranFeatAcc1NamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(ChIPNamesPlot))

## ranLoc
ggObj3_combined_ChIP <- mclapply(seq_along(ChIPNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_ChIP[[x]][[3]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_ChIP[[x]], ymax_list_ChIP[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_ChIP[[x]][[3]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_ChIP[[x]][[3]])[1]/quantiles),
                   labels = c(paste0("-", flankName),
                              "Start",
                              "End",
                              paste0("+", flankName))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_ChIP[[x]][[3]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = ChIPNamesPlot[x]) +
  annotation_custom(legendLabs_ranLocAcc1[[1]]) +
  annotation_custom(legendLabs_ranLocAcc1[[2]]) +
  annotation_custom(legendLabs_ranLocAcc1[[3]]) +
  annotation_custom(legendLabs_ranLocAcc1[[4]]) +
  annotation_custom(legendLabs_ranLocAcc1[[5]]) +
  annotation_custom(legendLabs_ranLocAcc1[[6]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = ChIPColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 0.5, size = 30)) +
  ggtitle(bquote(.(ranLocAcc1NamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(ChIPNamesPlot))

## feature
ggObj4_combined_ChIP <- mclapply(seq_along(ChIPNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_ChIP[[x]][[4]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_ChIP[[x]], ymax_list_ChIP[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_ChIP[[x]][[4]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_ChIP[[x]][[4]])[1]/quantiles),
                   labels = c(paste0("-", flankName),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankName))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_ChIP[[x]][[4]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = ChIPNamesPlot[x]) +
  annotation_custom(legendLabs_featureAcc2[[1]]) +
  annotation_custom(legendLabs_featureAcc2[[2]]) +
  annotation_custom(legendLabs_featureAcc2[[3]]) +
  annotation_custom(legendLabs_featureAcc2[[4]]) +
  annotation_custom(legendLabs_featureAcc2[[5]]) +
  annotation_custom(legendLabs_featureAcc2[[6]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = ChIPColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 0.5, size = 30)) +
  ggtitle(bquote(.(featureAcc2NamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(ChIPNamesPlot))

## ranFeat
ggObj5_combined_ChIP <- mclapply(seq_along(ChIPNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_ChIP[[x]][[5]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_ChIP[[x]], ymax_list_ChIP[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_ChIP[[x]][[5]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_ChIP[[x]][[5]])[1]/quantiles),
                   labels = c(paste0("-", flankName),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankName))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_ChIP[[x]][[5]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = ChIPNamesPlot[x]) +
  annotation_custom(legendLabs_ranFeatAcc2[[1]]) +
  annotation_custom(legendLabs_ranFeatAcc2[[2]]) +
  annotation_custom(legendLabs_ranFeatAcc2[[3]]) +
  annotation_custom(legendLabs_ranFeatAcc2[[4]]) +
  annotation_custom(legendLabs_ranFeatAcc2[[5]]) +
  annotation_custom(legendLabs_ranFeatAcc2[[6]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = ChIPColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 0.5, size = 30)) +
  ggtitle(bquote(.(ranFeatAcc2NamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(ChIPNamesPlot))

## ranLoc
ggObj6_combined_ChIP <- mclapply(seq_along(ChIPNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_ChIP[[x]][[6]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColours) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColours) +
  scale_y_continuous(limits = c(ymin_list_ChIP[[x]], ymax_list_ChIP[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_ChIP[[x]][[6]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_ChIP[[x]][[6]])[1]/quantiles),
                   labels = c(paste0("-", flankName),
                              "Start",
                              "End",
                              paste0("+", flankName))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_ChIP[[x]][[6]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = ChIPNamesPlot[x]) +
  annotation_custom(legendLabs_ranLocAcc2[[1]]) +
  annotation_custom(legendLabs_ranLocAcc2[[2]]) +
  annotation_custom(legendLabs_ranLocAcc2[[3]]) +
  annotation_custom(legendLabs_ranLocAcc2[[4]]) +
  annotation_custom(legendLabs_ranLocAcc2[[5]]) +
  annotation_custom(legendLabs_ranLocAcc2[[6]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = ChIPColours[x]),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
        plot.title = element_text(hjust = 0.5, size = 30)) +
  ggtitle(bquote(.(ranLocAcc2NamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
}, mc.cores = length(ChIPNamesPlot))

ggObjGA_combined <- grid.arrange(grobs = c(
                                           ggObj1_combined_ChIP,
                                           ggObj2_combined_ChIP,
                                           ggObj3_combined_ChIP,
                                           ggObj4_combined_ChIP,
                                           ggObj5_combined_ChIP,
                                           ggObj6_combined_ChIP
                                          ),
                                 layout_matrix = cbind(
                                                       1:length(c(ChIPNamesPlot)),
                                                       (length(c(ChIPNamesPlot))+1):(length(c(ChIPNamesPlot))*2),
                                                       ((length(c(ChIPNamesPlot))*2)+1):(length(c(ChIPNamesPlot))*3),
                                                       ((length(c(ChIPNamesPlot))*3)+1):(length(c(ChIPNamesPlot))*4),
                                                       ((length(c(ChIPNamesPlot))*4)+1):(length(c(ChIPNamesPlot))*5),
                                                       ((length(c(ChIPNamesPlot))*5)+1):(length(c(ChIPNamesPlot))*6)
                                                      ))
ggsave(paste0(plotDir,
              "ChIP_", align, "_avgProfiles_around_featuresAcc1_ortho_", quantiles, "quantiles",
              "_", genomeRegion, "_by_", varType, "_freq_in_", orderRegion,
              "_of_Acc1_Chr_genes_in_", refbase, "_",
              paste0(chrName, collapse = "_"), ".pdf"),
       plot = ggObjGA_combined,
       height = 6.5*length(c(ChIPNamesPlot)), width = 7*6, limitsize = FALSE)

#### Free up memory by removing no longer required objects
rm(
   ChIP_featureAcc1Mats, ChIP_ranLocAcc1Mats,
   ChIP_mats_quantiles,
   wideDFfeature_list_ChIP,
   tidyDFfeature_list_ChIP,
   summaryDFfeature_list_ChIP,
   summaryDFfeature_ChIP
  ) 
gc()
#####
