#!/applications/R/R-4.0.0/bin/Rscript

# author: Andy Tock
# contact: ajt200@cam.ac.uk
# date: 31.08.2021

# Calculate and plot metaprofiles of ChIP-seq
# (feature windowed means and 95% confidence intervals, CIs)
# for each group of features, defined either by
# decreasing absolute change in DNA methylation or randomly

# Usage:
# /applications/R/R-4.0.0/bin/Rscript TEsuperfam_DNAmethDiff_quantile_metaprofiles_absolute_change.R 'cmt3_BSseq_Rep1' 'WT_BSseq_Rep2_2013,WT_BSseq_Rep3_2013' 'kss_BSseq_Rep1' 'WT_BSseq_Rep2_2013,WT_BSseq_Rep3_2013' both 'Chr1,Chr2,Chr3,Chr4,Chr5' 2000 2000 2kb '2 kb' 10 10bp bodies CHG Gypsy_LTR 6 t2t-col.20210610 '0.02,0.96'

#libName1 <- unlist(strsplit("cmt3_BSseq_Rep1",
#                            split = ","))
#controlName1 <- unlist(strsplit("WT_BSseq_Rep2_2013,WT_BSseq_Rep3_2013",
#                                split = ","))
#libName2 <- unlist(strsplit("kss_BSseq_Rep1",
#                            split = ","))
#controlName2 <- unlist(strsplit("WT_BSseq_Rep2_2013,WT_BSseq_Rep3_2013",
#                                split = ","))
#align <- "both"
#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split = ","))
#bodyLength <- 2000
#upstream <- 2000
#downstream <- 2000
#flankName <- "2kb"
#flankNamePlot <- "2 kb"
#binSize <- 10
#binName <- "10bp"
#sortRegion <- "bodies"
#context <- "CHG"
#superfamName <- "Gypsy_LTR"
#quantiles <- 6
#refbase <- "t2t-col.20210610"
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
libName1 <- unlist(strsplit(args[1],
                            split = ","))
controlName1 <- unlist(strsplit(args[2],
                                split = ","))
libName2 <- unlist(strsplit(args[3],
                            split = ","))
controlName2 <- unlist(strsplit(args[4],
                                split = ","))
align <- args[5]
chrName <- unlist(strsplit(args[6],
                           split = ","))
bodyLength <- as.numeric(args[7])
upstream <- as.numeric(args[8])
downstream <- as.numeric(args[8])
flankName <- args[9]
flankNamePlot <- args[10]
binSize <- as.numeric(args[11])
binName <- args[12]
sortRegion <- args[13]
context <- args[14]
superfamName <- args[15]
quantiles <- as.numeric(args[16])
refbase <- args[17]
legendPos <- as.numeric(unlist(strsplit(args[18],
                                        split = ",")))

options(stringsAsFactors = F)
library(rtracklayer)
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
                 "/quantiles_by_", context, "_absolute_change_at_", superfamName, "_", sortRegion, "/")
plotDir <- paste0(outDir, "plots/")
system(paste0("[ -d ", outDir, " ] || mkdir -p ", outDir))
system(paste0("[ -d ", plotDir, " ] || mkdir -p ", plotDir))

# Define plot titles
featureAcc1NamePlot <- paste0(superfamName, " quantiles (", sub("_.+", "", libName1), " ", context, ")")
ranFeatAcc1NamePlot <- "Random quantiles (cmt3)"
ranLocAcc1NamePlot <- "RanLoc quantiles (cmt3)"
featureAcc2NamePlot <- paste0(superfamName, " quantiles (", sub("_.+", "", libName2), " ", context, ")")
ranFeatAcc2NamePlot <- "Random quantiles (kss)"
ranLocAcc2NamePlot <- "RanLoc quantiles (kss)"

# Define quantile colours
quantileColoursAcc1 <- rev(viridis(quantiles))
quantileColoursAcc1[1] <- "orange"
quantileColoursAcc2 <- rev(plasma(quantiles))
quantileColoursAcc2[1] <- "gold"

# Define feature start and end labels for plotting
featureStartLab <- "Start"
featureEndLab <- "End"

# Load table of featuresAcc1_ortho_DF grouped into quantiles
featuresAcc1_ortho_DF <- fread(
  paste0(paste0(chrName, collapse = "_"), "/",
         "feature_", quantiles, "quantiles",
         "_by_", context, "_change_in_",
         paste0(libName1, collapse = "_"), "_vs_", paste0(controlName1, collapse = "_"),
         "_in_", sortRegion, "_of_", superfamName, "_in_", refbase, "_", paste0(chrName, collapse = "_"), ".tsv"),
  header = T, sep = "\t")

# Load features to confirm feature (row) ordering in "featuresAcc1_ortho_DF" is the same
# as in "features" (which was used for generating the coverage matrices)
featuresAcc1_ortho_DF_bed <- read.table(
  paste0(paste0(chrName, collapse = "_"), "/",
         "feature_", quantiles, "quantiles",
         "_by_", context, "_change_in_",
         paste0(libName1, collapse = "_"), "_vs_", paste0(controlName1, collapse = "_"),
         "_in_", sortRegion, "_of_", superfamName, "_in_", refbase, "_", paste0(chrName, collapse = "_"), ".bed"),
  header = F, sep = "\t")
colnames(featuresAcc1_ortho_DF_bed) <- c("chr", "start", "end", "ID", "score", "strand")

stopifnot(identical(as.character(featuresAcc1_ortho_DF$start),
                    as.character(featuresAcc1_ortho_DF_bed$start+1)))
stopifnot(identical(as.character(featuresAcc1_ortho_DF$end),
                    as.character(featuresAcc1_ortho_DF_bed$end)))
rm(featuresAcc1_ortho_DF_bed); gc()

# Get row indices for each feature quantile
quantileIndicesAcc1 <- lapply(1:quantiles, function(k) {
  which(featuresAcc1_ortho_DF$ac_quantile == paste0("Quantile ", k))
})

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
chrs <- chrName[which(chrName %in% featuresAcc1_ortho_DF$seqid)]
randomPCIndicesAcc1 <- lapply(1:quantiles, function(k) {
  randomPCIndicesAcc1k <- NULL
  for(i in 1:length(chrs)) {
  # Define seed so that random selections are reproducible
  set.seed(93750174)
    randomPCfeatureskChr <- selectRandomFeatures(features = featuresAcc1_ortho_DF[featuresAcc1_ortho_DF$seqid == chrs[i],],
                                                 n = dim(featuresAcc1_ortho_DF[featuresAcc1_ortho_DF$ac_quantile == paste0("Quantile ", k) &
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


# Load table of featuresAcc2_ortho_DF grouped into quantiles
featuresAcc2_ortho_DF <- fread(
  paste0(paste0(chrName, collapse = "_"), "/",
         "feature_", quantiles, "quantiles",
         "_by_", context, "_change_in_",
         paste0(libName2, collapse = "_"), "_vs_", paste0(controlName2, collapse = "_"),
         "_in_", sortRegion, "_of_", superfamName, "_in_", refbase, "_", paste0(chrName, collapse = "_"), ".tsv"),
  header = T, sep = "\t")

# Load features to confirm feature (row) ordering in "featuresAcc2_ortho_DF" is the same
# as in "features" (which was used for generating the coverage matrices)
featuresAcc2_ortho_DF_bed <- read.table(
  paste0(paste0(chrName, collapse = "_"), "/",
         "feature_", quantiles, "quantiles",
         "_by_", context, "_change_in_",
         paste0(libName2, collapse = "_"), "_vs_", paste0(controlName2, collapse = "_"),
         "_in_", sortRegion, "_of_", superfamName, "_in_", refbase, "_", paste0(chrName, collapse = "_"), ".bed"),
  header = F, sep = "\t")
colnames(featuresAcc2_ortho_DF_bed) <- c("chr", "start", "end", "ID", "score", "strand")

stopifnot(identical(as.character(featuresAcc2_ortho_DF$start),
                    as.character(featuresAcc2_ortho_DF_bed$start+1)))
stopifnot(identical(as.character(featuresAcc2_ortho_DF$end),
                    as.character(featuresAcc2_ortho_DF_bed$end)))
rm(featuresAcc2_ortho_DF_bed); gc()

# Get row indices for each feature quantile
quantileIndicesAcc2 <- lapply(1:quantiles, function(k) {
  which(featuresAcc2_ortho_DF$ac_quantile == paste0("Quantile ", k))
})

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
chrs <- chrName[which(chrName %in% featuresAcc2_ortho_DF$seqid)]
randomPCIndicesAcc2 <- lapply(1:quantiles, function(k) {
  randomPCIndicesAcc2k <- NULL
  for(i in 1:length(chrs)) {
  # Define seed so that random selections are reproducible
  set.seed(93750174)
    randomPCfeatureskChr <- selectRandomFeatures(features = featuresAcc2_ortho_DF[featuresAcc2_ortho_DF$seqid == chrs[i],],
                                                 n = dim(featuresAcc2_ortho_DF[featuresAcc2_ortho_DF$ac_quantile == paste0("Quantile ", k) &
                                                                               featuresAcc2_ortho_DF$seqid == chrs[i],])[1])
    randomPCIndicesAcc2k <- c(randomPCIndicesAcc2k, as.integer(rownames(randomPCfeatureskChr)))
  }
  randomPCIndicesAcc2k
})
# Confirm per-chromosome feature numbers are the same for quantiles and random groupings
lapply(seq_along(1:quantiles), function(k) {
  sapply(seq_along(chrs), function(x) {
    if(!identical(dim(featuresAcc2_ortho_DF[randomPCIndicesAcc2[[k]],][featuresAcc2_ortho_DF[randomPCIndicesAcc2[[k]],]$chr == chrs[x],]),
                  dim(featuresAcc2_ortho_DF[quantileIndicesAcc2[[k]],][featuresAcc2_ortho_DF[quantileIndicesAcc2[[k]],]$chr == chrs[x],])))    {
      stop("Quantile features and random features do not consist of the same number of features per chromosome")
    }
  })
})



# Load feature matrices for each chromatin dataset, calculate log2(ChIP/control),
ChIPNames <- c(
               "cmt3_DMC1_V5_Rep1_ChIP",
               "kss_DMC1_V5_Rep1_ChIP",
               "Col_DMC1_V5_Rep1_ChIP",
               "Col_DMC1_V5_Rep2_ChIP",
               "kss_SPO11oligos_Rep1",
               "kss_SPO11oligos_Rep2",
               "met1_SPO11oligos_Rep1",
               "met1_SPO11oligos_Rep2",
               "met1_SPO11oligos_Rep3",
               "WT_SPO11oligos_Rep1",
               "WT_SPO11oligos_Rep2",
               "WT_SPO11oligos_Rep3"
              )
ChIPNamesDir <- c(
                  rep(paste0("20190917_dh580_Athaliana_ChIPseq_DMC1/fastq_pooled/snakemake_ChIPseq_", refbase), 4),
                  rep(paste0("160518_Kyuha_SPO11oligos/kss/snakemake_SPO11oligos_", refbase), 2),
                  rep(paste0("160518_Kyuha_SPO11oligos/met1/snakemake_SPO11oligos_", refbase), 3),
                  rep(paste0("160518_Kyuha_SPO11oligos/WT/snakemake_SPO11oligos_", refbase), 3)
                 )
log2ChIPNamesPlot <- c(
                       "(cmt3/WT) DMC1",
                       "(kss/WT) DMC1",
                       "(WT Rep1/input) DMC1",
                       "(WT Rep2/input) DMC1",
                       "(kss Rep1/WT) SPO11-1",
                       "(kss Rep2/WT) SPO11-1",
                       "(met1 Rep1/WT) SPO11-1",
                       "(met1 Rep2/WT) SPO11-1",
                       "(met1 Rep3/WT) SPO11-1",
                       "(WT Rep1/gDNA) SPO11-1",
                       "(WT Rep2/gDNA) SPO11-1",
                       "(WT Rep3/gDNA) SPO11-1"
                      )
ChIPNamesPlot <- c(
                   "cmt3 DMC1 Rep1",
                   "kss DMC1 Rep1",
                   "WT DMC1 Rep1",
                   "WT DMC1 Rep2",
                   "kss SPO11-1 Rep1",
                   "kss SPO11-1 Rep2",
                   "met1 SPO11-1 Rep1",
                   "met1 SPO11-1 Rep2",
                   "met1 SPO11-1 Rep3",
                   "WT SPO11-1 Rep1",
                   "WT SPO11-1 Rep2",
                   "WT SPO11-1 Rep3"
                  )
log2ChIPColours <- c(
                     rep("black", length(log2ChIPNamesPlot))
                    )
ChIPColours <- log2ChIPColours
ChIPDirs <- sapply(seq_along(ChIPNames), function(x) {
  paste0("/home/ajt200/analysis/",
         ChIPNamesDir[x],
         "/mapped/TEprofiles/matrices/")
})

controlNames <- c(
                  "Col_DMC1_V5_Rep1_ChIP",
                  "Col_DMC1_V5_Rep2_ChIP",
                  "Col_REC8_Myc_Rep1_input",
                  "WT_SPO11oligos_Rep1",
                  "WT_SPO11oligos_Rep2",
                  "WT_SPO11oligos_Rep3",
                  "WT_gDNA_Rep1_R1"
                 )
controlNamesDir <- c(
                     rep(paste0("20190917_dh580_Athaliana_ChIPseq_DMC1/fastq_pooled/snakemake_ChIPseq_", refbase), 3),
                     rep(paste0("160518_Kyuha_SPO11oligos/WT/snakemake_SPO11oligos_", refbase), 3),
                     rep(paste0("150701_Natasha_gDNA/WT/R1/snakemake_SPO11oligos_", refbase), 1)
                    )
controlNamesPlot <- c(
                      "WT DMC1 Rep1",
                      "WT DMC1 Rep2",
                      "WT REC8-Myc Rep1 input",
                      "WT SPO11-1 Rep1",
                      "WT SPO11-1 Rep2",
                      "WT SPO11-1 Rep3",
                      "WT gDNA Rep1 R1"
                     )
controlColours <- c(
                    rep("black", length(controlNamesPlot))
                   )
controlDirs <- sapply(seq_along(controlNames), function(x) {
  paste0("/home/ajt200/analysis/",
         controlNamesDir[x],
         "/mapped/TEprofiles/matrices/")
})

## ChIP
# featureAcc1
ChIP_featureAcc1Mats <- mclapply(seq_along(ChIPNames), function(x) {
  as.matrix(read.table(paste0(ChIPDirs[x],
                              ChIPNames[x],
                              "_MappedOn_", refbase, "_lowXM_", align, "_sort_norm_TEs_", superfamName, "_in_",
                              paste0(chrName, collapse = "_"),
                              "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                       header = F, skip = 3))
}, mc.cores = length(ChIPNames))

# ranLocAcc1
ChIP_ranLocAcc1Mats <- mclapply(seq_along(ChIPNames), function(x) {
  as.matrix(read.table(paste0(ChIPDirs[x],
                              ChIPNames[x],
                              "_MappedOn_", refbase, "_lowXM_", align, "_sort_norm_TEs_", superfamName, "_in_",
                              paste0(chrName, collapse = "_"),
                              "_ranLoc_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                       header = F, skip = 3))
}, mc.cores = length(ChIPNames))

# featureAcc2
ChIP_featureAcc2Mats <- mclapply(seq_along(ChIPNames), function(x) {
  as.matrix(read.table(paste0(ChIPDirs[x],
                              ChIPNames[x],
                              "_MappedOn_", refbase, "_lowXM_", align, "_sort_norm_TEs_", superfamName, "_in_",
                              paste0(chrName, collapse = "_"),
                              "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                       header = F, skip = 3))
}, mc.cores = length(ChIPNames))

# ranLocAcc2
ChIP_ranLocAcc2Mats <- mclapply(seq_along(ChIPNames), function(x) {
  as.matrix(read.table(paste0(ChIPDirs[x],
                              ChIPNames[x],
                              "_MappedOn_", refbase, "_lowXM_", align, "_sort_norm_TEs_", superfamName, "_in_",
                              paste0(chrName, collapse = "_"),
                              "_ranLoc_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                       header = F, skip = 3))
}, mc.cores = length(ChIPNames))


## control
# featureAcc1
control_featureAcc1Mats <- mclapply(seq_along(controlNames), function(x) {
  as.matrix(read.table(paste0(controlDirs[x],
                              controlNames[x],
                              "_MappedOn_", refbase, "_lowXM_", align, "_sort_norm_TEs_", superfamName, "_in_",
                              paste0(chrName, collapse = "_"),
                              "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                       header = F, skip = 3))
}, mc.cores = length(controlNames))

# ranLocAcc1
control_ranLocAcc1Mats <- mclapply(seq_along(controlNames), function(x) {
  as.matrix(read.table(paste0(controlDirs[x],
                              controlNames[x],
                              "_MappedOn_", refbase, "_lowXM_", align, "_sort_norm_TEs_", superfamName, "_in_",
                              paste0(chrName, collapse = "_"),
                              "_ranLoc_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                       header = F, skip = 3))
}, mc.cores = length(controlNames))

# featureAcc2
control_featureAcc2Mats <- mclapply(seq_along(controlNames), function(x) {
  as.matrix(read.table(paste0(controlDirs[x],
                              controlNames[x],
                              "_MappedOn_", refbase, "_lowXM_", align, "_sort_norm_TEs_", superfamName, "_in_",
                              paste0(chrName, collapse = "_"),
                              "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                       header = F, skip = 3))
}, mc.cores = length(controlNames))

# ranLocAcc2
control_ranLocAcc2Mats <- mclapply(seq_along(controlNames), function(x) {
  as.matrix(read.table(paste0(controlDirs[x],
                              controlNames[x],
                              "_MappedOn_", refbase, "_lowXM_", align, "_sort_norm_TEs_", superfamName, "_in_",
                              paste0(chrName, collapse = "_"),
                              "_ranLoc_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                       header = F, skip = 3))
}, mc.cores = length(controlNames))


# Conditionally calculate log2(ChIP/control)
# for each matrix depending on library
# featureAcc1
log2ChIP_featureAcc1Mats <- mclapply(seq_along(ChIP_featureAcc1Mats), function(x) {
  if ( ChIPNames[x] %in% c("cmt3_DMC1_V5_Rep1_ChIP", "kss_DMC1_V5_Rep1_ChIP") ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[which(grepl("Col_DMC1_V5_Rep2_ChIP", controlNames))], " for log2((ChIP+1)/(control+1))"))
    log2((ChIP_featureAcc1Mats[[x]]+1)/(control_featureAcc1Mats[[which(grepl("Col_DMC1_V5_Rep2_ChIP", controlNames))]]+1))
  } else if ( ChIPNames[x] %in% c("Col_DMC1_V5_Rep1_ChIP", "Col_DMC1_V5_Rep2_ChIP") ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[which(grepl("Col_REC8_Myc_Rep1_input", controlNames))], " for log2((ChIP+1)/(control+1))"))
    log2((ChIP_featureAcc1Mats[[x]]+1)/(control_featureAcc1Mats[[which(grepl("Col_REC8_Myc_Rep1_input", controlNames))]]+1))
  } else if ( ChIPNames[x] %in% c("kss_SPO11oligos_Rep1", "kss_SPO11oligos_Rep2", "met1_SPO11oligos_Rep1", "met1_SPO11oligos_Rep2", "met1_SPO11oligos_Rep3") ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[which(grepl("WT_SPO11oligos_Rep3", controlNames))], " for log2((ChIP+1)/(control+1))"))
    log2((ChIP_featureAcc1Mats[[x]]+1)/(control_featureAcc1Mats[[which(grepl("WT_SPO11oligos_Rep3", controlNames))]]+1))
  } else if ( ChIPNames[x] %in% c("WT_SPO11oligos_Rep1", "WT_SPO11oligos_Rep2", "WT_SPO11oligos_Rep3") ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[which(grepl("WT_gDNA_Rep1_R1", controlNames))], " for log2((ChIP+1)/(control+1))"))
    log2((ChIP_featureAcc1Mats[[x]]+1)/(control_featureAcc1Mats[[which(grepl("WT_gDNA_Rep1_R1", controlNames))]]+1))
  }
}, mc.cores = length(ChIP_featureAcc1Mats))

# ranLocAcc1
log2ChIP_ranLocAcc1Mats <- mclapply(seq_along(ChIP_ranLocAcc1Mats), function(x) {
  if ( ChIPNames[x] %in% c("cmt3_DMC1_V5_Rep1_ChIP", "kss_DMC1_V5_Rep1_ChIP") ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[which(grepl("Col_DMC1_V5_Rep2_ChIP", controlNames))], " for log2((ChIP+1)/(control+1))"))
    log2((ChIP_ranLocAcc1Mats[[x]]+1)/(control_ranLocAcc1Mats[[which(grepl("Col_DMC1_V5_Rep2_ChIP", controlNames))]]+1))
  } else if ( ChIPNames[x] %in% c("Col_DMC1_V5_Rep1_ChIP", "Col_DMC1_V5_Rep2_ChIP") ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[which(grepl("Col_REC8_Myc_Rep1_input", controlNames))], " for log2((ChIP+1)/(control+1))"))
    log2((ChIP_ranLocAcc1Mats[[x]]+1)/(control_ranLocAcc1Mats[[which(grepl("Col_REC8_Myc_Rep1_input", controlNames))]]+1))
  } else if ( ChIPNames[x] %in% c("kss_SPO11oligos_Rep1", "kss_SPO11oligos_Rep2", "met1_SPO11oligos_Rep1", "met1_SPO11oligos_Rep2", "met1_SPO11oligos_Rep3") ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[which(grepl("WT_SPO11oligos_Rep3", controlNames))], " for log2((ChIP+1)/(control+1))"))
    log2((ChIP_ranLocAcc1Mats[[x]]+1)/(control_ranLocAcc1Mats[[which(grepl("WT_SPO11oligos_Rep3", controlNames))]]+1))
  } else if ( ChIPNames[x] %in% c("WT_SPO11oligos_Rep1", "WT_SPO11oligos_Rep2", "WT_SPO11oligos_Rep3") ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[which(grepl("WT_gDNA_Rep1_R1", controlNames))], " for log2((ChIP+1)/(control+1))"))
    log2((ChIP_ranLocAcc1Mats[[x]]+1)/(control_ranLocAcc1Mats[[which(grepl("WT_gDNA_Rep1_R1", controlNames))]]+1))
  }
}, mc.cores = length(ChIP_ranLocAcc1Mats))

# Conditionally calculate log2(ChIP/control)
# for each matrix depending on library
# featureAcc2
log2ChIP_featureAcc2Mats <- mclapply(seq_along(ChIP_featureAcc2Mats), function(x) {
  if ( ChIPNames[x] %in% c("cmt3_DMC1_V5_Rep1_ChIP", "kss_DMC1_V5_Rep1_ChIP") ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[which(grepl("Col_DMC1_V5_Rep2_ChIP", controlNames))], " for log2((ChIP+1)/(control+1))"))
    log2((ChIP_featureAcc2Mats[[x]]+1)/(control_featureAcc2Mats[[which(grepl("Col_DMC1_V5_Rep2_ChIP", controlNames))]]+1))
  } else if ( ChIPNames[x] %in% c("Col_DMC1_V5_Rep1_ChIP", "Col_DMC1_V5_Rep2_ChIP") ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[which(grepl("Col_REC8_Myc_Rep1_input", controlNames))], " for log2((ChIP+1)/(control+1))"))
    log2((ChIP_featureAcc2Mats[[x]]+1)/(control_featureAcc2Mats[[which(grepl("Col_REC8_Myc_Rep1_input", controlNames))]]+1))
  } else if ( ChIPNames[x] %in% c("kss_SPO11oligos_Rep1", "kss_SPO11oligos_Rep2", "met1_SPO11oligos_Rep1", "met1_SPO11oligos_Rep2", "met1_SPO11oligos_Rep3") ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[which(grepl("WT_SPO11oligos_Rep3", controlNames))], " for log2((ChIP+1)/(control+1))"))
    log2((ChIP_featureAcc2Mats[[x]]+1)/(control_featureAcc2Mats[[which(grepl("WT_SPO11oligos_Rep3", controlNames))]]+1))
  } else if ( ChIPNames[x] %in% c("WT_SPO11oligos_Rep1", "WT_SPO11oligos_Rep2", "WT_SPO11oligos_Rep3") ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[which(grepl("WT_gDNA_Rep1_R1", controlNames))], " for log2((ChIP+1)/(control+1))"))
    log2((ChIP_featureAcc2Mats[[x]]+1)/(control_featureAcc2Mats[[which(grepl("WT_gDNA_Rep1_R1", controlNames))]]+1))
  }
}, mc.cores = length(ChIP_featureAcc2Mats))

# ranLocAcc2
log2ChIP_ranLocAcc2Mats <- mclapply(seq_along(ChIP_ranLocAcc2Mats), function(x) {
  if ( ChIPNames[x] %in% c("cmt3_DMC1_V5_Rep1_ChIP", "kss_DMC1_V5_Rep1_ChIP") ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[which(grepl("Col_DMC1_V5_Rep2_ChIP", controlNames))], " for log2((ChIP+1)/(control+1))"))
    log2((ChIP_ranLocAcc2Mats[[x]]+1)/(control_ranLocAcc2Mats[[which(grepl("Col_DMC1_V5_Rep2_ChIP", controlNames))]]+1))
  } else if ( ChIPNames[x] %in% c("Col_DMC1_V5_Rep1_ChIP", "Col_DMC1_V5_Rep2_ChIP") ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[which(grepl("Col_REC8_Myc_Rep1_input", controlNames))], " for log2((ChIP+1)/(control+1))"))
    log2((ChIP_ranLocAcc2Mats[[x]]+1)/(control_ranLocAcc2Mats[[which(grepl("Col_REC8_Myc_Rep1_input", controlNames))]]+1))
  } else if ( ChIPNames[x] %in% c("kss_SPO11oligos_Rep1", "kss_SPO11oligos_Rep2", "met1_SPO11oligos_Rep1", "met1_SPO11oligos_Rep2", "met1_SPO11oligos_Rep3") ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[which(grepl("WT_SPO11oligos_Rep3", controlNames))], " for log2((ChIP+1)/(control+1))"))
    log2((ChIP_ranLocAcc2Mats[[x]]+1)/(control_ranLocAcc2Mats[[which(grepl("WT_SPO11oligos_Rep3", controlNames))]]+1))
  } else if ( ChIPNames[x] %in% c("WT_SPO11oligos_Rep1", "WT_SPO11oligos_Rep2", "WT_SPO11oligos_Rep3") ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[which(grepl("WT_gDNA_Rep1_R1", controlNames))], " for log2((ChIP+1)/(control+1))"))
    log2((ChIP_ranLocAcc2Mats[[x]]+1)/(control_ranLocAcc2Mats[[which(grepl("WT_gDNA_Rep1_R1", controlNames))]]+1))
  }
}, mc.cores = length(ChIP_ranLocAcc2Mats))


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
         ChIP_featureAcc2Mats[[x]][quantileIndicesAcc2[[k]],]
       }),
       # featureAcc2 random groupings
       lapply(1:quantiles, function(k) {
         ChIP_featureAcc2Mats[[x]][randomPCIndicesAcc2[[k]],]
       }),
       # ranLocAcc2 groupings
       lapply(1:quantiles, function(k) {
         ChIP_ranLocAcc2Mats[[x]][quantileIndicesAcc2[[k]],]
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

## Repeated ymin and ymax values for as many ChIPNames (and defined as list)
## for consistency with above definitions and convenience
#ymin_list_ChIP <- as.list(rep(min(unlist(ymin_list_ChIP)), length(ymin_list_ChIP)))
#ymax_list_ChIP <- as.list(rep(max(unlist(ymax_list_ChIP)), length(ymax_list_ChIP)))

# Define legend labels
legendLabs_featureAcc1 <- lapply(seq_along(quantileNames), function(x) {
  grobTree(textGrob(bquote(.(quantileNames[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = quantileColoursAcc1[x], fontsize = 18)))
})
legendLabs_ranFeatAcc1 <- lapply(seq_along(randomPCNames), function(x) {
  grobTree(textGrob(bquote(.(randomPCNames[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = quantileColoursAcc1[x], fontsize = 18)))
})
legendLabs_ranLocAcc1 <- lapply(seq_along(randomPCNames), function(x) {
  grobTree(textGrob(bquote(.(randomPCNames[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = quantileColoursAcc1[x], fontsize = 18)))
})
legendLabs_featureAcc2 <- lapply(seq_along(quantileNames), function(x) {
  grobTree(textGrob(bquote(.(quantileNames[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = quantileColoursAcc2[x], fontsize = 18)))
})
legendLabs_ranFeatAcc2 <- lapply(seq_along(randomPCNames), function(x) {
  grobTree(textGrob(bquote(.(randomPCNames[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = quantileColoursAcc2[x], fontsize = 18)))
})
legendLabs_ranLocAcc2 <- lapply(seq_along(randomPCNames), function(x) {
  grobTree(textGrob(bquote(.(randomPCNames[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = quantileColoursAcc2[x], fontsize = 18)))
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
  scale_colour_manual(values = quantileColoursAcc1) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColoursAcc1) +
  scale_y_continuous(limits = c(ymin_list_ChIP[[x]], ymax_list_ChIP[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_ChIP[[x]][[1]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_ChIP[[x]][[1]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
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
  scale_colour_manual(values = quantileColoursAcc1) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColoursAcc1) +
  scale_y_continuous(limits = c(ymin_list_ChIP[[x]], ymax_list_ChIP[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_ChIP[[x]][[2]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_ChIP[[x]][[2]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
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
  scale_colour_manual(values = quantileColoursAcc1) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColoursAcc1) +
  scale_y_continuous(limits = c(ymin_list_ChIP[[x]], ymax_list_ChIP[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_ChIP[[x]][[3]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_ChIP[[x]][[3]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              "Start",
                              "End",
                              paste0("+", flankNamePlot))) +
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
  scale_colour_manual(values = quantileColoursAcc2) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColoursAcc2) +
  scale_y_continuous(limits = c(ymin_list_ChIP[[x]], ymax_list_ChIP[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_ChIP[[x]][[4]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_ChIP[[x]][[4]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
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
  scale_colour_manual(values = quantileColoursAcc2) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColoursAcc2) +
  scale_y_continuous(limits = c(ymin_list_ChIP[[x]], ymax_list_ChIP[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_ChIP[[x]][[5]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_ChIP[[x]][[5]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
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
  scale_colour_manual(values = quantileColoursAcc2) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColoursAcc2) +
  scale_y_continuous(limits = c(ymin_list_ChIP[[x]], ymax_list_ChIP[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_ChIP[[x]][[6]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_ChIP[[x]][[6]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              "Start",
                              "End",
                              paste0("+", flankNamePlot))) +
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
              "ChIP_", align, "_avgProfiles_around_features_", quantiles, "quantiles",
              "_by_", context, "_absolute_change_in_", libName1[1], "_", libName2[1],
              "_at_", superfamName, "_", sortRegion,
              "_in_", refbase, "_",
              paste0(chrName, collapse = "_"), ".pdf"),
       plot = ggObjGA_combined,
       height = 6.5*length(c(ChIPNamesPlot)), width = 7*6, limitsize = FALSE)

#### Free up memory by removing no longer required objects
rm(
   ChIP_featureAcc1Mats, ChIP_ranLocAcc1Mats,
   ChIP_featureAcc2Mats, ChIP_ranLocAcc2Mats,
   ChIP_mats_quantiles,
   wideDFfeature_list_ChIP,
   tidyDFfeature_list_ChIP,
   summaryDFfeature_list_ChIP,
   summaryDFfeature_ChIP
  ) 
gc()
#####


# control
# Add column names
for(x in seq_along(control_featureAcc1Mats)) {
  colnames(control_featureAcc1Mats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                           paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                           paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
  colnames(control_ranLocAcc1Mats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                          paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                          paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
  colnames(control_featureAcc2Mats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                           paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                           paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
  colnames(control_ranLocAcc2Mats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                          paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                          paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
}

# Subdivide coverage matrices into above-defined quantiles and random groupings
control_mats_quantiles <- mclapply(seq_along(control_featureAcc1Mats), function(x) {
  list(
       # featureAcc1 quantiles
       lapply(1:quantiles, function(k) {
         control_featureAcc1Mats[[x]][quantileIndicesAcc1[[k]],]
       }),
       # featureAcc1 random groupings
       lapply(1:quantiles, function(k) {
         control_featureAcc1Mats[[x]][randomPCIndicesAcc1[[k]],]
       }),
       # ranLocAcc1 groupings
       lapply(1:quantiles, function(k) {
         control_ranLocAcc1Mats[[x]][quantileIndicesAcc1[[k]],]
       }),
       # featureAcc2 quantiles
       lapply(1:quantiles, function(k) {
         control_featureAcc2Mats[[x]][quantileIndicesAcc2[[k]],]
       }),
       # featureAcc2 random groupings
       lapply(1:quantiles, function(k) {
         control_featureAcc2Mats[[x]][randomPCIndicesAcc2[[k]],]
       }),
       # ranLocAcc2 groupings
       lapply(1:quantiles, function(k) {
         control_ranLocAcc2Mats[[x]][quantileIndicesAcc2[[k]],]
       })
      ) 
}, mc.cores = length(control_featureAcc1Mats))

# Transpose matrix and convert into dataframe
# in which first column is window name
wideDFfeature_list_control <- mclapply(seq_along(control_mats_quantiles), function(x) {
  lapply(seq_along(control_mats_quantiles[[x]]), function(y) {
    lapply(seq_along(control_mats_quantiles[[x]][[y]]), function(k) {
      data.frame(window = colnames(control_mats_quantiles[[x]][[y]][[k]]),
                 t(control_mats_quantiles[[x]][[y]][[k]]))
    })
  })
}, mc.cores = length(control_mats_quantiles))

# Convert into tidy data.frame (long format)
tidyDFfeature_list_control  <- mclapply(seq_along(wideDFfeature_list_control), function(x) {
  lapply(seq_along(control_mats_quantiles[[x]]), function(y) {
    lapply(seq_along(control_mats_quantiles[[x]][[y]]), function(k) {
      gather(data  = wideDFfeature_list_control[[x]][[y]][[k]],
             key   = feature,
             value = coverage,
             -window)
    })
  }) 
}, mc.cores = length(wideDFfeature_list_control))

# Order levels of factor "window" so that sequential levels
# correspond to sequential windows
for(x in seq_along(tidyDFfeature_list_control)) {
  for(y in seq_along(control_mats_quantiles[[x]])) {
    for(k in seq_along(control_mats_quantiles[[x]][[y]])) {
      tidyDFfeature_list_control[[x]][[y]][[k]]$window <- factor(tidyDFfeature_list_control[[x]][[y]][[k]]$window,
                                                                  levels = as.character(wideDFfeature_list_control[[x]][[y]][[k]]$window))
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
summaryDFfeature_list_control  <- mclapply(seq_along(tidyDFfeature_list_control), function(x) {
  lapply(seq_along(control_mats_quantiles[[x]]), function(y) {
    lapply(seq_along(control_mats_quantiles[[x]][[y]]), function(k) {
      data.frame(window = as.character(wideDFfeature_list_control[[x]][[y]][[k]]$window),
                 n      = tapply(X     = tidyDFfeature_list_control[[x]][[y]][[k]]$coverage,
                                 INDEX = tidyDFfeature_list_control[[x]][[y]][[k]]$window,
                                 FUN   = length),
                 mean   = tapply(X     = tidyDFfeature_list_control[[x]][[y]][[k]]$coverage,
                                 INDEX = tidyDFfeature_list_control[[x]][[y]][[k]]$window,
                                 FUN   = mean,
                                 na.rm = TRUE),
                 sd     = tapply(X     = tidyDFfeature_list_control[[x]][[y]][[k]]$coverage,
                                 INDEX = tidyDFfeature_list_control[[x]][[y]][[k]]$window,
                                 FUN   = sd,
                                 na.rm = TRUE))
    })
  })
}, mc.cores = length(tidyDFfeature_list_control))

for(x in seq_along(summaryDFfeature_list_control)) {
  for(y in seq_along(control_mats_quantiles[[x]])) {
    for(k in seq_along(control_mats_quantiles[[x]][[y]])) {
      summaryDFfeature_list_control[[x]][[y]][[k]]$window <- factor(summaryDFfeature_list_control[[x]][[y]][[k]]$window,
                                                                     levels = as.character(wideDFfeature_list_control[[x]][[y]][[k]]$window))
      summaryDFfeature_list_control[[x]][[y]][[k]]$winNo <- factor(1:dim(summaryDFfeature_list_control[[x]][[y]][[k]])[1])
      summaryDFfeature_list_control[[x]][[y]][[k]]$sem <- summaryDFfeature_list_control[[x]][[y]][[k]]$sd/sqrt(summaryDFfeature_list_control[[x]][[y]][[k]]$n-1)
      summaryDFfeature_list_control[[x]][[y]][[k]]$CI_lower <- summaryDFfeature_list_control[[x]][[y]][[k]]$mean -
        qt(0.975, df = summaryDFfeature_list_control[[x]][[y]][[k]]$n-1)*summaryDFfeature_list_control[[x]][[y]][[k]]$sem
      summaryDFfeature_list_control[[x]][[y]][[k]]$CI_upper <- summaryDFfeature_list_control[[x]][[y]][[k]]$mean +
        qt(0.975, df = summaryDFfeature_list_control[[x]][[y]][[k]]$n-1)*summaryDFfeature_list_control[[x]][[y]][[k]]$sem
    }
  }
}

quantileNames <- paste0(rep("Quantile ", quantiles), 1:quantiles)
randomPCNames <- paste0(rep("Random ", quantiles), 1:quantiles)
for(x in seq_along(summaryDFfeature_list_control)) {
  # featureAcc1 quantiles
  names(summaryDFfeature_list_control[[x]][[1]]) <- quantileNames
  # featureAcc1 random groupings
  names(summaryDFfeature_list_control[[x]][[2]]) <- randomPCNames
  # ranLocAcc1 groupings
  names(summaryDFfeature_list_control[[x]][[3]]) <- randomPCNames
  # featureAcc2 quantiles
  names(summaryDFfeature_list_control[[x]][[4]]) <- quantileNames
  # featureAcc2 random groupings
  names(summaryDFfeature_list_control[[x]][[5]]) <- randomPCNames
  # ranLocAcc2 groupings
  names(summaryDFfeature_list_control[[x]][[6]]) <- randomPCNames
}

# Convert list of lists of lists of feature quantiles summaryDFfeature_list_control into
# a list of lists of single data.frames containing all feature quantiles for plotting
summaryDFfeature_control  <- mclapply(seq_along(summaryDFfeature_list_control), function(x) {
  lapply(seq_along(control_mats_quantiles[[x]]), function(y) {
    bind_rows(summaryDFfeature_list_control[[x]][[y]], .id = "quantile")
  })
}, mc.cores = length(summaryDFfeature_list_control))
for(x in seq_along(summaryDFfeature_control)) {
  # featureAcc1 quantiles
  summaryDFfeature_control[[x]][[1]]$quantile <- factor(summaryDFfeature_control[[x]][[1]]$quantile,
                                                     levels = names(summaryDFfeature_list_control[[x]][[1]]))
  # featureAcc1 random groupings
  summaryDFfeature_control[[x]][[2]]$quantile <- factor(summaryDFfeature_control[[x]][[2]]$quantile,
                                                     levels = names(summaryDFfeature_list_control[[x]][[2]]))
  # ranLocAcc1 groupings
  summaryDFfeature_control[[x]][[3]]$quantile <- factor(summaryDFfeature_control[[x]][[3]]$quantile,
                                                     levels = names(summaryDFfeature_list_control[[x]][[3]]))
  # featureAcc2 quantiles
  summaryDFfeature_control[[x]][[4]]$quantile <- factor(summaryDFfeature_control[[x]][[4]]$quantile,
                                                     levels = names(summaryDFfeature_list_control[[x]][[4]]))
  # featureAcc2 random groupings
  summaryDFfeature_control[[x]][[5]]$quantile <- factor(summaryDFfeature_control[[x]][[5]]$quantile,
                                                     levels = names(summaryDFfeature_list_control[[x]][[5]]))
  # ranLocAcc2 groupings
  summaryDFfeature_control[[x]][[6]]$quantile <- factor(summaryDFfeature_control[[x]][[6]]$quantile,
                                                     levels = names(summaryDFfeature_list_control[[x]][[6]]))
}

# Define y-axis limits
ymin_list_control <- lapply(seq_along(summaryDFfeature_control), function(x) {
  min(c(summaryDFfeature_control[[x]][[1]]$CI_lower,
        summaryDFfeature_control[[x]][[2]]$CI_lower,
        summaryDFfeature_control[[x]][[3]]$CI_lower,
        summaryDFfeature_control[[x]][[4]]$CI_lower,
        summaryDFfeature_control[[x]][[5]]$CI_lower,
        summaryDFfeature_control[[x]][[6]]$CI_lower))
})
ymax_list_control <- lapply(seq_along(summaryDFfeature_control), function(x) {
  max(c(summaryDFfeature_control[[x]][[1]]$CI_upper,
        summaryDFfeature_control[[x]][[2]]$CI_upper,
        summaryDFfeature_control[[x]][[3]]$CI_upper,
        summaryDFfeature_control[[x]][[4]]$CI_upper,
        summaryDFfeature_control[[x]][[5]]$CI_upper,
        summaryDFfeature_control[[x]][[6]]$CI_upper))
})

## Repeated ymin and ymax values for as many controlNames (and defined as list)
## for consistency with above definitions and convenience
#ymin_list_control <- as.list(rep(min(unlist(ymin_list_control)), length(ymin_list_control)))
#ymax_list_control <- as.list(rep(max(unlist(ymax_list_control)), length(ymax_list_control)))

# Define legend labels
legendLabs_featureAcc1 <- lapply(seq_along(quantileNames), function(x) {
  grobTree(textGrob(bquote(.(quantileNames[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = quantileColoursAcc1[x], fontsize = 18)))
})
legendLabs_ranFeatAcc1 <- lapply(seq_along(randomPCNames), function(x) {
  grobTree(textGrob(bquote(.(randomPCNames[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = quantileColoursAcc1[x], fontsize = 18)))
})
legendLabs_ranLocAcc1 <- lapply(seq_along(randomPCNames), function(x) {
  grobTree(textGrob(bquote(.(randomPCNames[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = quantileColoursAcc1[x], fontsize = 18)))
})
legendLabs_featureAcc2 <- lapply(seq_along(quantileNames), function(x) {
  grobTree(textGrob(bquote(.(quantileNames[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = quantileColoursAcc2[x], fontsize = 18)))
})
legendLabs_ranFeatAcc2 <- lapply(seq_along(randomPCNames), function(x) {
  grobTree(textGrob(bquote(.(randomPCNames[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = quantileColoursAcc2[x], fontsize = 18)))
})
legendLabs_ranLocAcc2 <- lapply(seq_along(randomPCNames), function(x) {
  grobTree(textGrob(bquote(.(randomPCNames[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = quantileColoursAcc2[x], fontsize = 18)))
})

# Plot average profiles with 95% CI ribbon
## feature
ggObj1_combined_control <- mclapply(seq_along(controlNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_control[[x]][[1]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColoursAcc1) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColoursAcc1) +
  scale_y_continuous(limits = c(ymin_list_control[[x]], ymax_list_control[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_control[[x]][[1]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_control[[x]][[1]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_control[[x]][[1]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = controlNamesPlot[x]) +
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
        axis.title = element_text(size = 30, colour = controlColours[x]),
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
}, mc.cores = length(controlNamesPlot))

## ranFeat
ggObj2_combined_control <- mclapply(seq_along(controlNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_control[[x]][[2]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColoursAcc1) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColoursAcc1) +
  scale_y_continuous(limits = c(ymin_list_control[[x]], ymax_list_control[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_control[[x]][[2]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_control[[x]][[2]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_control[[x]][[2]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = controlNamesPlot[x]) +
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
        axis.title = element_text(size = 30, colour = controlColours[x]),
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
}, mc.cores = length(controlNamesPlot))

## ranLoc
ggObj3_combined_control <- mclapply(seq_along(controlNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_control[[x]][[3]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColoursAcc1) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColoursAcc1) +
  scale_y_continuous(limits = c(ymin_list_control[[x]], ymax_list_control[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_control[[x]][[3]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_control[[x]][[3]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              "Start",
                              "End",
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_control[[x]][[3]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = controlNamesPlot[x]) +
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
        axis.title = element_text(size = 30, colour = controlColours[x]),
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
}, mc.cores = length(controlNamesPlot))

## feature
ggObj4_combined_control <- mclapply(seq_along(controlNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_control[[x]][[4]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColoursAcc2) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColoursAcc2) +
  scale_y_continuous(limits = c(ymin_list_control[[x]], ymax_list_control[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_control[[x]][[4]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_control[[x]][[4]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_control[[x]][[4]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = controlNamesPlot[x]) +
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
        axis.title = element_text(size = 30, colour = controlColours[x]),
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
}, mc.cores = length(controlNamesPlot))

## ranFeat
ggObj5_combined_control <- mclapply(seq_along(controlNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_control[[x]][[5]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColoursAcc2) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColoursAcc2) +
  scale_y_continuous(limits = c(ymin_list_control[[x]], ymax_list_control[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_control[[x]][[5]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_control[[x]][[5]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_control[[x]][[5]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = controlNamesPlot[x]) +
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
        axis.title = element_text(size = 30, colour = controlColours[x]),
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
}, mc.cores = length(controlNamesPlot))

## ranLoc
ggObj6_combined_control <- mclapply(seq_along(controlNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_control[[x]][[6]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColoursAcc2) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColoursAcc2) +
  scale_y_continuous(limits = c(ymin_list_control[[x]], ymax_list_control[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_control[[x]][[6]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_control[[x]][[6]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              "Start",
                              "End",
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_control[[x]][[6]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = controlNamesPlot[x]) +
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
        axis.title = element_text(size = 30, colour = controlColours[x]),
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
}, mc.cores = length(controlNamesPlot))

ggObjGA_combined <- grid.arrange(grobs = c(
                                           ggObj1_combined_control,
                                           ggObj2_combined_control,
                                           ggObj3_combined_control,
                                           ggObj4_combined_control,
                                           ggObj5_combined_control,
                                           ggObj6_combined_control
                                          ),
                                 layout_matrix = cbind(
                                                       1:length(c(controlNamesPlot)),
                                                       (length(c(controlNamesPlot))+1):(length(c(controlNamesPlot))*2),
                                                       ((length(c(controlNamesPlot))*2)+1):(length(c(controlNamesPlot))*3),
                                                       ((length(c(controlNamesPlot))*3)+1):(length(c(controlNamesPlot))*4),
                                                       ((length(c(controlNamesPlot))*4)+1):(length(c(controlNamesPlot))*5),
                                                       ((length(c(controlNamesPlot))*5)+1):(length(c(controlNamesPlot))*6)
                                                      ))
ggsave(paste0(plotDir,
              "control_", align, "_avgProfiles_around_features_", quantiles, "quantiles",
              "_by_", context, "_absolute_change_in_", libName1[1], "_", libName2[1],
              "_at_", superfamName, "_", sortRegion,
              "_in_", refbase, "_",
              paste0(chrName, collapse = "_"), ".pdf"),
       plot = ggObjGA_combined,
       height = 6.5*length(c(controlNamesPlot)), width = 7*6, limitsize = FALSE)

#### Free up memory by removing no longer required objects
rm(
   control_featureAcc1Mats, control_ranLocAcc1Mats,
   control_featureAcc2Mats, control_ranLocAcc2Mats,
   control_mats_quantiles,
   wideDFfeature_list_control,
   tidyDFfeature_list_control,
   summaryDFfeature_list_control,
   summaryDFfeature_control
  ) 
gc()
#####


# log2ChIP
# Add column names
for(x in seq_along(log2ChIP_featureAcc1Mats)) {
  colnames(log2ChIP_featureAcc1Mats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                               paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                               paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
  colnames(log2ChIP_ranLocAcc1Mats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                              paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                              paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
  colnames(log2ChIP_featureAcc2Mats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                               paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                               paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
  colnames(log2ChIP_ranLocAcc2Mats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                              paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                              paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
}

# Subdivide coverage matrices into above-defined quantiles and random groupings
log2ChIP_mats_quantiles <- mclapply(seq_along(log2ChIP_featureAcc1Mats), function(x) {
  list(
       # featureAcc1 quantiles
       lapply(1:quantiles, function(k) {
         log2ChIP_featureAcc1Mats[[x]][quantileIndicesAcc1[[k]],]
       }),
       # featureAcc1 random groupings
       lapply(1:quantiles, function(k) {
         log2ChIP_featureAcc1Mats[[x]][randomPCIndicesAcc1[[k]],]
       }),
       # ranLocAcc1 groupings
       lapply(1:quantiles, function(k) {
         log2ChIP_ranLocAcc1Mats[[x]][quantileIndicesAcc1[[k]],]
       }),
       # featureAcc2 quantiles
       lapply(1:quantiles, function(k) {
         log2ChIP_featureAcc2Mats[[x]][quantileIndicesAcc2[[k]],]
       }),
       # featureAcc2 random groupings
       lapply(1:quantiles, function(k) {
         log2ChIP_featureAcc2Mats[[x]][randomPCIndicesAcc2[[k]],]
       }),
       # ranLocAcc2 groupings
       lapply(1:quantiles, function(k) {
         log2ChIP_ranLocAcc2Mats[[x]][quantileIndicesAcc2[[k]],]
       })
      ) 
}, mc.cores = length(log2ChIP_featureAcc1Mats))

# Transpose matrix and convert into dataframe
# in which first column is window name
wideDFfeature_list_log2ChIP <- mclapply(seq_along(log2ChIP_mats_quantiles), function(x) {
  lapply(seq_along(log2ChIP_mats_quantiles[[x]]), function(y) {
    lapply(seq_along(log2ChIP_mats_quantiles[[x]][[y]]), function(k) {
      data.frame(window = colnames(log2ChIP_mats_quantiles[[x]][[y]][[k]]),
                 t(log2ChIP_mats_quantiles[[x]][[y]][[k]]))
    })
  })
}, mc.cores = length(log2ChIP_mats_quantiles))

# Convert into tidy data.frame (long format)
tidyDFfeature_list_log2ChIP  <- mclapply(seq_along(wideDFfeature_list_log2ChIP), function(x) {
  lapply(seq_along(log2ChIP_mats_quantiles[[x]]), function(y) {
    lapply(seq_along(log2ChIP_mats_quantiles[[x]][[y]]), function(k) {
      gather(data  = wideDFfeature_list_log2ChIP[[x]][[y]][[k]],
             key   = feature,
             value = coverage,
             -window)
    })
  }) 
}, mc.cores = length(wideDFfeature_list_log2ChIP))

# Order levels of factor "window" so that sequential levels
# correspond to sequential windows
for(x in seq_along(tidyDFfeature_list_log2ChIP)) {
  for(y in seq_along(log2ChIP_mats_quantiles[[x]])) {
    for(k in seq_along(log2ChIP_mats_quantiles[[x]][[y]])) {
      tidyDFfeature_list_log2ChIP[[x]][[y]][[k]]$window <- factor(tidyDFfeature_list_log2ChIP[[x]][[y]][[k]]$window,
                                                                  levels = as.character(wideDFfeature_list_log2ChIP[[x]][[y]][[k]]$window))
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
summaryDFfeature_list_log2ChIP  <- mclapply(seq_along(tidyDFfeature_list_log2ChIP), function(x) {
  lapply(seq_along(log2ChIP_mats_quantiles[[x]]), function(y) {
    lapply(seq_along(log2ChIP_mats_quantiles[[x]][[y]]), function(k) {
      data.frame(window = as.character(wideDFfeature_list_log2ChIP[[x]][[y]][[k]]$window),
                 n      = tapply(X     = tidyDFfeature_list_log2ChIP[[x]][[y]][[k]]$coverage,
                                 INDEX = tidyDFfeature_list_log2ChIP[[x]][[y]][[k]]$window,
                                 FUN   = length),
                 mean   = tapply(X     = tidyDFfeature_list_log2ChIP[[x]][[y]][[k]]$coverage,
                                 INDEX = tidyDFfeature_list_log2ChIP[[x]][[y]][[k]]$window,
                                 FUN   = mean,
                                 na.rm = TRUE),
                 sd     = tapply(X     = tidyDFfeature_list_log2ChIP[[x]][[y]][[k]]$coverage,
                                 INDEX = tidyDFfeature_list_log2ChIP[[x]][[y]][[k]]$window,
                                 FUN   = sd,
                                 na.rm = TRUE))
    })
  })
}, mc.cores = length(tidyDFfeature_list_log2ChIP))

for(x in seq_along(summaryDFfeature_list_log2ChIP)) {
  for(y in seq_along(log2ChIP_mats_quantiles[[x]])) {
    for(k in seq_along(log2ChIP_mats_quantiles[[x]][[y]])) {
      summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$window <- factor(summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$window,
                                                                     levels = as.character(wideDFfeature_list_log2ChIP[[x]][[y]][[k]]$window))
      summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$winNo <- factor(1:dim(summaryDFfeature_list_log2ChIP[[x]][[y]][[k]])[1])
      summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$sem <- summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$sd/sqrt(summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$n-1)
      summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$CI_lower <- summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$mean -
        qt(0.975, df = summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$n-1)*summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$sem
      summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$CI_upper <- summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$mean +
        qt(0.975, df = summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$n-1)*summaryDFfeature_list_log2ChIP[[x]][[y]][[k]]$sem
    }
  }
}

quantileNames <- paste0(rep("Quantile ", quantiles), 1:quantiles)
randomPCNames <- paste0(rep("Random ", quantiles), 1:quantiles)
for(x in seq_along(summaryDFfeature_list_log2ChIP)) {
  # featureAcc1 quantiles
  names(summaryDFfeature_list_log2ChIP[[x]][[1]]) <- quantileNames
  # featureAcc1 random groupings
  names(summaryDFfeature_list_log2ChIP[[x]][[2]]) <- randomPCNames
  # ranLocAcc1 groupings
  names(summaryDFfeature_list_log2ChIP[[x]][[3]]) <- randomPCNames
  # featureAcc2 quantiles
  names(summaryDFfeature_list_log2ChIP[[x]][[4]]) <- quantileNames
  # featureAcc2 random groupings
  names(summaryDFfeature_list_log2ChIP[[x]][[5]]) <- randomPCNames
  # ranLocAcc2 groupings
  names(summaryDFfeature_list_log2ChIP[[x]][[6]]) <- randomPCNames
}

# Convert list of lists of lists of feature quantiles summaryDFfeature_list_log2ChIP into
# a list of lists of single data.frames containing all feature quantiles for plotting
summaryDFfeature_log2ChIP  <- mclapply(seq_along(summaryDFfeature_list_log2ChIP), function(x) {
  lapply(seq_along(log2ChIP_mats_quantiles[[x]]), function(y) {
    bind_rows(summaryDFfeature_list_log2ChIP[[x]][[y]], .id = "quantile")
  })
}, mc.cores = length(summaryDFfeature_list_log2ChIP))
for(x in seq_along(summaryDFfeature_log2ChIP)) {
  # featureAcc1 quantiles
  summaryDFfeature_log2ChIP[[x]][[1]]$quantile <- factor(summaryDFfeature_log2ChIP[[x]][[1]]$quantile,
                                                         levels = names(summaryDFfeature_list_log2ChIP[[x]][[1]]))
  # featureAcc1 random groupings
  summaryDFfeature_log2ChIP[[x]][[2]]$quantile <- factor(summaryDFfeature_log2ChIP[[x]][[2]]$quantile,
                                                         levels = names(summaryDFfeature_list_log2ChIP[[x]][[2]]))
  # ranLocAcc1 groupings
  summaryDFfeature_log2ChIP[[x]][[3]]$quantile <- factor(summaryDFfeature_log2ChIP[[x]][[3]]$quantile,
                                                         levels = names(summaryDFfeature_list_log2ChIP[[x]][[3]]))
  # featureAcc2 quantiles
  summaryDFfeature_log2ChIP[[x]][[4]]$quantile <- factor(summaryDFfeature_log2ChIP[[x]][[4]]$quantile,
                                                         levels = names(summaryDFfeature_list_log2ChIP[[x]][[4]]))
  # featureAcc2 random groupings
  summaryDFfeature_log2ChIP[[x]][[5]]$quantile <- factor(summaryDFfeature_log2ChIP[[x]][[5]]$quantile,
                                                         levels = names(summaryDFfeature_list_log2ChIP[[x]][[5]]))
  # ranLocAcc2 groupings
  summaryDFfeature_log2ChIP[[x]][[6]]$quantile <- factor(summaryDFfeature_log2ChIP[[x]][[6]]$quantile,
                                                         levels = names(summaryDFfeature_list_log2ChIP[[x]][[6]]))
}

# Define y-axis limits
ymin_list_log2ChIP <- lapply(seq_along(summaryDFfeature_log2ChIP), function(x) {
  min(c(summaryDFfeature_log2ChIP[[x]][[1]]$CI_lower,
        summaryDFfeature_log2ChIP[[x]][[2]]$CI_lower,
        summaryDFfeature_log2ChIP[[x]][[3]]$CI_lower,
        summaryDFfeature_log2ChIP[[x]][[4]]$CI_lower,
        summaryDFfeature_log2ChIP[[x]][[5]]$CI_lower,
        summaryDFfeature_log2ChIP[[x]][[6]]$CI_lower))
})
ymax_list_log2ChIP <- lapply(seq_along(summaryDFfeature_log2ChIP), function(x) {
  max(c(summaryDFfeature_log2ChIP[[x]][[1]]$CI_upper,
        summaryDFfeature_log2ChIP[[x]][[2]]$CI_upper,
        summaryDFfeature_log2ChIP[[x]][[3]]$CI_upper,
        summaryDFfeature_log2ChIP[[x]][[4]]$CI_upper,
        summaryDFfeature_log2ChIP[[x]][[5]]$CI_upper,
        summaryDFfeature_log2ChIP[[x]][[6]]$CI_upper))
})

## Repeated ymin and ymax values for as many log2ChIPNames (and defined as list)
## for consistency with above definitions and convenience
#ymin_list_log2ChIP <- as.list(rep(min(unlist(ymin_list_log2ChIP)), length(ymin_list_log2ChIP)))
#ymax_list_log2ChIP <- as.list(rep(max(unlist(ymax_list_log2ChIP)), length(ymax_list_log2ChIP)))

# Define legend labels
legendLabs_featureAcc1 <- lapply(seq_along(quantileNames), function(x) {
  grobTree(textGrob(bquote(.(quantileNames[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = quantileColoursAcc1[x], fontsize = 18)))
})
legendLabs_ranFeatAcc1 <- lapply(seq_along(randomPCNames), function(x) {
  grobTree(textGrob(bquote(.(randomPCNames[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = quantileColoursAcc1[x], fontsize = 18)))
})
legendLabs_ranLocAcc1 <- lapply(seq_along(randomPCNames), function(x) {
  grobTree(textGrob(bquote(.(randomPCNames[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = quantileColoursAcc1[x], fontsize = 18)))
})
legendLabs_featureAcc2 <- lapply(seq_along(quantileNames), function(x) {
  grobTree(textGrob(bquote(.(quantileNames[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = quantileColoursAcc2[x], fontsize = 18)))
})
legendLabs_ranFeatAcc2 <- lapply(seq_along(randomPCNames), function(x) {
  grobTree(textGrob(bquote(.(randomPCNames[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = quantileColoursAcc2[x], fontsize = 18)))
})
legendLabs_ranLocAcc2 <- lapply(seq_along(randomPCNames), function(x) {
  grobTree(textGrob(bquote(.(randomPCNames[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = quantileColoursAcc2[x], fontsize = 18)))
})

# Plot average profiles with 95% CI ribbon
## feature
ggObj1_combined_log2ChIP <- mclapply(seq_along(log2ChIPNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_log2ChIP[[x]][[1]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColoursAcc1) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColoursAcc1) +
  scale_y_continuous(limits = c(ymin_list_log2ChIP[[x]], ymax_list_log2ChIP[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_log2ChIP[[x]][[1]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_log2ChIP[[x]][[1]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[x]][[1]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = bquote("Log"[2] * .(log2ChIPNamesPlot[x]))) +
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
        axis.title = element_text(size = 30, colour = log2ChIPColours[x]),
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
}, mc.cores = length(log2ChIPNamesPlot))

## ranFeat
ggObj2_combined_log2ChIP <- mclapply(seq_along(log2ChIPNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_log2ChIP[[x]][[2]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColoursAcc1) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColoursAcc1) +
  scale_y_continuous(limits = c(ymin_list_log2ChIP[[x]], ymax_list_log2ChIP[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_log2ChIP[[x]][[2]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_log2ChIP[[x]][[2]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[x]][[2]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = bquote("Log"[2] * .(log2ChIPNamesPlot[x]))) +
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
        axis.title = element_text(size = 30, colour = log2ChIPColours[x]),
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
}, mc.cores = length(log2ChIPNamesPlot))

## ranLoc
ggObj3_combined_log2ChIP <- mclapply(seq_along(log2ChIPNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_log2ChIP[[x]][[3]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColoursAcc1) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColoursAcc1) +
  scale_y_continuous(limits = c(ymin_list_log2ChIP[[x]], ymax_list_log2ChIP[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_log2ChIP[[x]][[3]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_log2ChIP[[x]][[3]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              "Start",
                              "End",
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[x]][[3]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = bquote("Log"[2] * .(log2ChIPNamesPlot[x]))) +
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
        axis.title = element_text(size = 30, colour = log2ChIPColours[x]),
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
}, mc.cores = length(log2ChIPNamesPlot))

## feature
ggObj4_combined_log2ChIP <- mclapply(seq_along(log2ChIPNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_log2ChIP[[x]][[4]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColoursAcc2) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColoursAcc2) +
  scale_y_continuous(limits = c(ymin_list_log2ChIP[[x]], ymax_list_log2ChIP[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_log2ChIP[[x]][[4]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_log2ChIP[[x]][[4]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[x]][[4]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = bquote("Log"[2] * .(log2ChIPNamesPlot[x]))) +
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
        axis.title = element_text(size = 30, colour = log2ChIPColours[x]),
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
}, mc.cores = length(log2ChIPNamesPlot))

## ranFeat
ggObj5_combined_log2ChIP <- mclapply(seq_along(log2ChIPNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_log2ChIP[[x]][[5]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColoursAcc2) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColoursAcc2) +
  scale_y_continuous(limits = c(ymin_list_log2ChIP[[x]], ymax_list_log2ChIP[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_log2ChIP[[x]][[5]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_log2ChIP[[x]][[5]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[x]][[5]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = bquote("Log"[2] * .(log2ChIPNamesPlot[x]))) +
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
        axis.title = element_text(size = 30, colour = log2ChIPColours[x]),
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
}, mc.cores = length(log2ChIPNamesPlot))

## ranLoc
ggObj6_combined_log2ChIP <- mclapply(seq_along(log2ChIPNamesPlot), function(x) {
  summaryDFfeature <- summaryDFfeature_log2ChIP[[x]][[6]]
  ggplot(data = summaryDFfeature,
         mapping = aes(x = winNo,
                       y = mean,
                       group = quantile)
        ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = quantile),
            size = 1) +
  scale_colour_manual(values = quantileColoursAcc2) +
  geom_ribbon(data = summaryDFfeature,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = quantile),
              alpha = 0.4) +
  scale_fill_manual(values = quantileColoursAcc2) +
  scale_y_continuous(limits = c(ymin_list_log2ChIP[[x]], ymax_list_log2ChIP[[x]]),
                     labels = function(x) sprintf("%6.3f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_log2ChIP[[x]][[6]])[1]/quantiles)-(downstream/binSize),
                              dim(summaryDFfeature_log2ChIP[[x]][[6]])[1]/quantiles),
                   labels = c(paste0("-", flankNamePlot),
                              "Start",
                              "End",
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[x]][[6]])[1]/quantiles)-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = bquote("Log"[2] * .(log2ChIPNamesPlot[x]))) +
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
        axis.title = element_text(size = 30, colour = log2ChIPColours[x]),
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
}, mc.cores = length(log2ChIPNamesPlot))

ggObjGA_combined <- grid.arrange(grobs = c(
                                           ggObj1_combined_log2ChIP,
                                           ggObj2_combined_log2ChIP,
                                           ggObj3_combined_log2ChIP,
                                           ggObj4_combined_log2ChIP,
                                           ggObj5_combined_log2ChIP,
                                           ggObj6_combined_log2ChIP
                                          ),
                                 layout_matrix = cbind(
                                                       1:length(c(log2ChIPNamesPlot)),
                                                       (length(c(log2ChIPNamesPlot))+1):(length(c(log2ChIPNamesPlot))*2),
                                                       ((length(c(log2ChIPNamesPlot))*2)+1):(length(c(log2ChIPNamesPlot))*3),
                                                       ((length(c(log2ChIPNamesPlot))*3)+1):(length(c(log2ChIPNamesPlot))*4),
                                                       ((length(c(log2ChIPNamesPlot))*4)+1):(length(c(log2ChIPNamesPlot))*5),
                                                       ((length(c(log2ChIPNamesPlot))*5)+1):(length(c(log2ChIPNamesPlot))*6)
                                                      ))
ggsave(paste0(plotDir,
              "log2ChIPcontrol_", align, "_avgProfiles_around_features_", quantiles, "quantiles",
              "_by_", context, "_absolute_change_in_", libName1[1], "_", libName2[1],
              "_at_", superfamName, "_", sortRegion,
              "_in_", refbase, "_",
              paste0(chrName, collapse = "_"), ".pdf"),
       plot = ggObjGA_combined,
       height = 6.5*length(c(log2ChIPNamesPlot)), width = 7*6, limitsize = FALSE)

#### Free up memory by removing no longer required objects
rm(
   log2ChIP_featureAcc1Mats, log2ChIP_ranLocAcc1Mats,
   log2ChIP_featureAcc2Mats, log2ChIP_ranLocAcc2Mats,
   log2ChIP_mats_quantiles,
   wideDFfeature_list_log2ChIP,
   tidyDFfeature_list_log2ChIP,
   summaryDFfeature_list_log2ChIP,
   summaryDFfeature_log2ChIP
  ) 
gc()
#####

