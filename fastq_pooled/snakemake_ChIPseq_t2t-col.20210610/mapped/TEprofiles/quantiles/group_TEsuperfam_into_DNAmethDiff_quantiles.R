#!/applications/R/R-4.0.0/bin/Rscript

# author: Andy Tock
# contact: ajt200@cam.ac.uk
# date: 27.08.2021

# Group TEs within a given superfamily into quantiles according to
# decreasing absolute change, log2 fold change, or relative change
# in context-specific DNA methylation levels in a given mutant genotype
# relative to wild type

# Usage:
# /applications/R/R-4.0.0/bin/Rscript group_TEsuperfam_into_DNAmethDiff_quantiles.R 'cmt3_BSseq_Rep1' 'WT_BSseq_Rep2_2013,WT_BSseq_Rep3_2013' 'BSseq_leaf_Stroud_Jacobsen_2013_Cell_2014_NSMB/snakemake_BSseq_t2t-col.20210610' 'BSseq_leaf_Stroud_Jacobsen_2013_Cell_2014_NSMB/snakemake_BSseq_t2t-col.20210610' 6 2000 2000 2kb '2 kb' 10 10bp t2t-col.20210610 'Chr1,Chr2,Chr3,Chr4,Chr5' bodies CHG Gypsy_LTR

#libName <- unlist(strsplit("cmt3_BSseq_Rep1",
#                           split = ","))
#controlName <- unlist(strsplit("WT_BSseq_Rep2_2013,WT_BSseq_Rep3_2013",
#                               split = ","))
#libDir <- "BSseq_leaf_Stroud_Jacobsen_2013_Cell_2014_NSMB/snakemake_BSseq_t2t-col.20210610"
#controlDir <- "BSseq_leaf_Stroud_Jacobsen_2013_Cell_2014_NSMB/snakemake_BSseq_t2t-col.20210610"
#quantiles <- 6
#regionBodyLength <- 2000
#upstream <- 2000
#downstream <- 2000
#flankName <- "2kb"
#flankNamePlot <- "2 kb"
#binSize <- 10
#binName <- "10bp"
#refbase <- "t2t-col.20210610"
#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split = ","))
#sortRegion <- "bodies"
#context <- "CHG"
#superfamName <- "Gypsy_LTR"

args <- commandArgs(trailingOnly = T)
libName <- unlist(strsplit(args[1],
                           split = ","))
controlName <- unlist(strsplit(args[2],
                               split = ","))
libDir <- args[3]
controlDir <- args[4]
quantiles <- as.integer(args[5])
regionBodyLength <- as.integer(args[6])
upstream <- as.integer(args[7])
downstream <- as.integer(args[7])
flankName <- args[8]
flankNamePlot <- args[9]
binSize <- as.integer(args[10])
binName <- args[11]
refbase <- args[12]
chrName <- unlist(strsplit(args[13],
                           split = ","))
sortRegion <- args[14]
context <- args[15]
superfamName <- args[16]

options(stringsAsFactors = F)
library(parallel)
library(dplyr)

# Genomic definitions
fai <- read.table(paste0("/home/ajt200/analysis/nanopore/", refbase, "/", refbase, ".fa.fai"), header = F)
chrs <- fai$V1[which(fai$V1 %in% chrName)]
chrLens <- fai$V2[which(fai$V1 %in% chrName)]

# Load table of representative feature coordinates in TSV and BED format
featuresTSV <- lapply(seq_along(chrName), function(x) {
  read.table(paste0("/home/ajt200/analysis/nanopore/", refbase,
                    "/annotation/TEs_EDTA/", refbase,
                    "_TEs_", superfamName, "_",
                    chrName[x], ".tsv"),
             header = T)
})
if(length(chrName) > 1) {
  featuresTSV <- do.call(rbind, featuresTSV)
} else {
  featuresTSV <- featuresTSV[[1]]
}

# Confirm BED files used for profile calculation have same
# row order as TSV file to be used for quantile definition
featuresBED <- lapply(seq_along(chrName), function(x) {
  read.table(paste0("/home/ajt200/analysis/nanopore/", refbase,
                    "/annotation/TEs_EDTA/", refbase,
                    "_TEs_", superfamName, "_",
                    chrName[x], ".bed"),
             header = F)
})
if(length(chrName) > 1) {
  featuresBED <- do.call(rbind, featuresBED)
} else {
  featuresBED <- featuresBED[[1]]
}

# Convert 0-based start coordinates (BED)
# into 1-based start coordinates
featuresBED[,2] <- as.integer(featuresBED[,2]+1)
colnames(featuresBED) <- c("chr", "start", "end", "name", "score", "strand")
stopifnot(identical(featuresTSV$Name, sub("\\d+_", "", featuresBED$name)))
stopifnot(identical(featuresTSV$start, featuresBED$start))
stopifnot(identical(featuresTSV$end, featuresBED$end))

rm(featuresBED); gc()

## Load table of coordinates for randomly positioned loci
ranLocBED <- lapply(seq_along(chrName), function(x) {
  read.table(paste0("/home/ajt200/analysis/nanopore/", refbase,
                    "/annotation/TEs_EDTA/", refbase,
                    "_TEs_", superfamName, "_",
                    chrName[x], ".bed"),
             header = F)
})
if(length(chrName) > 1) {
  ranLocBED <- do.call(rbind, ranLocBED)
} else {
  ranLocBED <- ranLocBED[[1]]
}

# Convert 0-based start coordinates (BED)
# into 1-based start coordinates
ranLocBED[,2] <- as.integer(ranLocBED[,2]+1)
colnames(ranLocBED) <- c("chr", "start", "end", "name", "score", "strand")

# Calculate context-specific DNA methylation differentials

#{sample}_MappedOn_{refbase}_{context}_TEs_{superfamName}_in_{chrName}_matrix_bin{binName}_flank{flankName}.tab

## lib
# feature

lib_featureMat <- lapply(seq_along(libName), function(w) {
  mclapply(seq_along(chrName), function(x) {
    as.matrix(read.table(paste0("/home/ajt200/analysis/",
                                libDir,
                                "/coverage/TEprofiles/matrices/",
                                libName[w],
                                "_MappedOn_", refbase, "_", context,
                                "_TEs_", superfamName, "_in_", chrName[x],
                                "_matrix_bin", binName, "_flank", flankName, ".tab"),
                         header = F, skip = 3))
  }, mc.cores = length(chrName), mc.preschedule = F)
})

control_featureMat <- lapply(seq_along(controlName), function(w) {
  mclapply(seq_along(chrName), function(x) {
    as.matrix(read.table(paste0("/home/ajt200/analysis/",
                                controlDir,
                                "/coverage/TEprofiles/matrices/",
                                controlName[w],
                                "_MappedOn_", refbase, "_", context,
                                "_TEs_", superfamName, "_in_", chrName[x],
                                "_matrix_bin", binName, "_flank", flankName, ".tab"),
                         header = F, skip = 3))
  }, mc.cores = length(chrName), mc.preschedule = F)
})

# Concatenate feature rows from successive chromosomes
for(w in 1:length(control_featureMat)) {
  if(length(chrName) > 1) {
    control_featureMat[[w]] <- do.call(rbind, control_featureMat[[w]])
  } else {
    control_featureMat[[w]] <- control_featureMat[[w]][[1]]
  }
}
# Calculate average values across replicates
if(length(control_featureMat) > 1) {
  control_featureMat <- 
  
if(length(chrName) > 1) {
  featuresBED <- do.call(rbind, featuresBED)
} else {
  featuresBED <- featuresBED[[1]]
}

## control
# feature
control_featureMat <- as.matrix(read.table(paste0("/home/ajt200/analysis/",
                                           controlDir,
                                           "/snakemake_ChIPseq/mapped/geneProfiles/matrices/",
                                           controlName,
                                           "_MappedOn_TAIR10_chr_all_lowXM_both_sort_norm_genes",
                                           "_matrix_bin", binName, "_flank", flankName, ".tab"),
                                        header = F, skip = 3))

# Calculate log2(ChIP/control) coverage values
log2ChIP_featureMat <- log2((ChIP_featureMat+1)/(control_featureMat+1))

# Calculate mean log2(ChIP/control) coverage values for each feature
if(sortRegion == "promoters") {
  log2ChIP_featureMat_sortRegion <- log2ChIP_featureMat[,(((upstream-1000)/binSize)+1):(upstream/binSize)]
} else if(sortRegion == "bodies") {
  log2ChIP_featureMat_sortRegion <- log2ChIP_featureMat[,((upstream/binSize)+1):((upstream+regionBodyLength)/binSize)]
} else if(sortRegion == "terminators") {
  log2ChIP_featureMat_sortRegion <- log2ChIP_featureMat[,(((upstream+regionBodyLength)/binSize)+1):(((upstream+regionBodyLength)/binSize)+(1000/binSize))]
} else if(sortRegion == "genes") {
  log2ChIP_featureMat_sortRegion <- log2ChIP_featureMat[,(((upstream-1000)/binSize)+1):(((upstream+regionBodyLength)/binSize)+(1000/binSize))]
} else {
  stop("sortRegion is none of promoters, bodies, terminators or genes")
}

log2ChIP_featureMat_sortRegionRowMeans <- rowMeans(log2ChIP_featureMat_sortRegion, na.rm = T)

# Add mean coverage values to features dataframe
features <- data.frame(features,
                       meanCov = log2ChIP_featureMat_sortRegionRowMeans)

# Group features into quantiles according to
# decreasing log2ChIP_featureMat_sortRegionRowMeans
outDir <- paste0("quantiles_by_", libName, "_in_", sortRegion, "/")
plotDir <- paste0(outDir, "/plots/")
system(paste0("[ -d ", outDir, " ] || mkdir -p ", outDir))
system(paste0("[ -d ", plotDir, " ] || mkdir -p ", plotDir))

features_DF <- data.frame(features,
                          quantile = as.character(""))
# Assign 0s to NA values only for coverage data
features_DF[,which(colnames(features_DF) == "meanCov")][
    which(is.na(features_DF[,which(colnames(features_DF) == "meanCov")]))] <- 0
quantilesStats <- data.frame()
for(k in 1:quantiles) {
  # First quantile should span 1 to greater than, e.g., 0.75
  if(k < quantiles) {
    features_DF[ !is.na(features_DF[,which(colnames(features_DF) == "meanCov")]) &
                 rank(features_DF[,which(colnames(features_DF) == "meanCov")]) /
                 length(features_DF[,which(colnames(features_DF) == "meanCov")]) <=
                 1-((k-1)/quantiles) &
                 rank(features_DF[,which(colnames(features_DF) == "meanCov")]) /
                 length(features_DF[,which(colnames(features_DF) == "meanCov")]) >
                 1-(k/quantiles), ]$quantile <- paste0("Quantile ", k)
  } else {
  # Final quantile should span 0 to, e.g., 0.25
    features_DF[ !is.na(features_DF[,which(colnames(features_DF) == "meanCov")]) &
                 rank(features_DF[,which(colnames(features_DF) == "meanCov")]) /
                 length(features_DF[,which(colnames(features_DF) == "meanCov")]) <=
                 1-((k-1)/quantiles) &
                 rank(features_DF[,which(colnames(features_DF) == "meanCov")]) /
                 length(features_DF[,which(colnames(features_DF) == "meanCov")]) >=
                 1-(k/quantiles), ]$quantile <- paste0("Quantile ", k)
  }
  write.table(features_DF[features_DF$quantile == paste0("Quantile ", k),],
              file = paste0(outDir,
                            "quantile", k, "_of_", quantiles,
                            "_by_log2_", libName, "_control",
                            "_in_", sortRegion, "_of_Araport11_genes.tsv"),
              quote = FALSE, sep = "\t", row.names = FALSE)
  stats <- data.frame(quantile = as.integer(k),
                      n = as.integer(dim(features_DF[features_DF$quantile == paste0("Quantile ", k),])[1]),
                      mean_width = as.integer(round(mean(
                        (features_DF[features_DF$quantile == paste0("Quantile ", k),]$end -
                         features_DF[features_DF$quantile == paste0("Quantile ", k),]$start) + 1, na.rm = T))),
                      total_width = as.integer(sum(
                        (features_DF[features_DF$quantile == paste0("Quantile ", k),]$end -
                         features_DF[features_DF$quantile == paste0("Quantile ", k),]$start) + 1, na.rm = T)),
                      mean_orderingFactor = as.numeric(mean(features_DF[features_DF$quantile == paste0("Quantile ", k),][,which(colnames(features_DF) == "meanCov")], na.rm = T)))
  quantilesStats <- rbind(quantilesStats, stats)
}
write.table(quantilesStats,
            file = paste0(outDir,
                          "summary_", quantiles, "quantiles",
                          "_by_log2_", libName, "_control",
                          "_in_", sortRegion, "_of_Araport11_genes.tsv"),
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(features_DF,
            file = paste0(outDir,
                          "features_", quantiles, "quantiles",
                          "_by_log2_", libName, "_control",
                          "_in_", sortRegion, "_of_Araport11_genes.tsv"),
            quote = FALSE, sep = "\t", row.names = FALSE)

# Divide ranLoc into quantiles based on feature quantile indices
ranLoc_DF <- data.frame(ranLoc,
                        random = as.character(""))
# Get row indices for each feature quantile
quantileIndices <- lapply(1:quantiles, function(k) {
  which(features_DF$quantile == paste0("Quantile ", k))
})
for(k in 1:quantiles) {
  ranLoc_DF[quantileIndices[[k]],]$random <- paste0("Random ", k)
}
write.table(ranLoc_DF,
            file = paste0(outDir,
                          "features_", quantiles, "quantiles",
                          "_by_log2_", libName, "_control",
                          "_in_", sortRegion, "_of_Araport11_genes_ranLoc.tsv"),
            quote = FALSE, sep = "\t", row.names = FALSE)
