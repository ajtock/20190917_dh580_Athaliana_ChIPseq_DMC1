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

outDir <- paste0(paste0(chrName, collapse = "_"), "/")
system(paste0("[ -d ", outDir, " ] || mkdir -p ", outDir))

# Load table of representative feature coordinates in TSV and BED format
featureTSV <- lapply(seq_along(chrName), function(x) {
  read.table(paste0("/home/ajt200/analysis/nanopore/", refbase,
                    "/annotation/TEs_EDTA/", refbase,
                    "_TEs_", superfamName, "_",
                    chrName[x], ".tsv"),
             header = T)
})
if(length(chrName) > 1) {
  featureTSV <- do.call(rbind, featureTSV)
} else {
  featureTSV <- featureTSV[[1]]
}

# Confirm BED files used for profile calculation have same
# row order as TSV file to be used for quantile definition
featureBED <- lapply(seq_along(chrName), function(x) {
  read.table(paste0("/home/ajt200/analysis/nanopore/", refbase,
                    "/annotation/TEs_EDTA/", refbase,
                    "_TEs_", superfamName, "_",
                    chrName[x], ".bed"),
             header = F)
})
if(length(chrName) > 1) {
  featureBED <- do.call(rbind, featureBED)
} else {
  featureBED <- featureBED[[1]]
}

# Convert 0-based start coordinates (BED)
# into 1-based start coordinates
featureBED[,2] <- as.integer(featureBED[,2]+1)
colnames(featureBED) <- c("chr", "start", "end", "name", "score", "strand")
stopifnot(identical(featureTSV$Name, sub("\\d+_", "", featureBED$name)))
stopifnot(identical(featureTSV$start, featureBED$start))
stopifnot(identical(featureTSV$end, featureBED$end))

rm(featureBED); gc()

feature <- featureTSV
rm(featureTSV); gc()

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

ranLoc <- ranLocBED
rm(ranLocBED); gc()

# Calculate context-specific DNA methylation differentials
# I.e., absolute change, log2 fold change, or relative change
# in context-specific DNA methylation levels in a given mutant genotype
# relative to wild type

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

# Concatenate feature rows from successive chromosomes
for(w in 1:length(lib_featureMat)) {
  if(length(chrName) > 1) {
    lib_featureMat[[w]] <- do.call(rbind, lib_featureMat[[w]])
  } else {
    lib_featureMat[[w]] <- lib_featureMat[[w]][[1]]
  }
}
# Calculate mean values across replicates
if(length(lib_featureMat) > 1) {
  lib_featureMat_cbind <- do.call(cbind, lib_featureMat)
  lib_featureMat_array <- array(data = lib_featureMat_cbind,
                                dim = c(dim(lib_featureMat[[1]]), length(lib_featureMat)))
  lib_featureMat <- colMeans(aperm(a = lib_featureMat_array,
                                   perm = c(3, 1, 2)),
                             na.rm = T)
} else {
  lib_featureMat <- lib_featureMat[[1]]
}

stopifnot(identical(nrow(feature), nrow(lib_featureMat)))

# ranLoc
lib_ranLocMat <- lapply(seq_along(libName), function(w) {
  mclapply(seq_along(chrName), function(x) {
    as.matrix(read.table(paste0("/home/ajt200/analysis/",
                                libDir,
                                "/coverage/TEprofiles/matrices/",
                                libName[w],
                                "_MappedOn_", refbase, "_", context,
                                "_TEs_", superfamName, "_in_", chrName[x],
                                "_ranLoc_matrix_bin", binName, "_flank", flankName, ".tab"),
                         header = F, skip = 3))
  }, mc.cores = length(chrName), mc.preschedule = F)
})

# Concatenate ranLoc rows from successive chromosomes
for(w in 1:length(lib_ranLocMat)) {
  if(length(chrName) > 1) {
    lib_ranLocMat[[w]] <- do.call(rbind, lib_ranLocMat[[w]])
  } else {
    lib_ranLocMat[[w]] <- lib_ranLocMat[[w]][[1]]
  }
}
# Calculate mean values across replicates
if(length(lib_ranLocMat) > 1) {
  lib_ranLocMat_cbind <- do.call(cbind, lib_ranLocMat)
  lib_ranLocMat_array <- array(data = lib_ranLocMat_cbind,
                               dim = c(dim(lib_ranLocMat[[1]]), length(lib_ranLocMat)))
  lib_ranLocMat <- colMeans(aperm(a = lib_ranLocMat_array,
                                  perm = c(3, 1, 2)),
                            na.rm = T)
} else {
  lib_ranLocMat <- lib_ranLocMat[[1]]
}

stopifnot(identical(nrow(ranLoc), nrow(lib_ranLocMat)))

## control
# feature
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
# Calculate mean values across replicates
if(length(control_featureMat) > 1) {
  control_featureMat_cbind <- do.call(cbind, control_featureMat)
  control_featureMat_array <- array(data = control_featureMat_cbind,
                                    dim = c(dim(control_featureMat[[1]]), length(control_featureMat)))
  control_featureMat <- colMeans(aperm(a = control_featureMat_array,
                                       perm = c(3, 1, 2)),
                                 na.rm = T)
} else {
  control_featureMat <- control_featureMat[[1]]
}

stopifnot(identical(nrow(feature), nrow(control_featureMat)))

# ranLoc
control_ranLocMat <- lapply(seq_along(controlName), function(w) {
  mclapply(seq_along(chrName), function(x) {
    as.matrix(read.table(paste0("/home/ajt200/analysis/",
                                controlDir,
                                "/coverage/TEprofiles/matrices/",
                                controlName[w],
                                "_MappedOn_", refbase, "_", context,
                                "_TEs_", superfamName, "_in_", chrName[x],
                                "_ranLoc_matrix_bin", binName, "_flank", flankName, ".tab"),
                         header = F, skip = 3))
  }, mc.cores = length(chrName), mc.preschedule = F)
})

# Concatenate ranLoc rows from successive chromosomes
for(w in 1:length(control_ranLocMat)) {
  if(length(chrName) > 1) {
    control_ranLocMat[[w]] <- do.call(rbind, control_ranLocMat[[w]])
  } else {
    control_ranLocMat[[w]] <- control_ranLocMat[[w]][[1]]
  }
}
# Calculate mean values across replicates
if(length(control_ranLocMat) > 1) {
  control_ranLocMat_cbind <- do.call(cbind, control_ranLocMat)
  control_ranLocMat_array <- array(data = control_ranLocMat_cbind,
                                   dim = c(dim(control_ranLocMat[[1]]), length(control_ranLocMat)))
  control_ranLocMat <- colMeans(aperm(a = control_ranLocMat_array,
                                      perm = c(3, 1, 2)),
                                na.rm = T)
} else {
  control_ranLocMat <- control_ranLocMat[[1]]
}

stopifnot(identical(nrow(ranLoc), nrow(control_ranLocMat)))


# Get sortRegion
if(sortRegion == "promoters") {
  lib_featureMat_sortRegion <- lib_featureMat[,(((upstream-1000)/binSize)+1):(upstream/binSize)]
  control_featureMat_sortRegion <- control_featureMat[,(((upstream-1000)/binSize)+1):(upstream/binSize)]

  lib_ranLocMat_sortRegion <- lib_ranLocMat[,(((upstream-1000)/binSize)+1):(upstream/binSize)]
  control_ranLocMat_sortRegion <- control_ranLocMat[,(((upstream-1000)/binSize)+1):(upstream/binSize)]
} else if(sortRegion == "bodies") {
  lib_featureMat_sortRegion <- lib_featureMat[,((upstream/binSize)+1):((upstream+regionBodyLength)/binSize)]
  control_featureMat_sortRegion <- control_featureMat[,((upstream/binSize)+1):((upstream+regionBodyLength)/binSize)]

  lib_ranLocMat_sortRegion <- lib_ranLocMat[,((upstream/binSize)+1):((upstream+regionBodyLength)/binSize)]
  control_ranLocMat_sortRegion <- control_ranLocMat[,((upstream/binSize)+1):((upstream+regionBodyLength)/binSize)]
} else if(sortRegion == "terminators") {
  lib_featureMat_sortRegion <- lib_featureMat[,(((upstream+regionBodyLength)/binSize)+1):(((upstream+regionBodyLength)/binSize)+(1000/binSize))]
  control_featureMat_sortRegion <- control_featureMat[,(((upstream+regionBodyLength)/binSize)+1):(((upstream+regionBodyLength)/binSize)+(1000/binSize))]

  lib_ranLocMat_sortRegion <- lib_ranLocMat[,(((upstream+regionBodyLength)/binSize)+1):(((upstream+regionBodyLength)/binSize)+(1000/binSize))]
  control_ranLocMat_sortRegion <- control_ranLocMat[,(((upstream+regionBodyLength)/binSize)+1):(((upstream+regionBodyLength)/binSize)+(1000/binSize))]
} else if(sortRegion == "TEs") {
  lib_featureMat_sortRegion <- lib_featureMat[,(((upstream-1000)/binSize)+1):(((upstream+regionBodyLength)/binSize)+(1000/binSize))]
  control_featureMat_sortRegion <- control_featureMat[,(((upstream-1000)/binSize)+1):(((upstream+regionBodyLength)/binSize)+(1000/binSize))]

  lib_ranLocMat_sortRegion <- lib_ranLocMat[,(((upstream-1000)/binSize)+1):(((upstream+regionBodyLength)/binSize)+(1000/binSize))]
  control_ranLocMat_sortRegion <- control_ranLocMat[,(((upstream-1000)/binSize)+1):(((upstream+regionBodyLength)/binSize)+(1000/binSize))]
} else {
  stop("sortRegion is none of promoters, bodies, terminators or TEs")
}

lib_featureMat_proportion <- rowMeans(lib_featureMat_sortRegion, na.rm = T)
control_featureMat_proportion <- rowMeans(control_featureMat_sortRegion, na.rm = T)
lib_ranLocMat_proportion <- rowMeans(lib_ranLocMat_sortRegion, na.rm = T)
control_ranLocMat_proportion <- rowMeans(control_ranLocMat_sortRegion, na.rm = T)

# feature
# Create new columns containing absolute change, log2 fold change, and relative change in context methylation proportion
feature$absolute_change <- as.numeric( control_featureMat_proportion - lib_featureMat_proportion )
# + 0.01 is an offset to prevent infinite log2 fold change (or equal relative change) values where proportions = 0
feature$log2_fold_change <- as.numeric( log2( (control_featureMat_proportion + 0.01) /
                                              (lib_featureMat_proportion + 0.01) ) )
feature$relative_change <- as.numeric( 1 - ( (lib_featureMat_proportion + 0.01) /
                                             (control_featureMat_proportion + 0.01) ) )

# Create new columns containing absolute change, log2 fold change, and relative change percentiles and quantiles
feature$ac_percentile <- as.numeric( rank(feature$absolute_change, na.last = "keep") /
                                     length(feature$absolute_change) )
feature$l2fc_percentile <- as.numeric( rank(feature$log2_fold_change, na.last = "keep") /
                                       length(feature$log2_fold_change) )
feature$rc_percentile <- as.numeric( rank(feature$relative_change, na.last = "keep") /
                                     length(feature$relative_change) )
feature$ac_quantile <- as.character("")
feature$l2fc_quantile <- as.character("")
feature$rc_quantile <- as.character("")

# Define quantiles
ac_quantilesStats <- data.frame()
l2fc_quantilesStats <- data.frame()
rc_quantilesStats <- data.frame()
for(k in 1:quantiles) {
  # absolute_change (ac)
  # First quantile should span 1 to greater than, e.g., 0.75 proportions of feature
  if(k < quantiles) {
     feature[ !is.na(feature$ac_percentile) &
                     feature$ac_percentile <= 1 - ( (k - 1) / quantiles ) &
                     feature$ac_percentile > 1 - ( k / quantiles ), ]$ac_quantile <- paste0("Quantile ", k)
  } else {
  # Final quantile should span 0 to, e.g., 0.25 proportions of feature
    feature[ !is.na(feature$ac_percentile) &
                    feature$ac_percentile <= 1 - ( (k - 1) / quantiles ) &
                    feature$ac_percentile >= 1 - ( k / quantiles ), ]$ac_quantile <- paste0("Quantile ", k)
  }
  ac_stats <- data.frame(quantile = as.integer(k),
                         n = as.integer(nrow(feature[feature$ac_quantile == paste0("Quantile ", k),])),
                         mean_width = as.integer(round(mean(
                           feature[feature$ac_quantile == paste0("Quantile ", k),]$end-feature[feature$ac_quantile == paste0("Quantile ", k),]$start+1, na.rm = T))),
                         total_width = as.integer(sum(
                           feature[feature$ac_quantile == paste0("Quantile ", k),]$end-feature[feature$ac_quantile == paste0("Quantile ", k),]$start+1, na.rm = T)),
                         mean_absolute_change = as.numeric(mean(
                           feature[feature$ac_quantile == paste0("Quantile ", k),]$absolute_change, na.rm = T)),
                         mean_log2_fold_change = as.numeric(mean(
                           feature[feature$ac_quantile == paste0("Quantile ", k),]$log2_fold_change, na.rm = T)),
                         mean_relative_change = as.numeric(mean(
                           feature[feature$ac_quantile == paste0("Quantile ", k),]$relative_change, na.rm = T)),
                         stringsAsFactors = FALSE)
  ac_quantilesStats <- rbind(ac_quantilesStats, ac_stats)
  # log2_fold_change
  # First quantile should span 1 to greater than, e.g., 0.75 proportions of feature
  if(k < quantiles) {
     feature[ !is.na(feature$l2fc_percentile) &
                     feature$l2fc_percentile <= 1 - ( (k - 1) / quantiles ) &
                     feature$l2fc_percentile > 1 - ( k / quantiles ), ]$l2fc_quantile <- paste0("Quantile ", k)
  } else {
  # Final quantile should span 0 to, e.g., 0.25 proportions of feature
    feature[ !is.na(feature$l2fc_percentile) &
                    feature$l2fc_percentile <= 1 - ( (k - 1) / quantiles ) &
                    feature$l2fc_percentile >= 1 - ( k / quantiles ), ]$l2fc_quantile <- paste0("Quantile ", k)
  }
  l2fc_stats <- data.frame(quantile = as.integer(k),
                           n = as.integer(nrow(feature[feature$l2fc_quantile == paste0("Quantile ", k),])),
                           mean_width = as.integer(round(mean(
                             feature[feature$l2fc_quantile == paste0("Quantile ", k),]$end-feature[feature$l2fc_quantile == paste0("Quantile ", k),]$start+1, na.rm = T))),
                           total_width = as.integer(sum(
                             feature[feature$l2fc_quantile == paste0("Quantile ", k),]$end-feature[feature$l2fc_quantile == paste0("Quantile ", k),]$start+1, na.rm = T)),
                           mean_absolute_change = as.numeric(mean(
                             feature[feature$l2fc_quantile == paste0("Quantile ", k),]$absolute_change, na.rm = T)),
                           mean_log2_fold_change = as.numeric(mean(
                             feature[feature$l2fc_quantile == paste0("Quantile ", k),]$log2_fold_change, na.rm = T)),
                           mean_relative_change = as.numeric(mean(
                             feature[feature$l2fc_quantile == paste0("Quantile ", k),]$relative_change, na.rm = T)),
                           stringsAsFactors = FALSE)
  l2fc_quantilesStats <- rbind(l2fc_quantilesStats, l2fc_stats)
  # relative_change
  # First quantile should span 1 to greater than, e.g., 0.75 proportions of feature
  if(k < quantiles) {
     feature[ !is.na(feature$rc_percentile) &
                     feature$rc_percentile <= 1 - ( (k - 1) / quantiles ) &
                     feature$rc_percentile > 1 - ( k / quantiles ), ]$rc_quantile <- paste0("Quantile ", k)
  } else {
  # Final quantile should span 0 to, e.g., 0.25 proportions of feature
    feature[ !is.na(feature$rc_percentile) &
                    feature$rc_percentile <= 1 - ( (k - 1) / quantiles ) &
                    feature$rc_percentile >= 1 - ( k / quantiles ), ]$rc_quantile <- paste0("Quantile ", k)
  }
  rc_stats <- data.frame(quantile = as.integer(k),
                         n = as.integer(nrow(feature[feature$rc_quantile == paste0("Quantile ", k),])),
                         mean_width = as.integer(round(mean(
                           feature[feature$rc_quantile == paste0("Quantile ", k),]$end-feature[feature$rc_quantile == paste0("Quantile ", k),]$start+1, na.rm = T))),
                         total_width = as.integer(sum(
                           feature[feature$rc_quantile == paste0("Quantile ", k),]$end-feature[feature$rc_quantile == paste0("Quantile ", k),]$start+1, na.rm = T)),
                         mean_absolute_change = as.numeric(mean(
                           feature[feature$rc_quantile == paste0("Quantile ", k),]$absolute_change, na.rm = T)),
                         mean_log2_fold_change = as.numeric(mean(
                           feature[feature$rc_quantile == paste0("Quantile ", k),]$log2_fold_change, na.rm = T)),
                         mean_relative_change = as.numeric(mean(
                           feature[feature$rc_quantile == paste0("Quantile ", k),]$relative_change, na.rm = T)),
                         stringsAsFactors = FALSE)
  rc_quantilesStats <- rbind(rc_quantilesStats, rc_stats)
}
write.table(ac_quantilesStats,
            file = paste0(outDir,
                          "feature_summary_", quantiles, "quantiles",
                          "_by_", sub("p", "", context), "_absolute_change_in_",
                          paste0(libName, collapse = "_"), "_vs_", paste0(controlName, collapse = "_"), 
                          "_in_", sortRegion, "_of_", superfamName, "_in_", refbase, "_", paste0(chrName, collapse = "_"), ".tsv"),
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(l2fc_quantilesStats,
            file = paste0(outDir,
                          "feature_summary_", quantiles, "quantiles",
                          "_by_", sub("p", "", context), "_log2_fold_change_in_",
                          paste0(libName, collapse = "_"), "_vs_", paste0(controlName, collapse = "_"), 
                          "_in_", sortRegion, "_of_", superfamName, "_in_", refbase, "_", paste0(chrName, collapse = "_"), ".tsv"),
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(rc_quantilesStats,
            file = paste0(outDir,
                          "feature_summary_", quantiles, "quantiles",
                          "_by_", sub("p", "", context), "_relative_change_in_",
                          paste0(libName, collapse = "_"), "_vs_", paste0(controlName, collapse = "_"), 
                          "_in_", sortRegion, "_of_", superfamName, "_in_", refbase, "_", paste0(chrName, collapse = "_"), ".tsv"),
            quote = FALSE, sep = "\t", row.names = FALSE)

# Export feature as annotation files
write.table(feature,
            file = paste0(outDir,
                          "feature_", quantiles, "quantiles",
                          "_by_", sub("p", "", context), "_change_in_",
                          paste0(libName, collapse = "_"), "_vs_", paste0(controlName, collapse = "_"), 
                          "_in_", sortRegion, "_of_", superfamName, "_in_", refbase, "_", paste0(chrName, collapse = "_"), ".tsv"),
            quote = FALSE, sep = "\t", row.names = FALSE)
feature_BED <- data.frame(chr = as.character(feature$seqid),
                          start = as.integer(feature$start-1),
                          end = as.integer(feature$end),
                          name = as.character(feature$Name),
                          score = as.numeric(feature$log2_fold_change),
                          strand = as.character(feature$strand),
                          stringsAsFactors = FALSE)
write.table(feature_BED,
            file = paste0(outDir,
                          "feature_", quantiles, "quantiles",
                          "_by_", sub("p", "", context), "_change_in_",
                          paste0(libName, collapse = "_"), "_vs_", paste0(controlName, collapse = "_"), 
                          "_in_", sortRegion, "_of_", superfamName, "_in_", refbase, "_", paste0(chrName, collapse = "_"), ".bed"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


# ranLoc
# Create new columns containing absolute change, log2 fold change, and relative change in context methylation proportion
ranLoc$absolute_change <- as.numeric( control_ranLocMat_proportion - lib_ranLocMat_proportion )
# + 0.01 is an offset to prevent infinite log2 fold change (or equal relative change) values where proportions = 0
ranLoc$log2_fold_change <- as.numeric( log2( (control_ranLocMat_proportion + 0.01) /
                                             (lib_ranLocMat_proportion + 0.01) ) )
ranLoc$relative_change <- as.numeric( 1 - ( (lib_ranLocMat_proportion + 0.01) /
                                            (control_ranLocMat_proportion + 0.01) ) )

# Create new columns containing absolute change, log2 fold change, and relative change percentiles and quantiles
ranLoc$ac_percentile <- as.numeric( rank(ranLoc$absolute_change, na.last = "keep") /
                                    length(ranLoc$absolute_change) )
ranLoc$l2fc_percentile <- as.numeric( rank(ranLoc$log2_fold_change, na.last = "keep") /
                                      length(ranLoc$log2_fold_change) )
ranLoc$rc_percentile <- as.numeric( rank(ranLoc$relative_change, na.last = "keep") /
                                    length(ranLoc$relative_change) )
ranLoc$ac_quantile <- as.character("")
ranLoc$l2fc_quantile <- as.character("")
ranLoc$rc_quantile <- as.character("")

# Define quantiles
ac_quantilesStats <- data.frame()
l2fc_quantilesStats <- data.frame()
rc_quantilesStats <- data.frame()
for(k in 1:quantiles) {
  # absolute_change (ac)
  # First quantile should span 1 to greater than, e.g., 0.75 proportions of ranLoc
  if(k < quantiles) {
     ranLoc[ !is.na(ranLoc$ac_percentile) &
                    ranLoc$ac_percentile <= 1 - ( (k - 1) / quantiles ) &
                    ranLoc$ac_percentile > 1 - ( k / quantiles ), ]$ac_quantile <- paste0("Quantile ", k)
  } else {
  # Final quantile should span 0 to, e.g., 0.25 proportions of ranLoc
    ranLoc[ !is.na(ranLoc$ac_percentile) &
                   ranLoc$ac_percentile <= 1 - ( (k - 1) / quantiles ) &
                   ranLoc$ac_percentile >= 1 - ( k / quantiles ), ]$ac_quantile <- paste0("Quantile ", k)
  }
  ac_stats <- data.frame(quantile = as.integer(k),
                         n = as.integer(nrow(ranLoc[ranLoc$ac_quantile == paste0("Quantile ", k),])),
                         mean_width = as.integer(round(mean(
                           ranLoc[ranLoc$ac_quantile == paste0("Quantile ", k),]$end-ranLoc[ranLoc$ac_quantile == paste0("Quantile ", k),]$start+1, na.rm = T))),
                         total_width = as.integer(sum(
                           ranLoc[ranLoc$ac_quantile == paste0("Quantile ", k),]$end-ranLoc[ranLoc$ac_quantile == paste0("Quantile ", k),]$start+1, na.rm = T)),
                         mean_absolute_change = as.numeric(mean(
                           ranLoc[ranLoc$ac_quantile == paste0("Quantile ", k),]$absolute_change, na.rm = T)),
                         mean_log2_fold_change = as.numeric(mean(
                           ranLoc[ranLoc$ac_quantile == paste0("Quantile ", k),]$log2_fold_change, na.rm = T)),
                         mean_relative_change = as.numeric(mean(
                           ranLoc[ranLoc$ac_quantile == paste0("Quantile ", k),]$relative_change, na.rm = T)),
                         stringsAsFactors = FALSE)
  ac_quantilesStats <- rbind(ac_quantilesStats, ac_stats)
  # log2_fold_change
  # First quantile should span 1 to greater than, e.g., 0.75 proportions of ranLoc
  if(k < quantiles) {
     ranLoc[ !is.na(ranLoc$l2fc_percentile) &
                    ranLoc$l2fc_percentile <= 1 - ( (k - 1) / quantiles ) &
                    ranLoc$l2fc_percentile > 1 - ( k / quantiles ), ]$l2fc_quantile <- paste0("Quantile ", k)
  } else {
  # Final quantile should span 0 to, e.g., 0.25 proportions of ranLoc
    ranLoc[ !is.na(ranLoc$l2fc_percentile) &
                   ranLoc$l2fc_percentile <= 1 - ( (k - 1) / quantiles ) &
                   ranLoc$l2fc_percentile >= 1 - ( k / quantiles ), ]$l2fc_quantile <- paste0("Quantile ", k)
  }
  l2fc_stats <- data.frame(quantile = as.integer(k),
                           n = as.integer(nrow(ranLoc[ranLoc$l2fc_quantile == paste0("Quantile ", k),])),
                           mean_width = as.integer(round(mean(
                             ranLoc[ranLoc$l2fc_quantile == paste0("Quantile ", k),]$end-ranLoc[ranLoc$l2fc_quantile == paste0("Quantile ", k),]$start+1, na.rm = T))),
                           total_width = as.integer(sum(
                             ranLoc[ranLoc$l2fc_quantile == paste0("Quantile ", k),]$end-ranLoc[ranLoc$l2fc_quantile == paste0("Quantile ", k),]$start+1, na.rm = T)),
                           mean_absolute_change = as.numeric(mean(
                             ranLoc[ranLoc$l2fc_quantile == paste0("Quantile ", k),]$absolute_change, na.rm = T)),
                           mean_log2_fold_change = as.numeric(mean(
                             ranLoc[ranLoc$l2fc_quantile == paste0("Quantile ", k),]$log2_fold_change, na.rm = T)),
                           mean_relative_change = as.numeric(mean(
                             ranLoc[ranLoc$l2fc_quantile == paste0("Quantile ", k),]$relative_change, na.rm = T)),
                           stringsAsFactors = FALSE)
  l2fc_quantilesStats <- rbind(l2fc_quantilesStats, l2fc_stats)
  # relative_change
  # First quantile should span 1 to greater than, e.g., 0.75 proportions of ranLoc
  if(k < quantiles) {
     ranLoc[ !is.na(ranLoc$rc_percentile) &
                    ranLoc$rc_percentile <= 1 - ( (k - 1) / quantiles ) &
                    ranLoc$rc_percentile > 1 - ( k / quantiles ), ]$rc_quantile <- paste0("Quantile ", k)
  } else {
  # Final quantile should span 0 to, e.g., 0.25 proportions of ranLoc
    ranLoc[ !is.na(ranLoc$rc_percentile) &
                   ranLoc$rc_percentile <= 1 - ( (k - 1) / quantiles ) &
                   ranLoc$rc_percentile >= 1 - ( k / quantiles ), ]$rc_quantile <- paste0("Quantile ", k)
  }
  rc_stats <- data.frame(quantile = as.integer(k),
                         n = as.integer(nrow(ranLoc[ranLoc$rc_quantile == paste0("Quantile ", k),])),
                         mean_width = as.integer(round(mean(
                           ranLoc[ranLoc$rc_quantile == paste0("Quantile ", k),]$end-ranLoc[ranLoc$rc_quantile == paste0("Quantile ", k),]$start+1, na.rm = T))),
                         total_width = as.integer(sum(
                           ranLoc[ranLoc$rc_quantile == paste0("Quantile ", k),]$end-ranLoc[ranLoc$rc_quantile == paste0("Quantile ", k),]$start+1, na.rm = T)),
                         mean_absolute_change = as.numeric(mean(
                           ranLoc[ranLoc$rc_quantile == paste0("Quantile ", k),]$absolute_change, na.rm = T)),
                         mean_log2_fold_change = as.numeric(mean(
                           ranLoc[ranLoc$rc_quantile == paste0("Quantile ", k),]$log2_fold_change, na.rm = T)),
                         mean_relative_change = as.numeric(mean(
                           ranLoc[ranLoc$rc_quantile == paste0("Quantile ", k),]$relative_change, na.rm = T)),
                         stringsAsFactors = FALSE)
  rc_quantilesStats <- rbind(rc_quantilesStats, rc_stats)
}
write.table(ac_quantilesStats,
            file = paste0(outDir,
                          "ranLoc_summary_", quantiles, "quantiles",
                          "_by_", sub("p", "", context), "_absolute_change_in_",
                          paste0(libName, collapse = "_"), "_vs_", paste0(controlName, collapse = "_"), 
                          "_in_", sortRegion, "_of_", superfamName, "_in_", refbase, "_", paste0(chrName, collapse = "_"), ".tsv"),
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(l2fc_quantilesStats,
            file = paste0(outDir,
                          "ranLoc_summary_", quantiles, "quantiles",
                          "_by_", sub("p", "", context), "_log2_fold_change_in_",
                          paste0(libName, collapse = "_"), "_vs_", paste0(controlName, collapse = "_"), 
                          "_in_", sortRegion, "_of_", superfamName, "_in_", refbase, "_", paste0(chrName, collapse = "_"), ".tsv"),
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(rc_quantilesStats,
            file = paste0(outDir,
                          "ranLoc_summary_", quantiles, "quantiles",
                          "_by_", sub("p", "", context), "_relative_change_in_",
                          paste0(libName, collapse = "_"), "_vs_", paste0(controlName, collapse = "_"), 
                          "_in_", sortRegion, "_of_", superfamName, "_in_", refbase, "_", paste0(chrName, collapse = "_"), ".tsv"),
            quote = FALSE, sep = "\t", row.names = FALSE)

# Export ranLoc as annotation files
write.table(ranLoc,
            file = paste0(outDir,
                          "ranLoc_", quantiles, "quantiles",
                          "_by_", sub("p", "", context), "_change_in_",
                          paste0(libName, collapse = "_"), "_vs_", paste0(controlName, collapse = "_"), 
                          "_in_", sortRegion, "_of_", superfamName, "_in_", refbase, "_", paste0(chrName, collapse = "_"), ".tsv"),
            quote = FALSE, sep = "\t", row.names = FALSE)
ranLoc_BED <- data.frame(chr = as.character(ranLoc$chr),
                         start = as.integer(ranLoc$start-1),
                         end = as.integer(ranLoc$end),
                         name = as.character(ranLoc$name),
                         score = as.numeric(ranLoc$log2_fold_change),
                         strand = as.character(ranLoc$strand),
                         stringsAsFactors = FALSE)
write.table(ranLoc_BED,
            file = paste0(outDir,
                          "ranLoc_", quantiles, "quantiles",
                          "_by_", sub("p", "", context), "_change_in_",
                          paste0(libName, collapse = "_"), "_vs_", paste0(controlName, collapse = "_"), 
                          "_in_", sortRegion, "_of_", superfamName, "_in_", refbase, "_", paste0(chrName, collapse = "_"), ".bed"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
