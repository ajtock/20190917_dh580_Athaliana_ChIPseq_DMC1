#!/applications/R/R-4.0.0/bin/Rscript

# author: Andy Tock
# contact: ajt200@cam.ac.uk
# date: 01.09.2021

# Calculate and plot metaprofiles of ChIP-seq
# (feature windowed means and 95% confidence intervals, CIs)
# for different feature sets

# Usage:
# /applications/R/R-4.0.0/bin/Rscript genes_crossovers_hypoDMRs4sets_Gypsy_CENAthila_nonCENAthila_CEN180_DNAmeth_6metaprofiles.R 'Chr1,Chr2,Chr3,Chr4,Chr5' 2000 '2kb' 10 10bp '0.02,0.96' 'WT_BSseq_Rep1_2014,WT_BSseq_Rep2_2013,WT_BSseq_Rep3_2013,cmt3_BSseq_Rep1,kss_BSseq_Rep1,met1_BSseq_Rep1' 'BSseq_leaf_Stroud_Jacobsen_2013_Cell_2014_NSMB/snakemake_BSseq_t2t-col.20210610,BSseq_leaf_Stroud_Jacobsen_2013_Cell_2014_NSMB/snakemake_BSseq_t2t-col.20210610,BSseq_leaf_Stroud_Jacobsen_2013_Cell_2014_NSMB/snakemake_BSseq_t2t-col.20210610,BSseq_leaf_Stroud_Jacobsen_2013_Cell_2014_NSMB/snakemake_BSseq_t2t-col.20210610,BSseq_leaf_Stroud_Jacobsen_2013_Cell_2014_NSMB/snakemake_BSseq_t2t-col.20210610,BSseq_leaf_Stroud_Jacobsen_2013_Cell_2014_NSMB/snakemake_BSseq_t2t-col.20210610' 'WT BS-seq Rep1,WT BS-seq Rep2,WT BS-seq Rep3,cmt3 BS-seq,kss BS-seq,met1 BS-seq' 'navy,blue,deepskyblue,green2,darkorange,magenta' CHG Gypsy_LTR 

#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split = ","))
#upstream <- 2000
#downstream <- 2000
#flankName <- "2kb"
#binSize <- 10
#binName <- "10bp"
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
#legendPos <- as.numeric(unlist(strsplit("0.02,0.40",
#                                        split = ",")))
#ChIPNames <- unlist(strsplit("WT_BSseq_Rep1_2014,WT_BSseq_Rep2_2013,WT_BSseq_Rep3_2013,cmt3_BSseq_Rep1,kss_BSseq_Rep1,met1_BSseq_Rep1",
#                             split = ","))
#ChIPNamesDir <- unlist(strsplit("BSseq_leaf_Stroud_Jacobsen_2013_Cell_2014_NSMB/snakemake_BSseq_t2t-col.20210610,BSseq_leaf_Stroud_Jacobsen_2013_Cell_2014_NSMB/snakemake_BSseq_t2t-col.20210610,BSseq_leaf_Stroud_Jacobsen_2013_Cell_2014_NSMB/snakemake_BSseq_t2t-col.20210610,BSseq_leaf_Stroud_Jacobsen_2013_Cell_2014_NSMB/snakemake_BSseq_t2t-col.20210610,BSseq_leaf_Stroud_Jacobsen_2013_Cell_2014_NSMB/snakemake_BSseq_t2t-col.20210610,BSseq_leaf_Stroud_Jacobsen_2013_Cell_2014_NSMB/snakemake_BSseq_t2t-col.20210610",
#                                split = ","))
#ChIPNamesPlot <- unlist(strsplit("WT BS-seq Rep1,WT BS-seq Rep2,WT BS-seq Rep3,cmt3 BS-seq,kss BS-seq,met1 BS-seq",
#                                 split = ","))
#ChIPColours <- unlist(strsplit("navy,blue,deepskyblue,green2,darkorange,magenta",
#                               split = ","))
#context <- "CHG"
#TEsf <- "Gypsy_LTR"

args <- commandArgs(trailingOnly = T)
chrName <- unlist(strsplit(args[1],
                           split = ","))
upstream <- as.numeric(args[2])
downstream <- as.numeric(args[2])
flankName <- args[3]
binSize <- as.numeric(args[4])
binName <- args[5]
legendPos <- as.numeric(unlist(strsplit(args[6],
                                        split = ",")))
ChIPNames <- unlist(strsplit(args[7],
                             split = ","))
ChIPNamesDir <- unlist(strsplit(args[8],
                                split = ","))
ChIPNamesPlot <- unlist(strsplit(args[9],
                                 split = ","))
ChIPColours <- unlist(strsplit(args[10],
                               split = ","))
context <- args[11]
TEsf <- args[12]

options(stringsAsFactors = F)
library(parallel)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(grid)
library(gridExtra)
library(extrafont)
#extrafont::loadfonts()

outDir <- paste0(paste0(chrName, collapse = "_"), "/")
plotDir <- paste0(outDir, "plots/")
system(paste0("[ -d ", outDir, " ] || mkdir -p ", outDir))
system(paste0("[ -d ", plotDir, " ] || mkdir -p ", plotDir))

geneNamePlot <- "Genes"
ranLocNamePlot <- "Random loci"
crossoverNamePlot <- "Crossovers"
cmt3_hypoCHG_DMRNamePlot <- "cmt3 hypoCHG DMRs"
kss_hypoCHG_DMRNamePlot <- "kss hypoCHG DMRs"
cmt3_hypoCHH_DMRNamePlot <- "cmt3 hypoCHH DMRs"
kss_hypoCHH_DMRNamePlot <- "kss hypoCHH DMRs"
TEsfNamePlot <- gsub("_", " ", TEsf)
CENAthilaNamePlot <- "CEN Athila"
nonCENAthilaNamePlot <- "NonCEN Athila"
CEN180NamePlot <- "CEN180"

featureNamePlot <- c(
                     geneNamePlot,
                     ranLocNamePlot,
                     crossoverNamePlot,
                     cmt3_hypoCHG_DMRNamePlot,
                     kss_hypoCHG_DMRNamePlot,
                     cmt3_hypoCHH_DMRNamePlot,
                     kss_hypoCHH_DMRNamePlot,
                     TEsfNamePlot,
                     CENAthilaNamePlot,
                     nonCENAthilaNamePlot,
                     CEN180NamePlot
                    )

gene_bodyLength <- 2200
crossover_bodyLength <- 1000
hypoDMR_bodyLength <- 1000
TEsf_bodyLength <- 2000
CEN180_bodyLength <- 180
#ranLoc_bodyLength <- 180
ranLoc_bodyLength <- 2200

# Define feature start and end labels for plotting
featureStartLab <- "Start"
featureEndLab <- "End"
geneStartLab <- "TSS"
geneEndLab <- "TTS"

# Load feature matrices for each dataset
ChIPDirs <- sapply(seq_along(ChIPNames), function(x) {
  paste0("/home/ajt200/analysis/",
         ChIPNamesDir[x],
         "/coverage/")
})

## ChIP
# gene
ChIP_geneMats <- mclapply(seq_along(ChIPNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(ChIPDirs[x], "geneProfiles/matrices/",
                                ChIPNames[x],
                                "_MappedOn_t2t-col.20210610_", context, "_genes_in_",
                                chrName[y], "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(ChIPNames))
# If genes from multiple chromosomes are to be analysed,
# concatenate the corresponding gene coverage matrices
ChIP_geneMats <- mclapply(seq_along(ChIP_geneMats), function(x) {
  if(length(chrName) > 1) {
    do.call(rbind, ChIP_geneMats[[x]])
  } else {
    ChIP_geneMats[[x]][[1]]
  }
}, mc.cores = length(ChIP_geneMats))

# crossover
ChIP_crossoverMats <- mclapply(seq_along(ChIPNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(ChIPDirs[x], "crossoverProfiles/matrices/",
                                ChIPNames[x],
                                "_MappedOn_t2t-col.20210610_", context, "_crossovers_in_",
                                chrName[y], "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(ChIPNames))
# If crossovers from multiple chromosomes are to be analysed,
# concatenate the corresponding crossover coverage matrices
ChIP_crossoverMats <- mclapply(seq_along(ChIP_crossoverMats), function(x) {
  if(length(chrName) > 1) {
    do.call(rbind, ChIP_crossoverMats[[x]])
  } else {
    ChIP_crossoverMats[[x]][[1]]
  }
}, mc.cores = length(ChIP_crossoverMats))

# cmt3_hypoCHG_DMR
ChIP_cmt3_hypoCHG_DMRMats <- mclapply(seq_along(ChIPNames), function(x) {
  as.matrix(read.table(paste0(ChIPDirs[x], "hypoDMRprofiles/matrices/",
                              ChIPNames[x],
                              "_MappedOn_t2t-col.20210610_", context, "_cmt3_hypoCHG_DMRs_in_",
                              paste0(chrName, collapse = "_"), "_genomewide_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                       header = F, skip = 3))
}, mc.cores = length(ChIPNames))

# kss_hypoCHG_DMR
ChIP_kss_hypoCHG_DMRMats <- mclapply(seq_along(ChIPNames), function(x) {
  as.matrix(read.table(paste0(ChIPDirs[x], "hypoDMRprofiles/matrices/",
                              ChIPNames[x],
                              "_MappedOn_t2t-col.20210610_", context, "_kss_hypoCHG_DMRs_in_",
                              paste0(chrName, collapse = "_"), "_genomewide_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                       header = F, skip = 3))
}, mc.cores = length(ChIPNames))

# cmt3_hypoCHH_DMR
ChIP_cmt3_hypoCHH_DMRMats <- mclapply(seq_along(ChIPNames), function(x) {
  as.matrix(read.table(paste0(ChIPDirs[x], "hypoDMRprofiles/matrices/",
                              ChIPNames[x],
                              "_MappedOn_t2t-col.20210610_", context, "_cmt3_hypoCHH_DMRs_in_",
                              paste0(chrName, collapse = "_"), "_genomewide_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                       header = F, skip = 3))
}, mc.cores = length(ChIPNames))

# kss_hypoCHH_DMR
ChIP_kss_hypoCHH_DMRMats <- mclapply(seq_along(ChIPNames), function(x) {
  as.matrix(read.table(paste0(ChIPDirs[x], "hypoDMRprofiles/matrices/",
                              ChIPNames[x],
                              "_MappedOn_t2t-col.20210610_", context, "_kss_hypoCHH_DMRs_in_",
                              paste0(chrName, collapse = "_"), "_genomewide_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                       header = F, skip = 3))
}, mc.cores = length(ChIPNames))

# TEsf
ChIP_TEsfMats <- mclapply(seq_along(ChIPNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(ChIPDirs[x], "TEprofiles/matrices/",
                                ChIPNames[x],
                                "_MappedOn_t2t-col.20210610_", context, "_TEs_", TEsf, "_in_",
                                chrName[y], "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(ChIPNames))
# If TEsfs from multiple chromosomes are to be analysed,
# concatenate the corresponding TEsf coverage matrices
ChIP_TEsfMats <- mclapply(seq_along(ChIP_TEsfMats), function(x) {
  if(length(chrName) > 1) {
    do.call(rbind, ChIP_TEsfMats[[x]])
  } else {
    ChIP_TEsfMats[[x]][[1]]
  }
}, mc.cores = length(ChIP_TEsfMats))

# CENAthila
ChIP_CENAthilaMats <- mclapply(seq_along(ChIPNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(ChIPDirs[x], "CENAthilaProfiles/matrices/",
                                ChIPNames[x],
                                "_MappedOn_t2t-col.20210610_", context, "_CENAthila_in_",
                                chrName[y], "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(ChIPNames))
# If CENAthilas from multiple chromosomes are to be analysed,
# concatenate the corresponding CENAthila coverage matrices
ChIP_CENAthilaMats <- mclapply(seq_along(ChIP_CENAthilaMats), function(x) {
  if(length(chrName) > 1) {
    do.call(rbind, ChIP_CENAthilaMats[[x]])
  } else {
    ChIP_CENAthilaMats[[x]][[1]]
  }
}, mc.cores = length(ChIP_CENAthilaMats))

# nonCENAthila
ChIP_nonCENAthilaMats <- mclapply(seq_along(ChIPNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(ChIPDirs[x], "CENAthilaProfiles/matrices/",
                                ChIPNames[x],
                                "_MappedOn_t2t-col.20210610_", context, "_nonCENAthila_in_",
                                chrName[y], "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(ChIPNames))
# If nonCENAthilas from multiple chromosomes are to be analysed,
# concatenate the corresponding nonCENAthila coverage matrices
ChIP_nonCENAthilaMats <- mclapply(seq_along(ChIP_nonCENAthilaMats), function(x) {
  if(length(chrName) > 1) {
    do.call(rbind, ChIP_nonCENAthilaMats[[x]])
  } else {
    ChIP_nonCENAthilaMats[[x]][[1]]
  }
}, mc.cores = length(ChIP_nonCENAthilaMats))

# CEN180
ChIP_CEN180Mats <- mclapply(seq_along(ChIPNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(ChIPDirs[x], "CEN180profiles/matrices/",
                                ChIPNames[x],
                                "_MappedOn_t2t-col.20210610_", context, "_CEN180_in_",
                                chrName[y], "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3,
                         colClasses = c(rep("NULL", 1000/binSize), rep(NA, ((1000*2)+(180))/binSize), rep("NULL", 1000/binSize))))
  })
}, mc.cores = length(ChIPNames))
# If CEN180s from multiple chromosomes are to be analysed,
# concatenate the corresponding CEN180 coverage matrices
ChIP_CEN180Mats <- mclapply(seq_along(ChIP_CEN180Mats), function(x) {
  if(length(chrName) > 1) {
    do.call(rbind, ChIP_CEN180Mats[[x]])
  } else {
    ChIP_CEN180Mats[[x]][[1]]
  }
}, mc.cores = length(ChIP_CEN180Mats))

## ranLoc
#ChIP_ranLocMats <- mclapply(seq_along(ChIPNames), function(x) {
#  lapply(seq_along(chrName), function(y) {
#    as.matrix(read.table(paste0(ChIPDirs[x], "CEN180profiles/matrices/",
#                                ChIPNames[x],
#                                "_MappedOn_t2t-col.20210610_", context, "_CEN180_in_",
#                                chrName[y], "_ranLoc_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
#                         header = F, skip = 3,
#                         colClasses = c(rep("NULL", 1000/binSize), rep(NA, ((1000*2)+(180))/binSize), rep("NULL", 1000/binSize))))
#  })
#}, mc.cores = length(ChIPNames))
## If ranLocs from multiple chromosomes are to be analysed,
## concatenate the corresponding ranLoc coverage matrices
#ChIP_ranLocMats <- mclapply(seq_along(ChIP_ranLocMats), function(x) {
#  if(length(chrName) > 1) {
#    do.call(rbind, ChIP_ranLocMats[[x]])
#  } else {
#    ChIP_ranLocMats[[x]][[1]]
#  }
#}, mc.cores = length(ChIP_ranLocMats))

# ranLoc
ChIP_ranLocMats <- mclapply(seq_along(ChIPNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(ChIPDirs[x], "geneProfiles/matrices/",
                                ChIPNames[x],
                                "_MappedOn_t2t-col.20210610_", context, "_genes_in_",
                                chrName[y], "_ranLoc_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(ChIPNames))
# If ranLocs from multiple chromosomes are to be analysed,
# concatenate the corresponding ranLoc coverage matrices
ChIP_ranLocMats <- mclapply(seq_along(ChIP_ranLocMats), function(x) {
  if(length(chrName) > 1) {
    do.call(rbind, ChIP_ranLocMats[[x]])
  } else {
    ChIP_ranLocMats[[x]][[1]]
  }
}, mc.cores = length(ChIP_ranLocMats))



# ChIP
# Add column names
for(x in seq_along(ChIPNames)) {
  print(ChIPNames[x])
  colnames(ChIP_geneMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                    paste0("t", ((upstream/binSize)+1):((upstream+gene_bodyLength)/binSize)),
                                    paste0("d", (((upstream+gene_bodyLength)/binSize)+1):(((upstream+gene_bodyLength)/binSize)+(downstream/binSize))))
  colnames(ChIP_crossoverMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                         paste0("t", ((upstream/binSize)+1):((upstream+crossover_bodyLength)/binSize)),
                                         paste0("d", (((upstream+crossover_bodyLength)/binSize)+1):(((upstream+crossover_bodyLength)/binSize)+(downstream/binSize))))
  colnames(ChIP_cmt3_hypoCHG_DMRMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                                paste0("t", ((upstream/binSize)+1):((upstream+hypoDMR_bodyLength)/binSize)),
                                                paste0("d", (((upstream+hypoDMR_bodyLength)/binSize)+1):(((upstream+hypoDMR_bodyLength)/binSize)+(downstream/binSize))))
  colnames(ChIP_kss_hypoCHG_DMRMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                               paste0("t", ((upstream/binSize)+1):((upstream+hypoDMR_bodyLength)/binSize)),
                                               paste0("d", (((upstream+hypoDMR_bodyLength)/binSize)+1):(((upstream+hypoDMR_bodyLength)/binSize)+(downstream/binSize))))
  colnames(ChIP_cmt3_hypoCHH_DMRMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                                paste0("t", ((upstream/binSize)+1):((upstream+hypoDMR_bodyLength)/binSize)),
                                                paste0("d", (((upstream+hypoDMR_bodyLength)/binSize)+1):(((upstream+hypoDMR_bodyLength)/binSize)+(downstream/binSize))))
  colnames(ChIP_kss_hypoCHH_DMRMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                               paste0("t", ((upstream/binSize)+1):((upstream+hypoDMR_bodyLength)/binSize)),
                                               paste0("d", (((upstream+hypoDMR_bodyLength)/binSize)+1):(((upstream+hypoDMR_bodyLength)/binSize)+(downstream/binSize))))
  colnames(ChIP_TEsfMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                    paste0("t", ((upstream/binSize)+1):((upstream+TEsf_bodyLength)/binSize)),
                                    paste0("d", (((upstream+TEsf_bodyLength)/binSize)+1):(((upstream+TEsf_bodyLength)/binSize)+(downstream/binSize))))
  colnames(ChIP_CENAthilaMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                         paste0("t", ((upstream/binSize)+1):((upstream+TEsf_bodyLength)/binSize)),
                                         paste0("d", (((upstream+TEsf_bodyLength)/binSize)+1):(((upstream+TEsf_bodyLength)/binSize)+(downstream/binSize))))
  colnames(ChIP_nonCENAthilaMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                            paste0("t", ((upstream/binSize)+1):((upstream+TEsf_bodyLength)/binSize)),
                                            paste0("d", (((upstream+TEsf_bodyLength)/binSize)+1):(((upstream+TEsf_bodyLength)/binSize)+(downstream/binSize))))
  colnames(ChIP_CEN180Mats[[x]]) <- c(paste0("u", 1:((upstream-1000)/binSize)),
                                      paste0("t", (((upstream-1000)/binSize)+1):(((upstream-1000)+CEN180_bodyLength)/binSize)),
                                      paste0("d", ((((upstream-1000)+CEN180_bodyLength)/binSize)+1):((((upstream-1000)+CEN180_bodyLength)/binSize)+((downstream-1000)/binSize))))
  colnames(ChIP_ranLocMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                      paste0("t", ((upstream/binSize)+1):((upstream+ranLoc_bodyLength)/binSize)),
                                      paste0("d", (((upstream+ranLoc_bodyLength)/binSize)+1):(((upstream+ranLoc_bodyLength)/binSize)+(downstream/binSize))))
}

# Create list of lists in which each element in the enclosing list corresponds to a library
# and the two elements in the nested list correspond to coverage matrices for features and random loci
ChIP_mats <- mclapply(seq_along(ChIPNames), function(x) {
  list(
       # gene
       ChIP_geneMats[[x]],
       # ranLoc
       ChIP_ranLocMats[[x]],
       # crossover
       ChIP_crossoverMats[[x]],
       # cmt3_hypoCHG_DMRs
       ChIP_cmt3_hypoCHG_DMRMats[[x]],
       # kss_hypoCHG_DMRs
       ChIP_kss_hypoCHG_DMRMats[[x]],
       # cmt3_hypoCHH_DMRs
       ChIP_cmt3_hypoCHH_DMRMats[[x]],
       # kss_hypoCHH_DMRs
       ChIP_kss_hypoCHH_DMRMats[[x]],
       # TEsfs
       ChIP_TEsfMats[[x]],
       # CENAthilas
       ChIP_CENAthilaMats[[x]],
       # nonCENAthilas
       ChIP_nonCENAthilaMats[[x]],
       # CEN180s
       ChIP_CEN180Mats[[x]]
      )
}, mc.cores = length(ChIPNames))

# Transpose matrix and convert into dataframe
# in which first column is window name
wideDFfeature_list_ChIP <- mclapply(seq_along(ChIP_mats), function(x) {
  lapply(seq_along(ChIP_mats[[x]]), function(y) {
    data.frame(window = colnames(ChIP_mats[[x]][[y]]),
               t(ChIP_mats[[x]][[y]]))
  })
}, mc.cores = length(ChIP_mats))

# Convert into tidy data.frame (long format)
tidyDFfeature_list_ChIP  <- mclapply(seq_along(wideDFfeature_list_ChIP), function(x) {
  lapply(seq_along(ChIP_mats[[x]]), function(y) {
    gather(data  = wideDFfeature_list_ChIP[[x]][[y]],
           key   = feature,
           value = coverage,
           -window)
  }) 
}, mc.cores = length(wideDFfeature_list_ChIP))

# Order levels of factor "window" so that sequential levels
# correspond to sequential windows
for(x in seq_along(tidyDFfeature_list_ChIP)) {
  for(y in seq_along(ChIP_mats[[x]])) {
    tidyDFfeature_list_ChIP[[x]][[y]]$window <- factor(tidyDFfeature_list_ChIP[[x]][[y]]$window,
                                                       levels = as.character(wideDFfeature_list_ChIP[[x]][[y]]$window))
  }
}

# Create summary data.frame in which each row corresponds to a window (Column 1),
# Column2 is the number of coverage values (features) per window,
# Column3 is the mean of coverage values per window,
# Column4 is the standard deviation of coverage values per window,
# Column5 is the standard error of the mean of coverage values per window,
# Column6 is the lower bound of the 95% confidence interval, and
# Column7 is the upper bound of the 95% confidence interval
#               n      = tapply(X     = tidyDFfeature_list_ChIP[[x]][[y]]$coverage[which(!is.na(tidyDFfeature_list_ChIP[[x]][[y]]$coverage))],
#                               INDEX = tidyDFfeature_list_ChIP[[x]][[y]]$window[which(!is.na(tidyDFfeature_list_ChIP[[x]][[y]]$coverage))],
summaryDFfeature_list_ChIP  <- mclapply(seq_along(tidyDFfeature_list_ChIP), function(x) {
  lapply(seq_along(ChIP_mats[[x]]), function(y) {
    data.frame(window = as.character(wideDFfeature_list_ChIP[[x]][[y]]$window),
               n      = tapply(X     = tidyDFfeature_list_ChIP[[x]][[y]]$coverage,
                               INDEX = tidyDFfeature_list_ChIP[[x]][[y]]$window,
                               FUN   = length),
               mean   = tapply(X     = tidyDFfeature_list_ChIP[[x]][[y]]$coverage,
                               INDEX = tidyDFfeature_list_ChIP[[x]][[y]]$window,
                               FUN   = mean,
                               na.rm = TRUE),
               sd     = tapply(X     = tidyDFfeature_list_ChIP[[x]][[y]]$coverage,
                               INDEX = tidyDFfeature_list_ChIP[[x]][[y]]$window,
                               FUN   = sd,
                               na.rm = TRUE))
  })
}, mc.cores = length(tidyDFfeature_list_ChIP))

for(x in seq_along(summaryDFfeature_list_ChIP)) {
  for(y in seq_along(ChIP_mats[[x]])) {
    summaryDFfeature_list_ChIP[[x]][[y]]$window <- factor(summaryDFfeature_list_ChIP[[x]][[y]]$window,
                                                          levels = as.character(wideDFfeature_list_ChIP[[x]][[y]]$window))
    summaryDFfeature_list_ChIP[[x]][[y]]$winNo <- factor(1:dim(summaryDFfeature_list_ChIP[[x]][[y]])[1])
    summaryDFfeature_list_ChIP[[x]][[y]]$sem <- summaryDFfeature_list_ChIP[[x]][[y]]$sd/sqrt(summaryDFfeature_list_ChIP[[x]][[y]]$n-1)
    summaryDFfeature_list_ChIP[[x]][[y]]$CI_lower <- summaryDFfeature_list_ChIP[[x]][[y]]$mean -
      qt(0.975, df = summaryDFfeature_list_ChIP[[x]][[y]]$n-1)*summaryDFfeature_list_ChIP[[x]][[y]]$sem
    summaryDFfeature_list_ChIP[[x]][[y]]$CI_upper <- summaryDFfeature_list_ChIP[[x]][[y]]$mean +
      qt(0.975, df = summaryDFfeature_list_ChIP[[x]][[y]]$n-1)*summaryDFfeature_list_ChIP[[x]][[y]]$sem
  }
}

# Convert list of lists summaryDFfeature_list_ChIP into
# a list of single data.frames containing all meta-profiles for plotting
geneTmp <- lapply(seq_along(summaryDFfeature_list_ChIP), function(x) {
  summaryDFfeature_list_ChIP[[x]][[1]]
})
ranLocTmp <- lapply(seq_along(summaryDFfeature_list_ChIP), function(x) {
  summaryDFfeature_list_ChIP[[x]][[2]]
})
crossoverTmp <- lapply(seq_along(summaryDFfeature_list_ChIP), function(x) {
  summaryDFfeature_list_ChIP[[x]][[3]]
})
cmt3_hypoCHG_DMRTmp <- lapply(seq_along(summaryDFfeature_list_ChIP), function(x) {
  summaryDFfeature_list_ChIP[[x]][[4]]
})
kss_hypoCHG_DMRTmp <- lapply(seq_along(summaryDFfeature_list_ChIP), function(x) {
  summaryDFfeature_list_ChIP[[x]][[5]]
})
cmt3_hypoCHH_DMRTmp <- lapply(seq_along(summaryDFfeature_list_ChIP), function(x) {
  summaryDFfeature_list_ChIP[[x]][[6]]
})
kss_hypoCHH_DMRTmp <- lapply(seq_along(summaryDFfeature_list_ChIP), function(x) {
  summaryDFfeature_list_ChIP[[x]][[7]]
})
TEsfTmp <- lapply(seq_along(summaryDFfeature_list_ChIP), function(x) {
  summaryDFfeature_list_ChIP[[x]][[8]]
})
CENAthilaTmp <- lapply(seq_along(summaryDFfeature_list_ChIP), function(x) {
  summaryDFfeature_list_ChIP[[x]][[9]]
})
nonCENAthilaTmp <- lapply(seq_along(summaryDFfeature_list_ChIP), function(x) {
  summaryDFfeature_list_ChIP[[x]][[10]]
})
CEN180Tmp <- lapply(seq_along(summaryDFfeature_list_ChIP), function(x) {
  summaryDFfeature_list_ChIP[[x]][[11]]
})

names(geneTmp) <- ChIPNamesPlot
names(ranLocTmp) <- ChIPNamesPlot
names(crossoverTmp) <- ChIPNamesPlot
names(cmt3_hypoCHG_DMRTmp) <- ChIPNamesPlot
names(kss_hypoCHG_DMRTmp) <- ChIPNamesPlot
names(cmt3_hypoCHH_DMRTmp) <- ChIPNamesPlot
names(kss_hypoCHH_DMRTmp) <- ChIPNamesPlot
names(TEsfTmp) <- ChIPNamesPlot
names(CENAthilaTmp) <- ChIPNamesPlot
names(nonCENAthilaTmp) <- ChIPNamesPlot
names(CEN180Tmp) <- ChIPNamesPlot
summaryDFfeature_ChIP <- list(
  bind_rows(geneTmp, .id = "libName"),
  bind_rows(ranLocTmp, .id = "libName"),
  bind_rows(crossoverTmp, .id = "libName"),
  bind_rows(cmt3_hypoCHG_DMRTmp, .id = "libName"),
  bind_rows(kss_hypoCHG_DMRTmp, .id = "libName"),
  bind_rows(cmt3_hypoCHH_DMRTmp, .id = "libName"),
  bind_rows(kss_hypoCHH_DMRTmp, .id = "libName"),
  bind_rows(TEsfTmp, .id = "libName"),
  bind_rows(CENAthilaTmp, .id = "libName"),
  bind_rows(nonCENAthilaTmp, .id = "libName"),
  bind_rows(CEN180Tmp, .id = "libName")
)
for(x in seq_along(summaryDFfeature_ChIP)) {
  summaryDFfeature_ChIP[[x]]$libName <- factor(summaryDFfeature_ChIP[[x]]$libName,
                                                   levels = ChIPNamesPlot)
}

# Define y-axis limits
ymin_ChIP <- min(c(summaryDFfeature_ChIP[[1]]$CI_lower,
                       summaryDFfeature_ChIP[[2]]$CI_lower,
                       summaryDFfeature_ChIP[[3]]$CI_lower,
                       summaryDFfeature_ChIP[[4]]$CI_lower,
                       summaryDFfeature_ChIP[[5]]$CI_lower,
                       summaryDFfeature_ChIP[[6]]$CI_lower,
                       summaryDFfeature_ChIP[[7]]$CI_lower,
                       summaryDFfeature_ChIP[[8]]$CI_lower,
                       summaryDFfeature_ChIP[[9]]$CI_lower,
                       summaryDFfeature_ChIP[[10]]$CI_lower,
                       summaryDFfeature_ChIP[[11]]$CI_lower),
                 na.rm = T)
ymax_ChIP <- max(c(summaryDFfeature_ChIP[[1]]$CI_upper,
                       summaryDFfeature_ChIP[[2]]$CI_upper,
                       summaryDFfeature_ChIP[[3]]$CI_upper,
                       summaryDFfeature_ChIP[[4]]$CI_upper,
                       summaryDFfeature_ChIP[[5]]$CI_upper,
                       summaryDFfeature_ChIP[[6]]$CI_upper,
                       summaryDFfeature_ChIP[[7]]$CI_upper,
                       summaryDFfeature_ChIP[[8]]$CI_upper,
                       summaryDFfeature_ChIP[[9]]$CI_upper,
                       summaryDFfeature_ChIP[[10]]$CI_upper,
                       summaryDFfeature_ChIP[[11]]$CI_upper),
                 na.rm = T)

# Define legend labels
legendLabs <- lapply(seq_along(ChIPNamesPlot), function(x) {
  grobTree(textGrob(bquote(.(ChIPNamesPlot[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = ChIPColours[x], fontsize = 18)))
})


# Plot average profiles with 95% CI ribbon
## gene
summaryDFfeature <- summaryDFfeature_ChIP[[1]]
ggObj1_combined_ChIP <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = ChIPColours) +
scale_y_continuous(limits = c(ymin_ChIP, ymax_ChIP),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_ChIP[[1]])[1]/length(ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_ChIP[[1]])[1]/length(ChIPNames)),
                 labels = c(paste0("-", flankName),
                            geneStartLab,
                            geneEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_ChIP[[1]])[1]/length(ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("m" * .(context))) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 3.5, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(featureNamePlot[[1]]) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## ranLoc
summaryDFfeature <- summaryDFfeature_ChIP[[2]]
ggObj2_combined_ChIP <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = ChIPColours) +
scale_y_continuous(limits = c(ymin_ChIP, ymax_ChIP),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_ChIP[[2]])[1]/length(ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_ChIP[[2]])[1]/length(ChIPNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_ChIP[[2]])[1]/length(ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("m" * .(context))) +
annotation_custom(legendLabs[[1]]) +
annotation_custom(legendLabs[[2]]) +
annotation_custom(legendLabs[[3]]) +
annotation_custom(legendLabs[[4]]) +
annotation_custom(legendLabs[[5]]) +
annotation_custom(legendLabs[[6]]) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 3.5, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(featureNamePlot[[2]]) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## crossover
summaryDFfeature <- summaryDFfeature_ChIP[[3]]
ggObj3_combined_ChIP <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = ChIPColours) +
scale_y_continuous(limits = c(ymin_ChIP, ymax_ChIP),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_ChIP[[3]])[1]/length(ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_ChIP[[3]])[1]/length(ChIPNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_ChIP[[3]])[1]/length(ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("m" * .(context))) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 3.5, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(featureNamePlot[[3]]) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## cmt3_hypoCHG_DMR
summaryDFfeature <- summaryDFfeature_ChIP[[4]]
ggObj4_combined_ChIP <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = ChIPColours) +
scale_y_continuous(limits = c(ymin_ChIP, ymax_ChIP),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_ChIP[[4]])[1]/length(ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_ChIP[[4]])[1]/length(ChIPNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_ChIP[[4]])[1]/length(ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("m" * .(context))) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 3.5, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(featureNamePlot[[4]]) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## kss_hypoCHG_DMR
summaryDFfeature <- summaryDFfeature_ChIP[[5]]
ggObj5_combined_ChIP <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = ChIPColours) +
scale_y_continuous(limits = c(ymin_ChIP, ymax_ChIP),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_ChIP[[5]])[1]/length(ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_ChIP[[5]])[1]/length(ChIPNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_ChIP[[5]])[1]/length(ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("m" * .(context))) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 3.5, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(featureNamePlot[[5]]) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## cmt3_hypoCHH_DMR
summaryDFfeature <- summaryDFfeature_ChIP[[6]]
ggObj6_combined_ChIP <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = ChIPColours) +
scale_y_continuous(limits = c(ymin_ChIP, ymax_ChIP),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_ChIP[[6]])[1]/length(ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_ChIP[[6]])[1]/length(ChIPNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_ChIP[[6]])[1]/length(ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("m" * .(context))) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 3.5, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(featureNamePlot[[6]]) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## kss_hypoCHH_DMR
summaryDFfeature <- summaryDFfeature_ChIP[[7]]
ggObj7_combined_ChIP <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = ChIPColours) +
scale_y_continuous(limits = c(ymin_ChIP, ymax_ChIP),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_ChIP[[7]])[1]/length(ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_ChIP[[7]])[1]/length(ChIPNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_ChIP[[7]])[1]/length(ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("m" * .(context))) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 3.5, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(featureNamePlot[[7]]) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## TEsf
summaryDFfeature <- summaryDFfeature_ChIP[[8]]
ggObj8_combined_ChIP <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = ChIPColours) +
scale_y_continuous(limits = c(ymin_ChIP, ymax_ChIP),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_ChIP[[8]])[1]/length(ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_ChIP[[8]])[1]/length(ChIPNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_ChIP[[8]])[1]/length(ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("m" * .(context))) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 3.5, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(featureNamePlot[[8]]) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## CENAthila
summaryDFfeature <- summaryDFfeature_ChIP[[9]]
ggObj9_combined_ChIP <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = ChIPColours) +
scale_y_continuous(limits = c(ymin_ChIP, ymax_ChIP),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_ChIP[[9]])[1]/length(ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_ChIP[[9]])[1]/length(ChIPNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_ChIP[[9]])[1]/length(ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("m" * .(context))) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 3.5, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(featureNamePlot[[9]]) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## nonCENAthila
summaryDFfeature <- summaryDFfeature_ChIP[[10]]
ggObj10_combined_ChIP <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = ChIPColours) +
scale_y_continuous(limits = c(ymin_ChIP, ymax_ChIP),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_ChIP[[10]])[1]/length(ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_ChIP[[10]])[1]/length(ChIPNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_ChIP[[10]])[1]/length(ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("m" * .(context))) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 3.5, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(featureNamePlot[[10]]) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## CEN180
summaryDFfeature <- summaryDFfeature_ChIP[[11]]
ggObj11_combined_ChIP <- ggplot(data = summaryDFfeature,
                                    mapping = aes(x = winNo,
                                                  y = mean,
                                                  group = libName)
                                   ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = ChIPColours) +
scale_y_continuous(limits = c(ymin_ChIP, ymax_ChIP),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            ((upstream-1000)/binSize)+1,
                            (dim(summaryDFfeature_ChIP[[11]])[1]/length(ChIPNames))-((downstream-1000)/binSize),
                            dim(summaryDFfeature_ChIP[[11]])[1]/length(ChIPNames)),
                 labels = c(paste0("-", "1kb"),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", "1kb"))) +
geom_vline(xintercept = c(((upstream-1000)/binSize)+1,
                          (dim(summaryDFfeature_ChIP[[11]])[1]/length(ChIPNames))-((downstream-1000)/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("m" * .(context))) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 3.5, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(featureNamePlot[[11]]) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

ggObjGA_combined <- grid.arrange(grobs = list(
                                              ggObj1_combined_ChIP,
                                              ggObj2_combined_ChIP,
                                              ggObj3_combined_ChIP,
                                              ggObj4_combined_ChIP,
                                              ggObj5_combined_ChIP,
                                              ggObj6_combined_ChIP,
                                              ggObj7_combined_ChIP,
                                              ggObj8_combined_ChIP,
                                              ggObj9_combined_ChIP,
                                              ggObj10_combined_ChIP,
                                              ggObj11_combined_ChIP
                                             ),
                                 layout_matrix = cbind(
                                                       1,
                                                       2,
                                                       3,
                                                       4,
                                                       5,
                                                       6,
                                                       7,
                                                       8,
                                                       9,
                                                       10,
                                                       11
                                                      ))
ggsave(paste0(plotDir,
              "DNAmeth_", context, "_",
              paste0(ChIPNames[c(1, 4)], collapse = "_"),
              "_avgProfiles_around",
              "_genes_ranLoc_COs_hypoDMRs4sets_", TEsf, "_CENAthila_nonCENAthila_CEN180_in_t2t-col.20210610_",
              paste0(chrName, collapse = "_"), ".pdf"),
       plot = ggObjGA_combined,
       height = 6.5, width = 7*11, limitsize = FALSE)
