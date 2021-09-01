#!/applications/R/R-4.0.0/bin/Rscript

# author: Andy Tock
# contact: ajt200@cam.ac.uk
# date: 01.09.2021

# Calculate and plot metaprofiles of ChIP-seq
# (feature windowed means and 95% confidence intervals, CIs)
# for different feature sets

# Usage:
# /applications/R/R-4.0.0/bin/Rscript genes_crossovers_hypoDMRs4sets_Gypsy_CENAthila_nonCENAthila_CEN180_3metaprofiles.R 'Chr1,Chr2,Chr3,Chr4,Chr5' both 2000 '2kb' 10 10bp '0.02,0.96' 'Col_DMC1_V5_Rep2_ChIP,cmt3_DMC1_V5_Rep1_ChIP,kss_DMC1_V5_Rep1_ChIP' '20190917_dh580_Athaliana_ChIPseq_DMC1/fastq_pooled/snakemake_ChIPseq_t2t-col.20210610,20190917_dh580_Athaliana_ChIPseq_DMC1/fastq_pooled/snakemake_ChIPseq_t2t-col.20210610,20190917_dh580_Athaliana_ChIPseq_DMC1/fastq_pooled/snakemake_ChIPseq_t2t-col.20210610' 'WT DMC1,cmt3 DMC1,kss DMC1' 'navy,green2,darkorange' 'ChIP' Gypsy_LTR

#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split = ","))
#align <- "both"
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
#ChIPNames <- unlist(strsplit("Col_DMC1_V5_Rep2_ChIP,cmt3_DMC1_V5_Rep1_ChIP,kss_DMC1_V5_Rep1_ChIP",
#                             split = ","))
#ChIPNamesDir <- unlist(strsplit("20190917_dh580_Athaliana_ChIPseq_DMC1/fastq_pooled/snakemake_ChIPseq_t2t-col.20210610,20190917_dh580_Athaliana_ChIPseq_DMC1/fastq_pooled/snakemake_ChIPseq_t2t-col.20210610,20190917_dh580_Athaliana_ChIPseq_DMC1/fastq_pooled/snakemake_ChIPseq_t2t-col.20210610",
#                                split = ","))
#ChIPNamesPlot <- unlist(strsplit("WT DMC1,cmt3 DMC1,kss DMC1",
#                                 split = ","))
#ChIPColours <- unlist(strsplit("navy,green2,darkorange",
#                               split = ","))
#yLabPlot <- "DMC1 ChIP"
#TEsf <- "Gypsy_LTR"

args <- commandArgs(trailingOnly = T)
chrName <- unlist(strsplit(args[1],
                           split = ","))
align <- args[2]
upstream <- as.numeric(args[3])
downstream <- as.numeric(args[3])
flankName <- args[4]
binSize <- as.numeric(args[5])
binName <- args[6]
legendPos <- as.numeric(unlist(strsplit(args[7],
                                        split = ",")))
ChIPNames <- unlist(strsplit(args[8],
                             split = ","))
ChIPNamesDir <- unlist(strsplit(args[9],
                                split = ","))
ChIPNamesPlot <- unlist(strsplit(args[10],
                                 split = ","))
ChIPColours <- unlist(strsplit(args[11],
                               split = ","))
yLabPlot <- args[12]
TEsf <- args[13]

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

# Load feature matrices for each chromatin dataset, calculate log2(ChIP/control),
log2ChIPNames <- ChIPNames
log2ChIPNamesPlot <- ChIPNamesPlot
log2ChIPColours <- ChIPColours
# Load feature matrices for each chromatin dataset, calculate log2(ChIP/control),
ChIPDirs <- sapply(seq_along(ChIPNames), function(x) {
  paste0("/home/ajt200/analysis/",
         ChIPNamesDir[x],
         "/mapped/")
})

controlNames <- c(
                  "Col_REC8_Myc_Rep1_input",
                  "WT_gDNA_Rep1_R1"
                 )
controlNamesDir <- c(
                     "20190917_dh580_Athaliana_ChIPseq_DMC1/fastq_pooled/snakemake_ChIPseq_t2t-col.20210610",
                     "150701_Natasha_gDNA/WT/R1/snakemake_SPO11oligos_t2t-col.20210610"
                    )
controlNamesPlot <- c(
                      "PE input (REC8)",
                      "SE gDNA"
                     )
controlColours <- c(
                    "grey10",
                    "grey90"
                   )
controlDirs <- sapply(seq_along(controlNames), function(x) {
  paste0("/home/ajt200/analysis/",
         controlNamesDir[x],
         "/mapped/")
})


## ChIP
# gene
ChIP_geneMats <- mclapply(seq_along(ChIPNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(ChIPDirs[x], "geneProfiles/matrices/",
                                ChIPNames[x],
                                "_MappedOn_t2t-col.20210610_lowXM_", align, "_sort_norm_genes_in_",
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
                                "_MappedOn_t2t-col.20210610_lowXM_", align, "_sort_norm_crossovers_in_",
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
                              "_MappedOn_t2t-col.20210610_lowXM_", align, "_sort_norm_cmt3_hypoCHG_DMRs_in_",
                              paste0(chrName, collapse = "_"), "_genomewide_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                       header = F, skip = 3))
}, mc.cores = length(ChIPNames))

# kss_hypoCHG_DMR
ChIP_kss_hypoCHG_DMRMats <- mclapply(seq_along(ChIPNames), function(x) {
  as.matrix(read.table(paste0(ChIPDirs[x], "hypoDMRprofiles/matrices/",
                              ChIPNames[x],
                              "_MappedOn_t2t-col.20210610_lowXM_", align, "_sort_norm_kss_hypoCHG_DMRs_in_",
                              paste0(chrName, collapse = "_"), "_genomewide_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                       header = F, skip = 3))
}, mc.cores = length(ChIPNames))

# cmt3_hypoCHH_DMR
ChIP_cmt3_hypoCHH_DMRMats <- mclapply(seq_along(ChIPNames), function(x) {
  as.matrix(read.table(paste0(ChIPDirs[x], "hypoDMRprofiles/matrices/",
                              ChIPNames[x],
                              "_MappedOn_t2t-col.20210610_lowXM_", align, "_sort_norm_cmt3_hypoCHH_DMRs_in_",
                              paste0(chrName, collapse = "_"), "_genomewide_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                       header = F, skip = 3))
}, mc.cores = length(ChIPNames))

# kss_hypoCHH_DMR
ChIP_kss_hypoCHH_DMRMats <- mclapply(seq_along(ChIPNames), function(x) {
  as.matrix(read.table(paste0(ChIPDirs[x], "hypoDMRprofiles/matrices/",
                              ChIPNames[x],
                              "_MappedOn_t2t-col.20210610_lowXM_", align, "_sort_norm_kss_hypoCHH_DMRs_in_",
                              paste0(chrName, collapse = "_"), "_genomewide_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                       header = F, skip = 3))
}, mc.cores = length(ChIPNames))

# TEsf
ChIP_TEsfMats <- mclapply(seq_along(ChIPNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(ChIPDirs[x], "TEprofiles/matrices/",
                                ChIPNames[x],
                                "_MappedOn_t2t-col.20210610_lowXM_", align, "_sort_norm_TEs_", TEsf, "_in_",
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
                                "_MappedOn_t2t-col.20210610_lowXM_", align, "_sort_norm_CENAthila_in_",
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
                                "_MappedOn_t2t-col.20210610_lowXM_", align, "_sort_norm_nonCENAthila_in_",
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
                                "_MappedOn_t2t-col.20210610_lowXM_", align, "_sort_norm_CEN180_in_",
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
#                                "_MappedOn_t2t-col.20210610_lowXM_", align, "_sort_norm_CEN180_in_",
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
                                "_MappedOn_t2t-col.20210610_lowXM_", align, "_sort_norm_genes_in_",
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


## control
# gene
control_geneMats <- mclapply(seq_along(controlNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(controlDirs[x], "geneProfiles/matrices/",
                                controlNames[x],
                                "_MappedOn_t2t-col.20210610_lowXM_", align, "_sort_norm_genes_in_",
                                chrName[y], "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(controlNames))
# If genes from multiple chromosomes are to be analysed,
# concatenate the corresponding gene coverage matrices
control_geneMats <- mclapply(seq_along(control_geneMats), function(x) {
  if(length(chrName) > 1) {
    do.call(rbind, control_geneMats[[x]])
  } else {
    control_geneMats[[x]][[1]]
  }
}, mc.cores = length(control_geneMats))

# crossover
control_crossoverMats <- mclapply(seq_along(controlNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(controlDirs[x], "crossoverProfiles/matrices/",
                                controlNames[x],
                                "_MappedOn_t2t-col.20210610_lowXM_", align, "_sort_norm_crossovers_in_",
                                chrName[y], "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(controlNames))
# If crossovers from multiple chromosomes are to be analysed,
# concatenate the corresponding crossover coverage matrices
control_crossoverMats <- mclapply(seq_along(control_crossoverMats), function(x) {
  if(length(chrName) > 1) {
    do.call(rbind, control_crossoverMats[[x]])
  } else {
    control_crossoverMats[[x]][[1]]
  }
}, mc.cores = length(control_crossoverMats))

# cmt3_hypoCHG_DMR
control_cmt3_hypoCHG_DMRMats <- mclapply(seq_along(controlNames), function(x) {
  as.matrix(read.table(paste0(controlDirs[x], "hypoDMRprofiles/matrices/",
                              controlNames[x],
                              "_MappedOn_t2t-col.20210610_lowXM_", align, "_sort_norm_cmt3_hypoCHG_DMRs_in_",
                              paste0(chrName, collapse = "_"), "_genomewide_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                       header = F, skip = 3))
}, mc.cores = length(controlNames))

# kss_hypoCHG_DMR
control_kss_hypoCHG_DMRMats <- mclapply(seq_along(controlNames), function(x) {
  as.matrix(read.table(paste0(controlDirs[x], "hypoDMRprofiles/matrices/",
                              controlNames[x],
                              "_MappedOn_t2t-col.20210610_lowXM_", align, "_sort_norm_kss_hypoCHG_DMRs_in_",
                              paste0(chrName, collapse = "_"), "_genomewide_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                       header = F, skip = 3))
}, mc.cores = length(controlNames))

# cmt3_hypoCHH_DMR
control_cmt3_hypoCHH_DMRMats <- mclapply(seq_along(controlNames), function(x) {
  as.matrix(read.table(paste0(controlDirs[x], "hypoDMRprofiles/matrices/",
                              controlNames[x],
                              "_MappedOn_t2t-col.20210610_lowXM_", align, "_sort_norm_cmt3_hypoCHH_DMRs_in_",
                              paste0(chrName, collapse = "_"), "_genomewide_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                       header = F, skip = 3))
}, mc.cores = length(controlNames))

# kss_hypoCHH_DMR
control_kss_hypoCHH_DMRMats <- mclapply(seq_along(controlNames), function(x) {
  as.matrix(read.table(paste0(controlDirs[x], "hypoDMRprofiles/matrices/",
                              controlNames[x],
                              "_MappedOn_t2t-col.20210610_lowXM_", align, "_sort_norm_kss_hypoCHH_DMRs_in_",
                              paste0(chrName, collapse = "_"), "_genomewide_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                       header = F, skip = 3))
}, mc.cores = length(controlNames))

# TEsf
control_TEsfMats <- mclapply(seq_along(controlNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(controlDirs[x], "TEprofiles/matrices/",
                                controlNames[x],
                                "_MappedOn_t2t-col.20210610_lowXM_", align, "_sort_norm_TEs_", TEsf, "_in_",
                                chrName[y], "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(controlNames))
# If TEsfs from multiple chromosomes are to be analysed,
# concatenate the corresponding TEsf coverage matrices
control_TEsfMats <- mclapply(seq_along(control_TEsfMats), function(x) {
  if(length(chrName) > 1) {
    do.call(rbind, control_TEsfMats[[x]])
  } else {
    control_TEsfMats[[x]][[1]]
  }
}, mc.cores = length(control_TEsfMats))

# CENAthila
control_CENAthilaMats <- mclapply(seq_along(controlNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(controlDirs[x], "CENAthilaProfiles/matrices/",
                                controlNames[x],
                                "_MappedOn_t2t-col.20210610_lowXM_", align, "_sort_norm_CENAthila_in_",
                                chrName[y], "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(controlNames))
# If CENAthilas from multiple chromosomes are to be analysed,
# concatenate the corresponding CENAthila coverage matrices
control_CENAthilaMats <- mclapply(seq_along(control_CENAthilaMats), function(x) {
  if(length(chrName) > 1) {
    do.call(rbind, control_CENAthilaMats[[x]])
  } else {
    control_CENAthilaMats[[x]][[1]]
  }
}, mc.cores = length(control_CENAthilaMats))

# nonCENAthila
control_nonCENAthilaMats <- mclapply(seq_along(controlNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(controlDirs[x], "CENAthilaProfiles/matrices/",
                                controlNames[x],
                                "_MappedOn_t2t-col.20210610_lowXM_", align, "_sort_norm_nonCENAthila_in_",
                                chrName[y], "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(controlNames))
# If nonCENAthilas from multiple chromosomes are to be analysed,
# concatenate the corresponding nonCENAthila coverage matrices
control_nonCENAthilaMats <- mclapply(seq_along(control_nonCENAthilaMats), function(x) {
  if(length(chrName) > 1) {
    do.call(rbind, control_nonCENAthilaMats[[x]])
  } else {
    control_nonCENAthilaMats[[x]][[1]]
  }
}, mc.cores = length(control_nonCENAthilaMats))

# CEN180
control_CEN180Mats <- mclapply(seq_along(controlNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(controlDirs[x], "CEN180profiles/matrices/",
                                controlNames[x],
                                "_MappedOn_t2t-col.20210610_lowXM_", align, "_sort_norm_CEN180_in_",
                                chrName[y], "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3,
                         colClasses = c(rep("NULL", 1000/binSize), rep(NA, ((1000*2)+(180))/binSize), rep("NULL", 1000/binSize))))
  })
}, mc.cores = length(controlNames))
# If CEN180s from multiple chromosomes are to be analysed,
# concatenate the corresponding CEN180 coverage matrices
control_CEN180Mats <- mclapply(seq_along(control_CEN180Mats), function(x) {
  if(length(chrName) > 1) {
    do.call(rbind, control_CEN180Mats[[x]])
  } else {
    control_CEN180Mats[[x]][[1]]
  }
}, mc.cores = length(control_CEN180Mats))

## ranLoc
#control_ranLocMats <- mclapply(seq_along(controlNames), function(x) {
#  lapply(seq_along(chrName), function(y) {
#    as.matrix(read.table(paste0(controlDirs[x], "CEN180profiles/matrices/",
#                                controlNames[x],
#                                "_MappedOn_t2t-col.20210610_lowXM_", align, "_sort_norm_CEN180_in_",
#                                chrName[y], "_ranLoc_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
#                         header = F, skip = 3,
#                         colClasses = c(rep("NULL", 1000/binSize), rep(NA, ((1000*2)+(180))/binSize), rep("NULL", 1000/binSize))))
#  })
#}, mc.cores = length(controlNames))
## If ranLocs from multiple chromosomes are to be analysed,
## concatenate the corresponding ranLoc coverage matrices
#control_ranLocMats <- mclapply(seq_along(control_ranLocMats), function(x) {
#  if(length(chrName) > 1) {
#    do.call(rbind, control_ranLocMats[[x]])
#  } else {
#    control_ranLocMats[[x]][[1]]
#  }
#}, mc.cores = length(control_ranLocMats))

# ranLoc
control_ranLocMats <- mclapply(seq_along(controlNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(controlDirs[x], "geneProfiles/matrices/",
                                controlNames[x],
                                "_MappedOn_t2t-col.20210610_lowXM_", align, "_sort_norm_genes_in_",
                                chrName[y], "_ranLoc_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(controlNames))
# If ranLocs from multiple chromosomes are to be analysed,
# concatenate the corresponding ranLoc coverage matrices
control_ranLocMats <- mclapply(seq_along(control_ranLocMats), function(x) {
  if(length(chrName) > 1) {
    do.call(rbind, control_ranLocMats[[x]])
  } else {
    control_ranLocMats[[x]][[1]]
  }
}, mc.cores = length(control_ranLocMats))


# Conditionally calculate log2(ChIP/control)
# for each matrix depending on library
# gene
log2ChIP_geneMats <- mclapply(seq_along(ChIP_geneMats), function(x) {
  if ( ChIPNames[x] %in% c("Col_DMC1_V5_Rep1_ChIP", "Col_DMC1_V5_Rep2_ChIP", "cmt3_DMC1_V5_Rep1_ChIP", "kss_DMC1_V5_Rep1_ChIP") ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[which(grepl("Col_REC8_Myc_Rep1_input", controlNames))], " for log2((ChIP+1)/(control+1))"))
    log2((ChIP_geneMats[[x]]+1)/(control_geneMats[[which(grepl("Col_REC8_Myc_Rep1_input", controlNames))]]+1))
  } else if ( grepl("SPO11oligos", ChIPNames[x]) ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[which(grepl("WT_gDNA_Rep1_R1", controlNames))], " for log2((ChIP+1)/(control+1))"))
    log2((ChIP_geneMats[[x]]+1)/(control_geneMats[[which(grepl("WT_gDNA_Rep1_R1", controlNames))]]+1))
  }
}, mc.cores = length(ChIP_geneMats))

# crossover
log2ChIP_crossoverMats <- mclapply(seq_along(ChIP_crossoverMats), function(x) {
  if ( ChIPNames[x] %in% c("Col_DMC1_V5_Rep1_ChIP", "Col_DMC1_V5_Rep2_ChIP", "cmt3_DMC1_V5_Rep1_ChIP", "kss_DMC1_V5_Rep1_ChIP") ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[which(grepl("Col_REC8_Myc_Rep1_input", controlNames))], " for log2((ChIP+1)/(control+1))"))
    log2((ChIP_crossoverMats[[x]]+1)/(control_crossoverMats[[which(grepl("Col_REC8_Myc_Rep1_input", controlNames))]]+1))
  } else if ( grepl("SPO11oligos", ChIPNames[x]) ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[which(grepl("WT_gDNA_Rep1_R1", controlNames))], " for log2((ChIP+1)/(control+1))"))
    log2((ChIP_crossoverMats[[x]]+1)/(control_crossoverMats[[which(grepl("WT_gDNA_Rep1_R1", controlNames))]]+1))
  }
}, mc.cores = length(ChIP_crossoverMats))

# cmt3_hypoCHG_DMR
log2ChIP_cmt3_hypoCHG_DMRMats <- mclapply(seq_along(ChIP_cmt3_hypoCHG_DMRMats), function(x) {
  if ( ChIPNames[x] %in% c("Col_DMC1_V5_Rep1_ChIP", "Col_DMC1_V5_Rep2_ChIP", "cmt3_DMC1_V5_Rep1_ChIP", "kss_DMC1_V5_Rep1_ChIP") ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[which(grepl("Col_REC8_Myc_Rep1_input", controlNames))], " for log2((ChIP+1)/(control+1))"))
    log2((ChIP_cmt3_hypoCHG_DMRMats[[x]]+1)/(control_cmt3_hypoCHG_DMRMats[[which(grepl("Col_REC8_Myc_Rep1_input", controlNames))]]+1))
  } else if ( grepl("SPO11oligos", ChIPNames[x]) ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[which(grepl("WT_gDNA_Rep1_R1", controlNames))], " for log2((ChIP+1)/(control+1))"))
    log2((ChIP_cmt3_hypoCHG_DMRMats[[x]]+1)/(control_cmt3_hypoCHG_DMRMats[[which(grepl("WT_gDNA_Rep1_R1", controlNames))]]+1))
  }
}, mc.cores = length(ChIP_cmt3_hypoCHG_DMRMats))

# kss_hypoCHG_DMR
log2ChIP_kss_hypoCHG_DMRMats <- mclapply(seq_along(ChIP_kss_hypoCHG_DMRMats), function(x) {
  if ( ChIPNames[x] %in% c("Col_DMC1_V5_Rep1_ChIP", "Col_DMC1_V5_Rep2_ChIP", "cmt3_DMC1_V5_Rep1_ChIP", "kss_DMC1_V5_Rep1_ChIP") ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[which(grepl("Col_REC8_Myc_Rep1_input", controlNames))], " for log2((ChIP+1)/(control+1))"))
    log2((ChIP_kss_hypoCHG_DMRMats[[x]]+1)/(control_kss_hypoCHG_DMRMats[[which(grepl("Col_REC8_Myc_Rep1_input", controlNames))]]+1))
  } else if ( grepl("SPO11oligos", ChIPNames[x]) ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[which(grepl("WT_gDNA_Rep1_R1", controlNames))], " for log2((ChIP+1)/(control+1))"))
    log2((ChIP_kss_hypoCHG_DMRMats[[x]]+1)/(control_kss_hypoCHG_DMRMats[[which(grepl("WT_gDNA_Rep1_R1", controlNames))]]+1))
  }
}, mc.cores = length(ChIP_kss_hypoCHG_DMRMats))

# cmt3_hypoCHH_DMR
log2ChIP_cmt3_hypoCHH_DMRMats <- mclapply(seq_along(ChIP_cmt3_hypoCHH_DMRMats), function(x) {
  if ( ChIPNames[x] %in% c("Col_DMC1_V5_Rep1_ChIP", "Col_DMC1_V5_Rep2_ChIP", "cmt3_DMC1_V5_Rep1_ChIP", "kss_DMC1_V5_Rep1_ChIP") ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[which(grepl("Col_REC8_Myc_Rep1_input", controlNames))], " for log2((ChIP+1)/(control+1))"))
    log2((ChIP_cmt3_hypoCHH_DMRMats[[x]]+1)/(control_cmt3_hypoCHH_DMRMats[[which(grepl("Col_REC8_Myc_Rep1_input", controlNames))]]+1))
  } else if ( grepl("SPO11oligos", ChIPNames[x]) ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[which(grepl("WT_gDNA_Rep1_R1", controlNames))], " for log2((ChIP+1)/(control+1))"))
    log2((ChIP_cmt3_hypoCHH_DMRMats[[x]]+1)/(control_cmt3_hypoCHH_DMRMats[[which(grepl("WT_gDNA_Rep1_R1", controlNames))]]+1))
  }
}, mc.cores = length(ChIP_cmt3_hypoCHH_DMRMats))

# kss_hypoCHH_DMR
log2ChIP_kss_hypoCHH_DMRMats <- mclapply(seq_along(ChIP_kss_hypoCHH_DMRMats), function(x) {
  if ( ChIPNames[x] %in% c("Col_DMC1_V5_Rep1_ChIP", "Col_DMC1_V5_Rep2_ChIP", "cmt3_DMC1_V5_Rep1_ChIP", "kss_DMC1_V5_Rep1_ChIP") ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[which(grepl("Col_REC8_Myc_Rep1_input", controlNames))], " for log2((ChIP+1)/(control+1))"))
    log2((ChIP_kss_hypoCHH_DMRMats[[x]]+1)/(control_kss_hypoCHH_DMRMats[[which(grepl("Col_REC8_Myc_Rep1_input", controlNames))]]+1))
  } else if ( grepl("SPO11oligos", ChIPNames[x]) ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[which(grepl("WT_gDNA_Rep1_R1", controlNames))], " for log2((ChIP+1)/(control+1))"))
    log2((ChIP_kss_hypoCHH_DMRMats[[x]]+1)/(control_kss_hypoCHH_DMRMats[[which(grepl("WT_gDNA_Rep1_R1", controlNames))]]+1))
  }
}, mc.cores = length(ChIP_kss_hypoCHH_DMRMats))

# TEsf
log2ChIP_TEsfMats <- mclapply(seq_along(ChIP_TEsfMats), function(x) {
  if ( ChIPNames[x] %in% c("Col_DMC1_V5_Rep1_ChIP", "Col_DMC1_V5_Rep2_ChIP", "cmt3_DMC1_V5_Rep1_ChIP", "kss_DMC1_V5_Rep1_ChIP") ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[which(grepl("Col_REC8_Myc_Rep1_input", controlNames))], " for log2((ChIP+1)/(control+1))"))
    log2((ChIP_TEsfMats[[x]]+1)/(control_TEsfMats[[which(grepl("Col_REC8_Myc_Rep1_input", controlNames))]]+1))
  } else if ( grepl("SPO11oligos", ChIPNames[x]) ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[which(grepl("WT_gDNA_Rep1_R1", controlNames))], " for log2((ChIP+1)/(control+1))"))
    log2((ChIP_TEsfMats[[x]]+1)/(control_TEsfMats[[which(grepl("WT_gDNA_Rep1_R1", controlNames))]]+1))
  }
}, mc.cores = length(ChIP_TEsfMats))

# CENAthila
log2ChIP_CENAthilaMats <- mclapply(seq_along(ChIP_CENAthilaMats), function(x) {
  if ( ChIPNames[x] %in% c("Col_DMC1_V5_Rep1_ChIP", "Col_DMC1_V5_Rep2_ChIP", "cmt3_DMC1_V5_Rep1_ChIP", "kss_DMC1_V5_Rep1_ChIP") ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[which(grepl("Col_REC8_Myc_Rep1_input", controlNames))], " for log2((ChIP+1)/(control+1))"))
    log2((ChIP_CENAthilaMats[[x]]+1)/(control_CENAthilaMats[[which(grepl("Col_REC8_Myc_Rep1_input", controlNames))]]+1))
  } else if ( grepl("SPO11oligos", ChIPNames[x]) ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[which(grepl("WT_gDNA_Rep1_R1", controlNames))], " for log2((ChIP+1)/(control+1))"))
    log2((ChIP_CENAthilaMats[[x]]+1)/(control_CENAthilaMats[[which(grepl("WT_gDNA_Rep1_R1", controlNames))]]+1))
  }
}, mc.cores = length(ChIP_CENAthilaMats))

# nonCENAthila
log2ChIP_nonCENAthilaMats <- mclapply(seq_along(ChIP_nonCENAthilaMats), function(x) {
  if ( ChIPNames[x] %in% c("Col_DMC1_V5_Rep1_ChIP", "Col_DMC1_V5_Rep2_ChIP", "cmt3_DMC1_V5_Rep1_ChIP", "kss_DMC1_V5_Rep1_ChIP") ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[which(grepl("Col_REC8_Myc_Rep1_input", controlNames))], " for log2((ChIP+1)/(control+1))"))
    log2((ChIP_nonCENAthilaMats[[x]]+1)/(control_nonCENAthilaMats[[which(grepl("Col_REC8_Myc_Rep1_input", controlNames))]]+1))
  } else if ( grepl("SPO11oligos", ChIPNames[x]) ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[which(grepl("WT_gDNA_Rep1_R1", controlNames))], " for log2((ChIP+1)/(control+1))"))
    log2((ChIP_nonCENAthilaMats[[x]]+1)/(control_nonCENAthilaMats[[which(grepl("WT_gDNA_Rep1_R1", controlNames))]]+1))
  }
}, mc.cores = length(ChIP_nonCENAthilaMats))

# CEN180
log2ChIP_CEN180Mats <- mclapply(seq_along(ChIP_CEN180Mats), function(x) {
  if ( ChIPNames[x] %in% c("Col_DMC1_V5_Rep1_ChIP", "Col_DMC1_V5_Rep2_ChIP", "cmt3_DMC1_V5_Rep1_ChIP", "kss_DMC1_V5_Rep1_ChIP") ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[which(grepl("Col_REC8_Myc_Rep1_input", controlNames))], " for log2((ChIP+1)/(control+1))"))
    log2((ChIP_CEN180Mats[[x]]+1)/(control_CEN180Mats[[which(grepl("Col_REC8_Myc_Rep1_input", controlNames))]]+1))
  } else if ( grepl("SPO11oligos", ChIPNames[x]) ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[which(grepl("WT_gDNA_Rep1_R1", controlNames))], " for log2((ChIP+1)/(control+1))"))
    log2((ChIP_CEN180Mats[[x]]+1)/(control_CEN180Mats[[which(grepl("WT_gDNA_Rep1_R1", controlNames))]]+1))
  }
}, mc.cores = length(ChIP_CEN180Mats))

# ranLoc
log2ChIP_ranLocMats <- mclapply(seq_along(ChIP_ranLocMats), function(x) {
  if ( ChIPNames[x] %in% c("Col_DMC1_V5_Rep1_ChIP", "Col_DMC1_V5_Rep2_ChIP", "cmt3_DMC1_V5_Rep1_ChIP", "kss_DMC1_V5_Rep1_ChIP") ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[which(grepl("Col_REC8_Myc_Rep1_input", controlNames))], " for log2((ChIP+1)/(control+1))"))
    log2((ChIP_ranLocMats[[x]]+1)/(control_ranLocMats[[which(grepl("Col_REC8_Myc_Rep1_input", controlNames))]]+1))
  } else if ( grepl("SPO11oligos", ChIPNames[x]) ) {
    print(paste0(ChIPNames[x], " library; using ", controlNames[which(grepl("WT_gDNA_Rep1_R1", controlNames))], " for log2((ChIP+1)/(control+1))"))
    log2((ChIP_ranLocMats[[x]]+1)/(control_ranLocMats[[which(grepl("WT_gDNA_Rep1_R1", controlNames))]]+1))
  }
}, mc.cores = length(ChIP_ranLocMats))


# log2ChIP
# Add column names
for(x in seq_along(log2ChIPNames)) {
  print(log2ChIPNames[x])
  colnames(log2ChIP_geneMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                        paste0("t", ((upstream/binSize)+1):((upstream+gene_bodyLength)/binSize)),
                                        paste0("d", (((upstream+gene_bodyLength)/binSize)+1):(((upstream+gene_bodyLength)/binSize)+(downstream/binSize))))
  colnames(log2ChIP_crossoverMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                             paste0("t", ((upstream/binSize)+1):((upstream+crossover_bodyLength)/binSize)),
                                             paste0("d", (((upstream+crossover_bodyLength)/binSize)+1):(((upstream+crossover_bodyLength)/binSize)+(downstream/binSize))))
  colnames(log2ChIP_cmt3_hypoCHG_DMRMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                                    paste0("t", ((upstream/binSize)+1):((upstream+hypoDMR_bodyLength)/binSize)),
                                                    paste0("d", (((upstream+hypoDMR_bodyLength)/binSize)+1):(((upstream+hypoDMR_bodyLength)/binSize)+(downstream/binSize))))
  colnames(log2ChIP_kss_hypoCHG_DMRMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                                   paste0("t", ((upstream/binSize)+1):((upstream+hypoDMR_bodyLength)/binSize)),
                                                   paste0("d", (((upstream+hypoDMR_bodyLength)/binSize)+1):(((upstream+hypoDMR_bodyLength)/binSize)+(downstream/binSize))))
  colnames(log2ChIP_cmt3_hypoCHH_DMRMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                                    paste0("t", ((upstream/binSize)+1):((upstream+hypoDMR_bodyLength)/binSize)),
                                                    paste0("d", (((upstream+hypoDMR_bodyLength)/binSize)+1):(((upstream+hypoDMR_bodyLength)/binSize)+(downstream/binSize))))
  colnames(log2ChIP_kss_hypoCHH_DMRMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                                   paste0("t", ((upstream/binSize)+1):((upstream+hypoDMR_bodyLength)/binSize)),
                                                   paste0("d", (((upstream+hypoDMR_bodyLength)/binSize)+1):(((upstream+hypoDMR_bodyLength)/binSize)+(downstream/binSize))))
  colnames(log2ChIP_TEsfMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                        paste0("t", ((upstream/binSize)+1):((upstream+TEsf_bodyLength)/binSize)),
                                        paste0("d", (((upstream+TEsf_bodyLength)/binSize)+1):(((upstream+TEsf_bodyLength)/binSize)+(downstream/binSize))))
  colnames(log2ChIP_CENAthilaMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                             paste0("t", ((upstream/binSize)+1):((upstream+TEsf_bodyLength)/binSize)),
                                             paste0("d", (((upstream+TEsf_bodyLength)/binSize)+1):(((upstream+TEsf_bodyLength)/binSize)+(downstream/binSize))))
  colnames(log2ChIP_nonCENAthilaMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                                paste0("t", ((upstream/binSize)+1):((upstream+TEsf_bodyLength)/binSize)),
                                                paste0("d", (((upstream+TEsf_bodyLength)/binSize)+1):(((upstream+TEsf_bodyLength)/binSize)+(downstream/binSize))))
  colnames(log2ChIP_CEN180Mats[[x]]) <- c(paste0("u", 1:((upstream-1000)/binSize)),
                                          paste0("t", (((upstream-1000)/binSize)+1):(((upstream-1000)+CEN180_bodyLength)/binSize)),
                                          paste0("d", ((((upstream-1000)+CEN180_bodyLength)/binSize)+1):((((upstream-1000)+CEN180_bodyLength)/binSize)+((downstream-1000)/binSize))))
  colnames(log2ChIP_ranLocMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                          paste0("t", ((upstream/binSize)+1):((upstream+ranLoc_bodyLength)/binSize)),
                                          paste0("d", (((upstream+ranLoc_bodyLength)/binSize)+1):(((upstream+ranLoc_bodyLength)/binSize)+(downstream/binSize))))
}

# Create list of lists in which each element in the enclosing list corresponds to a library
# and the two elements in the nested list correspond to coverage matrices for features and random loci
log2ChIP_mats <- mclapply(seq_along(log2ChIPNames), function(x) {
  list(
       # gene
       log2ChIP_geneMats[[x]],
       # ranLoc
       log2ChIP_ranLocMats[[x]],
       # crossover
       log2ChIP_crossoverMats[[x]],
       # cmt3_hypoCHG_DMRs
       log2ChIP_cmt3_hypoCHG_DMRMats[[x]],
       # kss_hypoCHG_DMRs
       log2ChIP_kss_hypoCHG_DMRMats[[x]],
       # cmt3_hypoCHH_DMRs
       log2ChIP_cmt3_hypoCHH_DMRMats[[x]],
       # kss_hypoCHH_DMRs
       log2ChIP_kss_hypoCHH_DMRMats[[x]],
       # TEsfs
       log2ChIP_TEsfMats[[x]],
       # CENAthilas
       log2ChIP_CENAthilaMats[[x]],
       # nonCENAthilas
       log2ChIP_nonCENAthilaMats[[x]],
       # CEN180s
       log2ChIP_CEN180Mats[[x]]
      )
}, mc.cores = length(log2ChIPNames))

# Transpose matrix and convert into dataframe
# in which first column is window name
wideDFfeature_list_log2ChIP <- mclapply(seq_along(log2ChIP_mats), function(x) {
  lapply(seq_along(log2ChIP_mats[[x]]), function(y) {
    data.frame(window = colnames(log2ChIP_mats[[x]][[y]]),
               t(log2ChIP_mats[[x]][[y]]))
  })
}, mc.cores = length(log2ChIP_mats))

# Convert into tidy data.frame (long format)
tidyDFfeature_list_log2ChIP  <- mclapply(seq_along(wideDFfeature_list_log2ChIP), function(x) {
  lapply(seq_along(log2ChIP_mats[[x]]), function(y) {
    gather(data  = wideDFfeature_list_log2ChIP[[x]][[y]],
           key   = feature,
           value = coverage,
           -window)
  }) 
}, mc.cores = length(wideDFfeature_list_log2ChIP))

# Order levels of factor "window" so that sequential levels
# correspond to sequential windows
for(x in seq_along(tidyDFfeature_list_log2ChIP)) {
  for(y in seq_along(log2ChIP_mats[[x]])) {
    tidyDFfeature_list_log2ChIP[[x]][[y]]$window <- factor(tidyDFfeature_list_log2ChIP[[x]][[y]]$window,
                                                           levels = as.character(wideDFfeature_list_log2ChIP[[x]][[y]]$window))
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
  lapply(seq_along(log2ChIP_mats[[x]]), function(y) {
    data.frame(window = as.character(wideDFfeature_list_log2ChIP[[x]][[y]]$window),
               n      = tapply(X     = tidyDFfeature_list_log2ChIP[[x]][[y]]$coverage,
                               INDEX = tidyDFfeature_list_log2ChIP[[x]][[y]]$window,
                               FUN   = length),
               mean   = tapply(X     = tidyDFfeature_list_log2ChIP[[x]][[y]]$coverage,
                               INDEX = tidyDFfeature_list_log2ChIP[[x]][[y]]$window,
                               FUN   = mean,
                               na.rm = TRUE),
               sd     = tapply(X     = tidyDFfeature_list_log2ChIP[[x]][[y]]$coverage,
                               INDEX = tidyDFfeature_list_log2ChIP[[x]][[y]]$window,
                               FUN   = sd,
                               na.rm = TRUE))
  })
}, mc.cores = length(tidyDFfeature_list_log2ChIP))

for(x in seq_along(summaryDFfeature_list_log2ChIP)) {
  for(y in seq_along(log2ChIP_mats[[x]])) {
    summaryDFfeature_list_log2ChIP[[x]][[y]]$window <- factor(summaryDFfeature_list_log2ChIP[[x]][[y]]$window,
                                                              levels = as.character(wideDFfeature_list_log2ChIP[[x]][[y]]$window))
    summaryDFfeature_list_log2ChIP[[x]][[y]]$winNo <- factor(1:dim(summaryDFfeature_list_log2ChIP[[x]][[y]])[1])
    summaryDFfeature_list_log2ChIP[[x]][[y]]$sem <- summaryDFfeature_list_log2ChIP[[x]][[y]]$sd/sqrt(summaryDFfeature_list_log2ChIP[[x]][[y]]$n-1)
    summaryDFfeature_list_log2ChIP[[x]][[y]]$CI_lower <- summaryDFfeature_list_log2ChIP[[x]][[y]]$mean -
      qt(0.975, df = summaryDFfeature_list_log2ChIP[[x]][[y]]$n-1)*summaryDFfeature_list_log2ChIP[[x]][[y]]$sem
    summaryDFfeature_list_log2ChIP[[x]][[y]]$CI_upper <- summaryDFfeature_list_log2ChIP[[x]][[y]]$mean +
      qt(0.975, df = summaryDFfeature_list_log2ChIP[[x]][[y]]$n-1)*summaryDFfeature_list_log2ChIP[[x]][[y]]$sem
  }
}

# Convert list of lists summaryDFfeature_list_log2ChIP into
# a list of single data.frames containing all meta-profiles for plotting
geneTmp <- lapply(seq_along(summaryDFfeature_list_log2ChIP), function(x) {
  summaryDFfeature_list_log2ChIP[[x]][[1]]
})
ranLocTmp <- lapply(seq_along(summaryDFfeature_list_log2ChIP), function(x) {
  summaryDFfeature_list_log2ChIP[[x]][[2]]
})
crossoverTmp <- lapply(seq_along(summaryDFfeature_list_log2ChIP), function(x) {
  summaryDFfeature_list_log2ChIP[[x]][[3]]
})
cmt3_hypoCHG_DMRTmp <- lapply(seq_along(summaryDFfeature_list_log2ChIP), function(x) {
  summaryDFfeature_list_log2ChIP[[x]][[4]]
})
kss_hypoCHG_DMRTmp <- lapply(seq_along(summaryDFfeature_list_log2ChIP), function(x) {
  summaryDFfeature_list_log2ChIP[[x]][[5]]
})
cmt3_hypoCHH_DMRTmp <- lapply(seq_along(summaryDFfeature_list_log2ChIP), function(x) {
  summaryDFfeature_list_log2ChIP[[x]][[6]]
})
kss_hypoCHH_DMRTmp <- lapply(seq_along(summaryDFfeature_list_log2ChIP), function(x) {
  summaryDFfeature_list_log2ChIP[[x]][[7]]
})
TEsfTmp <- lapply(seq_along(summaryDFfeature_list_log2ChIP), function(x) {
  summaryDFfeature_list_log2ChIP[[x]][[8]]
})
CENAthilaTmp <- lapply(seq_along(summaryDFfeature_list_log2ChIP), function(x) {
  summaryDFfeature_list_log2ChIP[[x]][[9]]
})
nonCENAthilaTmp <- lapply(seq_along(summaryDFfeature_list_log2ChIP), function(x) {
  summaryDFfeature_list_log2ChIP[[x]][[10]]
})
CEN180Tmp <- lapply(seq_along(summaryDFfeature_list_log2ChIP), function(x) {
  summaryDFfeature_list_log2ChIP[[x]][[11]]
})

names(geneTmp) <- log2ChIPNamesPlot
names(ranLocTmp) <- log2ChIPNamesPlot
names(crossoverTmp) <- log2ChIPNamesPlot
names(cmt3_hypoCHG_DMRTmp) <- log2ChIPNamesPlot
names(kss_hypoCHG_DMRTmp) <- log2ChIPNamesPlot
names(cmt3_hypoCHH_DMRTmp) <- log2ChIPNamesPlot
names(kss_hypoCHH_DMRTmp) <- log2ChIPNamesPlot
names(TEsfTmp) <- log2ChIPNamesPlot
names(CENAthilaTmp) <- log2ChIPNamesPlot
names(nonCENAthilaTmp) <- log2ChIPNamesPlot
names(CEN180Tmp) <- log2ChIPNamesPlot
summaryDFfeature_log2ChIP <- list(
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
for(x in seq_along(summaryDFfeature_log2ChIP)) {
  summaryDFfeature_log2ChIP[[x]]$libName <- factor(summaryDFfeature_log2ChIP[[x]]$libName,
                                                   levels = log2ChIPNamesPlot)
}

# Define y-axis limits
ymin_log2ChIP <- min(c(summaryDFfeature_log2ChIP[[1]]$CI_lower,
                       summaryDFfeature_log2ChIP[[2]]$CI_lower,
                       summaryDFfeature_log2ChIP[[3]]$CI_lower,
                       summaryDFfeature_log2ChIP[[4]]$CI_lower,
                       summaryDFfeature_log2ChIP[[5]]$CI_lower,
                       summaryDFfeature_log2ChIP[[6]]$CI_lower,
                       summaryDFfeature_log2ChIP[[7]]$CI_lower,
                       summaryDFfeature_log2ChIP[[8]]$CI_lower,
                       summaryDFfeature_log2ChIP[[9]]$CI_lower,
                       summaryDFfeature_log2ChIP[[10]]$CI_lower,
                       summaryDFfeature_log2ChIP[[11]]$CI_lower),
                 na.rm = T)
ymax_log2ChIP <- max(c(summaryDFfeature_log2ChIP[[1]]$CI_upper,
                       summaryDFfeature_log2ChIP[[2]]$CI_upper,
                       summaryDFfeature_log2ChIP[[3]]$CI_upper,
                       summaryDFfeature_log2ChIP[[4]]$CI_upper,
                       summaryDFfeature_log2ChIP[[5]]$CI_upper,
                       summaryDFfeature_log2ChIP[[6]]$CI_upper,
                       summaryDFfeature_log2ChIP[[7]]$CI_upper,
                       summaryDFfeature_log2ChIP[[8]]$CI_upper,
                       summaryDFfeature_log2ChIP[[9]]$CI_upper,
                       summaryDFfeature_log2ChIP[[10]]$CI_upper,
                       summaryDFfeature_log2ChIP[[11]]$CI_upper),
                 na.rm = T)

# Define legend labels
legendLabs <- lapply(seq_along(log2ChIPNamesPlot), function(x) {
  grobTree(textGrob(bquote(.(log2ChIPNamesPlot[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = log2ChIPColours[x], fontsize = 18)))
})


# Plot average profiles with 95% CI ribbon
## gene
summaryDFfeature <- summaryDFfeature_log2ChIP[[1]]
ggObj1_combined_log2ChIP <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = log2ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = log2ChIPColours) +
scale_y_continuous(limits = c(ymin_log2ChIP, ymax_log2ChIP),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[1]])[1]/length(log2ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_log2ChIP[[1]])[1]/length(log2ChIPNames)),
                 labels = c(paste0("-", flankName),
                            geneStartLab,
                            geneEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_log2ChIP[[1]])[1]/length(log2ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Log"[2] * "(" * .(yLabPlot) * "/control)")) +
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
summaryDFfeature <- summaryDFfeature_log2ChIP[[2]]
ggObj2_combined_log2ChIP <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = log2ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = log2ChIPColours) +
scale_y_continuous(limits = c(ymin_log2ChIP, ymax_log2ChIP),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[2]])[1]/length(log2ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_log2ChIP[[2]])[1]/length(log2ChIPNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_log2ChIP[[2]])[1]/length(log2ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Log"[2] * "(" * .(yLabPlot) * "/control)")) +
annotation_custom(legendLabs[[1]]) +
annotation_custom(legendLabs[[2]]) +
annotation_custom(legendLabs[[3]]) +
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
summaryDFfeature <- summaryDFfeature_log2ChIP[[3]]
ggObj3_combined_log2ChIP <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = log2ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = log2ChIPColours) +
scale_y_continuous(limits = c(ymin_log2ChIP, ymax_log2ChIP),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[3]])[1]/length(log2ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_log2ChIP[[3]])[1]/length(log2ChIPNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_log2ChIP[[3]])[1]/length(log2ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Log"[2] * "(" * .(yLabPlot) * "/control)")) +
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
summaryDFfeature <- summaryDFfeature_log2ChIP[[4]]
ggObj4_combined_log2ChIP <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = log2ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = log2ChIPColours) +
scale_y_continuous(limits = c(ymin_log2ChIP, ymax_log2ChIP),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[4]])[1]/length(log2ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_log2ChIP[[4]])[1]/length(log2ChIPNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_log2ChIP[[4]])[1]/length(log2ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Log"[2] * "(" * .(yLabPlot) * "/control)")) +
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
summaryDFfeature <- summaryDFfeature_log2ChIP[[5]]
ggObj5_combined_log2ChIP <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = log2ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = log2ChIPColours) +
scale_y_continuous(limits = c(ymin_log2ChIP, ymax_log2ChIP),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[5]])[1]/length(log2ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_log2ChIP[[5]])[1]/length(log2ChIPNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_log2ChIP[[5]])[1]/length(log2ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Log"[2] * "(" * .(yLabPlot) * "/control)")) +
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
summaryDFfeature <- summaryDFfeature_log2ChIP[[6]]
ggObj6_combined_log2ChIP <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = log2ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = log2ChIPColours) +
scale_y_continuous(limits = c(ymin_log2ChIP, ymax_log2ChIP),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[6]])[1]/length(log2ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_log2ChIP[[6]])[1]/length(log2ChIPNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_log2ChIP[[6]])[1]/length(log2ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Log"[2] * "(" * .(yLabPlot) * "/control)")) +
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
summaryDFfeature <- summaryDFfeature_log2ChIP[[7]]
ggObj7_combined_log2ChIP <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = log2ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = log2ChIPColours) +
scale_y_continuous(limits = c(ymin_log2ChIP, ymax_log2ChIP),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[7]])[1]/length(log2ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_log2ChIP[[7]])[1]/length(log2ChIPNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_log2ChIP[[7]])[1]/length(log2ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Log"[2] * "(" * .(yLabPlot) * "/control)")) +
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
summaryDFfeature <- summaryDFfeature_log2ChIP[[8]]
ggObj8_combined_log2ChIP <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = log2ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = log2ChIPColours) +
scale_y_continuous(limits = c(ymin_log2ChIP, ymax_log2ChIP),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[8]])[1]/length(log2ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_log2ChIP[[8]])[1]/length(log2ChIPNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_log2ChIP[[8]])[1]/length(log2ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Log"[2] * "(" * .(yLabPlot) * "/control)")) +
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
summaryDFfeature <- summaryDFfeature_log2ChIP[[9]]
ggObj9_combined_log2ChIP <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = log2ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = log2ChIPColours) +
scale_y_continuous(limits = c(ymin_log2ChIP, ymax_log2ChIP),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[9]])[1]/length(log2ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_log2ChIP[[9]])[1]/length(log2ChIPNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_log2ChIP[[9]])[1]/length(log2ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Log"[2] * "(" * .(yLabPlot) * "/control)")) +
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
summaryDFfeature <- summaryDFfeature_log2ChIP[[10]]
ggObj10_combined_log2ChIP <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = log2ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = log2ChIPColours) +
scale_y_continuous(limits = c(ymin_log2ChIP, ymax_log2ChIP),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[10]])[1]/length(log2ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_log2ChIP[[10]])[1]/length(log2ChIPNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_log2ChIP[[10]])[1]/length(log2ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Log"[2] * "(" * .(yLabPlot) * "/control)")) +
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
summaryDFfeature <- summaryDFfeature_log2ChIP[[11]]
ggObj11_combined_log2ChIP <- ggplot(data = summaryDFfeature,
                                    mapping = aes(x = winNo,
                                                  y = mean,
                                                  group = libName)
                                   ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = log2ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = log2ChIPColours) +
scale_y_continuous(limits = c(ymin_log2ChIP, ymax_log2ChIP),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            ((upstream-1000)/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[11]])[1]/length(log2ChIPNames))-((downstream-1000)/binSize),
                            dim(summaryDFfeature_log2ChIP[[11]])[1]/length(log2ChIPNames)),
                 labels = c(paste0("-", "1kb"),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", "1kb"))) +
geom_vline(xintercept = c(((upstream-1000)/binSize)+1,
                          (dim(summaryDFfeature_log2ChIP[[11]])[1]/length(log2ChIPNames))-((downstream-1000)/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Log"[2] * "(" * .(yLabPlot) * "/control)")) +
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
                                              ggObj1_combined_log2ChIP,
                                              ggObj2_combined_log2ChIP,
                                              ggObj3_combined_log2ChIP,
                                              ggObj4_combined_log2ChIP,
                                              ggObj5_combined_log2ChIP,
                                              ggObj6_combined_log2ChIP,
                                              ggObj7_combined_log2ChIP,
                                              ggObj8_combined_log2ChIP,
                                              ggObj9_combined_log2ChIP,
                                              ggObj10_combined_log2ChIP,
                                              ggObj11_combined_log2ChIP
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
              "log2ChIPcontrol_",
              paste0(log2ChIPNames, collapse = "_"),
              "_avgProfiles_around",
              "_genes_ranLoc_COs_hypoDMRs4sets_", TEsf, "_CENAthila_nonCENAthila_CEN180_in_t2t-col.20210610_",
              paste0(chrName, collapse = "_"), "_", align, ".pdf"),
       plot = ggObjGA_combined,
       height = 6.5, width = 7*11, limitsize = FALSE)


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
     y = bquote(.(yLabPlot))) +
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
     y = bquote(.(yLabPlot))) +
annotation_custom(legendLabs[[1]]) +
annotation_custom(legendLabs[[2]]) +
annotation_custom(legendLabs[[3]]) +
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
     y = bquote(.(yLabPlot))) +
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
     y = bquote(.(yLabPlot))) +
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
     y = bquote(.(yLabPlot))) +
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
     y = bquote(.(yLabPlot))) +
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
     y = bquote(.(yLabPlot))) +
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
     y = bquote(.(yLabPlot))) +
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
     y = bquote(.(yLabPlot))) +
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
     y = bquote(.(yLabPlot))) +
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
     y = bquote(.(yLabPlot))) +
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
              "ChIP_",
              paste0(ChIPNames, collapse = "_"),
              "_avgProfiles_around",
              "_genes_ranLoc_COs_hypoDMRs4sets_", TEsf, "_CENAthila_nonCENAthila_CEN180_in_t2t-col.20210610_",
              paste0(chrName, collapse = "_"), "_", align, ".pdf"),
       plot = ggObjGA_combined,
       height = 6.5, width = 7*11, limitsize = FALSE)


# control
# Add column names
for(x in seq_along(controlNames)) {
  print(controlNames[x])
  colnames(control_geneMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                        paste0("t", ((upstream/binSize)+1):((upstream+gene_bodyLength)/binSize)),
                                        paste0("d", (((upstream+gene_bodyLength)/binSize)+1):(((upstream+gene_bodyLength)/binSize)+(downstream/binSize))))
  colnames(control_crossoverMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                             paste0("t", ((upstream/binSize)+1):((upstream+crossover_bodyLength)/binSize)),
                                             paste0("d", (((upstream+crossover_bodyLength)/binSize)+1):(((upstream+crossover_bodyLength)/binSize)+(downstream/binSize))))
  colnames(control_cmt3_hypoCHG_DMRMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                                    paste0("t", ((upstream/binSize)+1):((upstream+hypoDMR_bodyLength)/binSize)),
                                                    paste0("d", (((upstream+hypoDMR_bodyLength)/binSize)+1):(((upstream+hypoDMR_bodyLength)/binSize)+(downstream/binSize))))
  colnames(control_kss_hypoCHG_DMRMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                                   paste0("t", ((upstream/binSize)+1):((upstream+hypoDMR_bodyLength)/binSize)),
                                                   paste0("d", (((upstream+hypoDMR_bodyLength)/binSize)+1):(((upstream+hypoDMR_bodyLength)/binSize)+(downstream/binSize))))
  colnames(control_cmt3_hypoCHH_DMRMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                                    paste0("t", ((upstream/binSize)+1):((upstream+hypoDMR_bodyLength)/binSize)),
                                                    paste0("d", (((upstream+hypoDMR_bodyLength)/binSize)+1):(((upstream+hypoDMR_bodyLength)/binSize)+(downstream/binSize))))
  colnames(control_kss_hypoCHH_DMRMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                                   paste0("t", ((upstream/binSize)+1):((upstream+hypoDMR_bodyLength)/binSize)),
                                                   paste0("d", (((upstream+hypoDMR_bodyLength)/binSize)+1):(((upstream+hypoDMR_bodyLength)/binSize)+(downstream/binSize))))
  colnames(control_TEsfMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                        paste0("t", ((upstream/binSize)+1):((upstream+TEsf_bodyLength)/binSize)),
                                        paste0("d", (((upstream+TEsf_bodyLength)/binSize)+1):(((upstream+TEsf_bodyLength)/binSize)+(downstream/binSize))))
  colnames(control_CENAthilaMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                             paste0("t", ((upstream/binSize)+1):((upstream+TEsf_bodyLength)/binSize)),
                                             paste0("d", (((upstream+TEsf_bodyLength)/binSize)+1):(((upstream+TEsf_bodyLength)/binSize)+(downstream/binSize))))
  colnames(control_nonCENAthilaMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                                paste0("t", ((upstream/binSize)+1):((upstream+TEsf_bodyLength)/binSize)),
                                                paste0("d", (((upstream+TEsf_bodyLength)/binSize)+1):(((upstream+TEsf_bodyLength)/binSize)+(downstream/binSize))))
  colnames(control_CEN180Mats[[x]]) <- c(paste0("u", 1:((upstream-1000)/binSize)),
                                          paste0("t", (((upstream-1000)/binSize)+1):(((upstream-1000)+CEN180_bodyLength)/binSize)),
                                          paste0("d", ((((upstream-1000)+CEN180_bodyLength)/binSize)+1):((((upstream-1000)+CEN180_bodyLength)/binSize)+((downstream-1000)/binSize))))
  colnames(control_ranLocMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                          paste0("t", ((upstream/binSize)+1):((upstream+ranLoc_bodyLength)/binSize)),
                                          paste0("d", (((upstream+ranLoc_bodyLength)/binSize)+1):(((upstream+ranLoc_bodyLength)/binSize)+(downstream/binSize))))
}

# Create list of lists in which each element in the enclosing list corresponds to a library
# and the two elements in the nested list correspond to coverage matrices for features and random loci
control_mats <- mclapply(seq_along(controlNames), function(x) {
  list(
       # gene
       control_geneMats[[x]],
       # ranLoc
       control_ranLocMats[[x]],
       # crossover
       control_crossoverMats[[x]],
       # cmt3_hypoCHG_DMRs
       control_cmt3_hypoCHG_DMRMats[[x]],
       # kss_hypoCHG_DMRs
       control_kss_hypoCHG_DMRMats[[x]],
       # cmt3_hypoCHH_DMRs
       control_cmt3_hypoCHH_DMRMats[[x]],
       # kss_hypoCHH_DMRs
       control_kss_hypoCHH_DMRMats[[x]],
       # TEsfs
       control_TEsfMats[[x]],
       # CENAthilas
       control_CENAthilaMats[[x]],
       # nonCENAthilas
       control_nonCENAthilaMats[[x]],
       # CEN180s
       control_CEN180Mats[[x]]
      )
}, mc.cores = length(controlNames))

# Transpose matrix and convert into dataframe
# in which first column is window name
wideDFfeature_list_control <- mclapply(seq_along(control_mats), function(x) {
  lapply(seq_along(control_mats[[x]]), function(y) {
    data.frame(window = colnames(control_mats[[x]][[y]]),
               t(control_mats[[x]][[y]]))
  })
}, mc.cores = length(control_mats))

# Convert into tidy data.frame (long format)
tidyDFfeature_list_control  <- mclapply(seq_along(wideDFfeature_list_control), function(x) {
  lapply(seq_along(control_mats[[x]]), function(y) {
    gather(data  = wideDFfeature_list_control[[x]][[y]],
           key   = feature,
           value = coverage,
           -window)
  }) 
}, mc.cores = length(wideDFfeature_list_control))

# Order levels of factor "window" so that sequential levels
# correspond to sequential windows
for(x in seq_along(tidyDFfeature_list_control)) {
  for(y in seq_along(control_mats[[x]])) {
    tidyDFfeature_list_control[[x]][[y]]$window <- factor(tidyDFfeature_list_control[[x]][[y]]$window,
                                                           levels = as.character(wideDFfeature_list_control[[x]][[y]]$window))
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
  lapply(seq_along(control_mats[[x]]), function(y) {
    data.frame(window = as.character(wideDFfeature_list_control[[x]][[y]]$window),
               n      = tapply(X     = tidyDFfeature_list_control[[x]][[y]]$coverage,
                               INDEX = tidyDFfeature_list_control[[x]][[y]]$window,
                               FUN   = length),
               mean   = tapply(X     = tidyDFfeature_list_control[[x]][[y]]$coverage,
                               INDEX = tidyDFfeature_list_control[[x]][[y]]$window,
                               FUN   = mean,
                               na.rm = TRUE),
               sd     = tapply(X     = tidyDFfeature_list_control[[x]][[y]]$coverage,
                               INDEX = tidyDFfeature_list_control[[x]][[y]]$window,
                               FUN   = sd,
                               na.rm = TRUE))
  })
}, mc.cores = length(tidyDFfeature_list_control))

for(x in seq_along(summaryDFfeature_list_control)) {
  for(y in seq_along(control_mats[[x]])) {
    summaryDFfeature_list_control[[x]][[y]]$window <- factor(summaryDFfeature_list_control[[x]][[y]]$window,
                                                              levels = as.character(wideDFfeature_list_control[[x]][[y]]$window))
    summaryDFfeature_list_control[[x]][[y]]$winNo <- factor(1:dim(summaryDFfeature_list_control[[x]][[y]])[1])
    summaryDFfeature_list_control[[x]][[y]]$sem <- summaryDFfeature_list_control[[x]][[y]]$sd/sqrt(summaryDFfeature_list_control[[x]][[y]]$n-1)
    summaryDFfeature_list_control[[x]][[y]]$CI_lower <- summaryDFfeature_list_control[[x]][[y]]$mean -
      qt(0.975, df = summaryDFfeature_list_control[[x]][[y]]$n-1)*summaryDFfeature_list_control[[x]][[y]]$sem
    summaryDFfeature_list_control[[x]][[y]]$CI_upper <- summaryDFfeature_list_control[[x]][[y]]$mean +
      qt(0.975, df = summaryDFfeature_list_control[[x]][[y]]$n-1)*summaryDFfeature_list_control[[x]][[y]]$sem
  }
}

# Convert list of lists summaryDFfeature_list_control into
# a list of single data.frames containing all meta-profiles for plotting
geneTmp <- lapply(seq_along(summaryDFfeature_list_control), function(x) {
  summaryDFfeature_list_control[[x]][[1]]
})
ranLocTmp <- lapply(seq_along(summaryDFfeature_list_control), function(x) {
  summaryDFfeature_list_control[[x]][[2]]
})
crossoverTmp <- lapply(seq_along(summaryDFfeature_list_control), function(x) {
  summaryDFfeature_list_control[[x]][[3]]
})
cmt3_hypoCHG_DMRTmp <- lapply(seq_along(summaryDFfeature_list_control), function(x) {
  summaryDFfeature_list_control[[x]][[4]]
})
kss_hypoCHG_DMRTmp <- lapply(seq_along(summaryDFfeature_list_control), function(x) {
  summaryDFfeature_list_control[[x]][[5]]
})
cmt3_hypoCHH_DMRTmp <- lapply(seq_along(summaryDFfeature_list_control), function(x) {
  summaryDFfeature_list_control[[x]][[6]]
})
kss_hypoCHH_DMRTmp <- lapply(seq_along(summaryDFfeature_list_control), function(x) {
  summaryDFfeature_list_control[[x]][[7]]
})
TEsfTmp <- lapply(seq_along(summaryDFfeature_list_control), function(x) {
  summaryDFfeature_list_control[[x]][[8]]
})
CENAthilaTmp <- lapply(seq_along(summaryDFfeature_list_control), function(x) {
  summaryDFfeature_list_control[[x]][[9]]
})
nonCENAthilaTmp <- lapply(seq_along(summaryDFfeature_list_control), function(x) {
  summaryDFfeature_list_control[[x]][[10]]
})
CEN180Tmp <- lapply(seq_along(summaryDFfeature_list_control), function(x) {
  summaryDFfeature_list_control[[x]][[11]]
})

names(geneTmp) <- controlNamesPlot
names(ranLocTmp) <- controlNamesPlot
names(crossoverTmp) <- controlNamesPlot
names(cmt3_hypoCHG_DMRTmp) <- controlNamesPlot
names(kss_hypoCHG_DMRTmp) <- controlNamesPlot
names(cmt3_hypoCHH_DMRTmp) <- controlNamesPlot
names(kss_hypoCHH_DMRTmp) <- controlNamesPlot
names(TEsfTmp) <- controlNamesPlot
names(CENAthilaTmp) <- controlNamesPlot
names(nonCENAthilaTmp) <- controlNamesPlot
names(CEN180Tmp) <- controlNamesPlot
summaryDFfeature_control <- list(
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
for(x in seq_along(summaryDFfeature_control)) {
  summaryDFfeature_control[[x]]$libName <- factor(summaryDFfeature_control[[x]]$libName,
                                                   levels = controlNamesPlot)
}

# Define y-axis limits
ymin_control <- min(c(summaryDFfeature_control[[1]]$CI_lower,
                       summaryDFfeature_control[[2]]$CI_lower,
                       summaryDFfeature_control[[3]]$CI_lower,
                       summaryDFfeature_control[[4]]$CI_lower,
                       summaryDFfeature_control[[5]]$CI_lower,
                       summaryDFfeature_control[[6]]$CI_lower,
                       summaryDFfeature_control[[7]]$CI_lower,
                       summaryDFfeature_control[[8]]$CI_lower,
                       summaryDFfeature_control[[9]]$CI_lower,
                       summaryDFfeature_control[[10]]$CI_lower,
                       summaryDFfeature_control[[11]]$CI_lower),
                 na.rm = T)
ymax_control <- max(c(summaryDFfeature_control[[1]]$CI_upper,
                       summaryDFfeature_control[[2]]$CI_upper,
                       summaryDFfeature_control[[3]]$CI_upper,
                       summaryDFfeature_control[[4]]$CI_upper,
                       summaryDFfeature_control[[5]]$CI_upper,
                       summaryDFfeature_control[[6]]$CI_upper,
                       summaryDFfeature_control[[7]]$CI_upper,
                       summaryDFfeature_control[[8]]$CI_upper,
                       summaryDFfeature_control[[9]]$CI_upper,
                       summaryDFfeature_control[[10]]$CI_upper,
                       summaryDFfeature_control[[11]]$CI_upper),
                 na.rm = T)

# Define legend labels
legendLabs <- lapply(seq_along(controlNamesPlot), function(x) {
  grobTree(textGrob(bquote(.(controlNamesPlot[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = controlColours[x], fontsize = 18)))
})


# Plot average profiles with 95% CI ribbon
## gene
summaryDFfeature <- summaryDFfeature_control[[1]]
ggObj1_combined_control <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = controlColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = controlColours) +
scale_y_continuous(limits = c(ymin_control, ymax_control),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_control[[1]])[1]/length(controlNames))-(downstream/binSize),
                            dim(summaryDFfeature_control[[1]])[1]/length(controlNames)),
                 labels = c(paste0("-", flankName),
                            geneStartLab,
                            geneEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_control[[1]])[1]/length(controlNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Norm. coverage")) +
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
summaryDFfeature <- summaryDFfeature_control[[2]]
ggObj2_combined_control <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = controlColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = controlColours) +
scale_y_continuous(limits = c(ymin_control, ymax_control),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_control[[2]])[1]/length(controlNames))-(downstream/binSize),
                            dim(summaryDFfeature_control[[2]])[1]/length(controlNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_control[[2]])[1]/length(controlNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Norm. coverage")) +
annotation_custom(legendLabs[[1]]) +
annotation_custom(legendLabs[[2]]) +
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
summaryDFfeature <- summaryDFfeature_control[[3]]
ggObj3_combined_control <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = controlColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = controlColours) +
scale_y_continuous(limits = c(ymin_control, ymax_control),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_control[[3]])[1]/length(controlNames))-(downstream/binSize),
                            dim(summaryDFfeature_control[[3]])[1]/length(controlNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_control[[3]])[1]/length(controlNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Norm. coverage")) +
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
summaryDFfeature <- summaryDFfeature_control[[4]]
ggObj4_combined_control <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = controlColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = controlColours) +
scale_y_continuous(limits = c(ymin_control, ymax_control),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_control[[4]])[1]/length(controlNames))-(downstream/binSize),
                            dim(summaryDFfeature_control[[4]])[1]/length(controlNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_control[[4]])[1]/length(controlNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Norm. coverage")) +
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
summaryDFfeature <- summaryDFfeature_control[[5]]
ggObj5_combined_control <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = controlColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = controlColours) +
scale_y_continuous(limits = c(ymin_control, ymax_control),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_control[[5]])[1]/length(controlNames))-(downstream/binSize),
                            dim(summaryDFfeature_control[[5]])[1]/length(controlNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_control[[5]])[1]/length(controlNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Norm. coverage")) +
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
summaryDFfeature <- summaryDFfeature_control[[6]]
ggObj6_combined_control <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = controlColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = controlColours) +
scale_y_continuous(limits = c(ymin_control, ymax_control),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_control[[6]])[1]/length(controlNames))-(downstream/binSize),
                            dim(summaryDFfeature_control[[6]])[1]/length(controlNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_control[[6]])[1]/length(controlNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Norm. coverage")) +
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
summaryDFfeature <- summaryDFfeature_control[[7]]
ggObj7_combined_control <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = controlColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = controlColours) +
scale_y_continuous(limits = c(ymin_control, ymax_control),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_control[[7]])[1]/length(controlNames))-(downstream/binSize),
                            dim(summaryDFfeature_control[[7]])[1]/length(controlNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_control[[7]])[1]/length(controlNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Norm. coverage")) +
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
summaryDFfeature <- summaryDFfeature_control[[8]]
ggObj8_combined_control <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = controlColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = controlColours) +
scale_y_continuous(limits = c(ymin_control, ymax_control),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_control[[8]])[1]/length(controlNames))-(downstream/binSize),
                            dim(summaryDFfeature_control[[8]])[1]/length(controlNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_control[[8]])[1]/length(controlNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Norm. coverage")) +
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
summaryDFfeature <- summaryDFfeature_control[[9]]
ggObj9_combined_control <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = controlColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = controlColours) +
scale_y_continuous(limits = c(ymin_control, ymax_control),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_control[[9]])[1]/length(controlNames))-(downstream/binSize),
                            dim(summaryDFfeature_control[[9]])[1]/length(controlNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_control[[9]])[1]/length(controlNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Norm. coverage")) +
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
summaryDFfeature <- summaryDFfeature_control[[10]]
ggObj10_combined_control <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = controlColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = controlColours) +
scale_y_continuous(limits = c(ymin_control, ymax_control),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_control[[10]])[1]/length(controlNames))-(downstream/binSize),
                            dim(summaryDFfeature_control[[10]])[1]/length(controlNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_control[[10]])[1]/length(controlNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Norm. coverage")) +
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
summaryDFfeature <- summaryDFfeature_control[[11]]
ggObj11_combined_control <- ggplot(data = summaryDFfeature,
                                    mapping = aes(x = winNo,
                                                  y = mean,
                                                  group = libName)
                                   ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = controlColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = controlColours) +
scale_y_continuous(limits = c(ymin_control, ymax_control),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            ((upstream-1000)/binSize)+1,
                            (dim(summaryDFfeature_control[[11]])[1]/length(controlNames))-((downstream-1000)/binSize),
                            dim(summaryDFfeature_control[[11]])[1]/length(controlNames)),
                 labels = c(paste0("-", "1kb"),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", "1kb"))) +
geom_vline(xintercept = c(((upstream-1000)/binSize)+1,
                          (dim(summaryDFfeature_control[[11]])[1]/length(controlNames))-((downstream-1000)/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Norm. coverage")) +
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
                                              ggObj1_combined_control,
                                              ggObj2_combined_control,
                                              ggObj3_combined_control,
                                              ggObj4_combined_control,
                                              ggObj5_combined_control,
                                              ggObj6_combined_control,
                                              ggObj7_combined_control,
                                              ggObj8_combined_control,
                                              ggObj9_combined_control,
                                              ggObj10_combined_control,
                                              ggObj11_combined_control
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
              "control_",
              paste0(controlNames, collapse = "_"),
              "_avgProfiles_around",
              "_genes_ranLoc_COs_hypoDMRs4sets_", TEsf, "_CENAthila_nonCENAthila_CEN180_in_t2t-col.20210610_",
              paste0(chrName, collapse = "_"), "_", align, ".pdf"),
       plot = ggObjGA_combined,
       height = 6.5, width = 7*11, limitsize = FALSE)


