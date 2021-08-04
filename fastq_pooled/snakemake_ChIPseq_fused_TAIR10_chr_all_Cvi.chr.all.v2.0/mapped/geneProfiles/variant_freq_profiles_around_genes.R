#!/applications/R/R-4.0.0/bin/Rscript

# author: Andy Tock
# contact: ajt200@cam.ac.uk
# date: 04.08.2021

# Profile varType variant frequency around genes and random loci

# Usage:
# /applications/R/R-4.0.0/bin/Rscript variant_freq_profiles_around_genes.R SNP_INDEL '/data/public_data/arabidopsis/MPIPZ_Jiao_Schneeberger_2020_NatCommun/Cvi' unique 'Chr1,Chr2,Chr3,Chr4,Chr5' 2200 2000 2kb 10 10bp genomewide fused_TAIR10_chr_all_Cvi.chr.all.v2.0

#varType <- "SNP_INDEL"
#dirName <- "/data/public_data/arabidopsis/MPIPZ_Jiao_Schneeberger_2020_NatCommun/Cvi"
#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split = ","))
#regionBodyLength <- 2200
#upstream <- 2000
#downstream <- 2000
#flankName <- "2kb"
#binSize <- 10
#binName <- "10bp"
#genomeRegion <- "genomewide"
#refbase <- "fused_TAIR10_chr_all_Cvi.chr.all.v2.0"

args <- commandArgs(trailingOnly = T)
varType <- args[1]
dirName <- args[2]
chrName <- unlist(strsplit(args[3],
                               split = ","))
regionBodyLength <- as.numeric(args[4])
upstream <- as.numeric(args[5])
downstream <- as.numeric(args[5])
flankName <- args[6]
binSize <- as.numeric(args[7])
binName <- args[8]
genomeRegion <- args[9]
refbase <- args[10]

options(stringsAsFactors = F)
library(data.table)
library(EnrichedHeatmap)
library(parallel)

matDir <- paste0("matrices_smoothed/")
system(paste0("[ -d ", matDir, " ] || mkdir ", matDir))

# Genomic definitions
fai <- read.table(paste0(dirName, "/", refbase, "/", refbase, ".fa.fai"), header = F)
chrs <- fai$V1[which(gsub(pattern = ".+_", replacement = "", x = fai$V1) %in% chrName)]
chrLens <- fai$V2[which(gsub(pattern = ".+_", replacement = "", x = fai$V1) %in% chrName)]
# Pericentromeres in TAIR10 as defined in Ziolkowski et al. 2017 Genes Dev. 31: Table S26 
pericen <- read.table(paste0(dirName, "/", refbase, "/", refbase, ".fa.pericentromeres"), header = T)
pericen <- pericen[which(gsub(pattern = ".+_", replacement = "", x = pericen$chr) %in% chrName),]

# Define genomic regions to be analysed (genomeRegionGR)
if(genomeRegion == "arm") {
  genomeRegionGR <- GRanges(seqnames = rep(chrs, 2),
                            ranges = IRanges(start = c(rep(1, length(chrs)),
                                                       pericen$end+1),
                                             end = c(pericen$start-1,
                                                     chrLens)),
                            strand = "*")
  genomeRegionGR <- genomeRegionGR[which(gsub(pattern = ".+_", replacement = "", x = seqnames(genomeRegionGR)@values) %in% chrName)]
} else if(genomeRegion == "peri") {
  genomeRegionGR <- GRanges(seqnames = chrs,
                            ranges = IRanges(start = pericen$start,
                                             end = pericen$end),
                            strand = "*")
  genomeRegionGR <- genomeRegionGR[which(gsub(pattern = ".+_", replacement = "", x = seqnames(genomeRegionGR)@values) %in% chrName)]
} else if(genomeRegion == "genomewide") {
  genomeRegionGR <- GRanges(seqnames = chrs,
                            ranges = IRanges(start = rep(1, length(chrs)),
                                             end = chrLens),
                            strand = "*")
  genomeRegionGR <- genomeRegionGR[which(gsub(pattern = ".+_", replacement = "", x = seqnames(genomeRegionGR)@values) %in% chrName)]
} else {
  stop("genomeRegion is not arm, peri or genomewide")
}

# Define genomic regions to be masked from analyses (genomeMaskGR)
if(genomeRegion == "arm") {
  genomeMaskGR <- GRanges(seqnames = chrs,
                          ranges = IRanges(start = pericen$start,
                                           end = pericen$end),
                          strand = "*")
  genomeMaskGR <- genomeMaskGR[which(gsub(pattern = ".+_", replacement = "", x = seqnames(genomeMaskGR)@values) %in% chrName)]
} else if(genomeRegion == "peri") {
  genomeMaskGR <- GRanges(seqnames = rep(chrs, 2),
                          ranges = IRanges(start = c(rep(1, length(chrs)),
                                                     pericen$end+1),
                                           end = c(pericen$start-1,
                                                   chrLens)),
                          strand = "*")
  genomeMaskGR <- genomeMaskGR[which(gsub(pattern = ".+_", replacement = "", x = seqnames(genomeMaskGR)@values) %in% chrName)]
} else if(genomeRegion == "genomewide") {
  genomeMaskGR <- GRanges()
  genomeMaskGR <- genomeMaskGR[which(gsub(pattern = ".+_", replacement = "", x = seqnames(genomeMaskGR)@values) %in% chrName)]
} else {
  stop("genomeRegion is not arm, peri or genomewide")
}

# Load featuresAcc1_ortho_DF_bed
featuresAcc1_ortho_DF_bed <- read.table(
  paste0("quantiles/", paste0(chrName, collapse = "_"), "/",
         "quantiles_", genomeRegion, "_by_SNP_INDEL_freq_in_bodies/",
         "featuresAcc1_ortho_6quantiles",
         "_", genomeRegion, "_by_SNP_INDEL_freq_in_bodies",
         "_of_Acc1_Chr_genes_in_", refbase, "_",
         paste0(chrName, collapse = "_"), ".bed"),
  header = FALSE, sep = "\t")
colnames(featuresAcc1_ortho_DF_bed) <- c("chr", "start", "end", "ID", "score", "strand")

# Convert into GRanges
# Convert featuresAcc1_ortho_DF_bed into GRanges
# Note addition of 1 to 0-based BED start coordinates
featuresAcc1_orthoGR <- GRanges(seqnames = featuresAcc1_ortho_DF_bed$chr,
                                ranges = IRanges(start = featuresAcc1_ortho_DF_bed$start+1,
                                                 end = featuresAcc1_ortho_DF_bed$end),
                                strand = featuresAcc1_ortho_DF_bed$strand,
                                featureID = featuresAcc1_ortho_DF_bed$ID)

# Load ranLocAcc1_ortho_DF_bed
ranLocAcc1_ortho_DF_bed <- read.table(
  paste0("quantiles/", paste0(chrName, collapse = "_"), "/",
         "quantiles_", genomeRegion, "_by_SNP_INDEL_freq_in_bodies/",
         "featuresAcc1_ortho_6quantiles",
         "_", genomeRegion, "_by_SNP_INDEL_freq_in_bodies",
         "_of_Acc1_Chr_genes_in_", refbase, "_",
         paste0(chrName, collapse = "_"), "_ranLocAcc1.bed"),
  header = FALSE, sep = "\t")
colnames(ranLocAcc1_ortho_DF_bed) <- c("chr", "start", "end", "ID", "score", "strand")

# Convert into GRanges
# Convert ranLocAcc1_ortho_DF_bed into GRanges
# Note addition of 1 to 0-based BED start coordinates
ranLocAcc1_orthoGR <- GRanges(seqnames = ranLocAcc1_ortho_DF_bed$chr,
                              ranges = IRanges(start = ranLocAcc1_ortho_DF_bed$start+1,
                                               end = ranLocAcc1_ortho_DF_bed$end),
                              strand = ranLocAcc1_ortho_DF_bed$strand,
                              featureID = ranLocAcc1_ortho_DF_bed$ID)


# Load table of varAcc1
varAcc1 <- fread(input = paste0(dirName, "/",
                                substr(x = refbase, start = 22, stop = 24), ".syri.out.gz"),
                  header = FALSE)
varAcc1 <- data.frame(varAcc1)
colnames(varAcc1) <- c("ref_chr", "ref_start", "ref_end", "ref_seq", "que_seq",
                       "que_chr", "que_start", "que_end", "uniqueID", "parentID",
                       "type", "copy_status")
varAcc1$ref_start <- as.numeric(varAcc1$ref_start)
varAcc1$ref_end <- as.numeric(varAcc1$ref_end)
varAcc1$que_start <- as.numeric(varAcc1$que_start)
varAcc1$que_end <- as.numeric(varAcc1$que_end)
varAcc1$ref_chr <- gsub(pattern = "Chr", replacement = "Col_Chr",
                        x = varAcc1$ref_chr)
varAcc1$que_chr <- gsub(pattern = "Chr", replacement = paste0(substr(x = refbase, start = 22, stop = 24), "_Chr"),
                        x = varAcc1$que_chr)

# Remove varAcc1 that are listed more than once relative to the reference genome (Columns 1-5)
varAcc1 <- varAcc1[-which(duplicated(varAcc1[,1:5])),]

# Extract rows corresponding to varAcc1 of type varType
if(varType == "SNP_INDEL") {
  varAcc1 <- varAcc1[which(varAcc1$type %in% c("SNP", "INS", "DEL")),]
} else if(!varType %in% c("all_variants")) {
  varAcc1 <- varAcc1[which(varAcc1$type == varType),]
}

# Flip variant coordinates in wrong orientation for GRanges conversion
varAcc1_end_start <- varAcc1[which(varAcc1$ref_end < varAcc1$ref_start),]
varAcc1 <- varAcc1[-which(varAcc1$ref_end < varAcc1$ref_start),]
varAcc1_start_end <- varAcc1_end_start
varAcc1_start_end$ref_start <- varAcc1_end_start$ref_end
varAcc1_start_end$ref_end <- varAcc1_end_start$ref_start
varAcc1 <- rbind(varAcc1,
                 varAcc1_start_end)

# Order varAcc1 according to chromosome and coordinates
varAcc1 <- varAcc1[ with(varAcc1, order(ref_chr, ref_start, ref_end)), ]

print(paste0("Number of varAcc1 of type ", varType, ":"))
print(nrow(varAcc1))

varAcc1 <- data.frame(varAcc1,
                      greptype = NA)

varAcc1[varAcc1$type == "SNP",]$greptype <- "A_SNP"
varAcc1[varAcc1$type == "INS",]$greptype <- "A_insertion"
varAcc1[varAcc1$type == "DEL",]$greptype <- "A_deletion"

# Convert into GRanges
varAcc1GR <- GRanges(seqnames = varAcc1$ref_chr,
                     ranges = IRanges(start = varAcc1$ref_start,
                                      end = varAcc1$ref_end),
                     strand = "*",
                     coverage = rep(1, nrow(varAcc1)))

# Subset varAcc1 by class (grep-ing for "A" within "greptype" field returns all varAcc1)

varAcc1class <- c(
                  "A",
                  "SNP",
                  "insertion",
                  "deletion",
                  "tion"
                 )
varAcc1ListGR <- mclapply(seq_along(varAcc1class), function(x) {
  classvarAcc1 <- varAcc1[grep(varAcc1class[x], varAcc1$greptype),]
  GRanges(seqnames = classvarAcc1$chr,
          ranges = IRanges(start = classvarAcc1$pos,
                           end = classvarAcc1$pos),
          strand = "*",
          coverage = rep(1, dim(classvarAcc1)[1]))
}, mc.cores = length(varAcc1class))
stopifnot(length(varAcc1ListGR[[5]]) ==
          length(varAcc1ListGR[[4]]) + length(varAcc1ListGR[[3]]))

## transition
varAcc1_transition <- varAcc1[(varAcc1$ref == "A" | varAcc1$ref == "G") & (varAcc1$alt == "G" | varAcc1$alt == "A") |
                        (varAcc1$ref == "C" | varAcc1$ref == "T") & (varAcc1$alt == "T" | varAcc1$alt == "C"),]
varAcc1_transition_GR <- GRanges(seqnames = varAcc1_transition$chr,
                              ranges = IRanges(start = varAcc1_transition$pos,
                                               end = varAcc1_transition$pos),
                              strand = "*",
                              coverage = rep(1, dim(varAcc1_transition)[1]))

## transversion
varAcc1_transversion <- varAcc1[(varAcc1$ref == "A" | varAcc1$ref == "G") & (varAcc1$alt == "C" | varAcc1$alt == "T") |
                          (varAcc1$ref == "C" | varAcc1$ref == "T") & (varAcc1$alt == "A" | varAcc1$alt == "G"),]
stopifnot((dim(varAcc1_transition)[1] +
          dim(varAcc1_transversion)[1]) ==
          dim(varAcc1[varAcc1$greptype == "A_SNP",])[1])
varAcc1_transversion_GR <- GRanges(seqnames = varAcc1_transversion$chr,
                                ranges = IRanges(start = varAcc1_transversion$pos,
                                                 end = varAcc1_transversion$pos),
                                strand = "*",
                                coverage = rep(1, dim(varAcc1_transversion)[1]))

# Add transitions and transversions GRanges to varAcc1ListGR
varAcc1ListGR <- c(varAcc1ListGR,
                varAcc1_transition_GR,
                varAcc1_transversion_GR)
varAcc1classNames <- c(
                   "all",
                   "SNP",
                   "insertion",
                   "deletion",
                   "indel",
                   "transition",
                   "transversion"
                  )

if(featureName %in% c("CEN180", "CENgapAllAthila")) {
  # Define matrix and column mean outfiles
  outDF <- lapply(seq_along(varAcc1classNames), function(x) {
    list(paste0(matDir, varAcc1classNames[x],
                "_CEN180_consensus_variants_MappedOn_t2t-col.20210610_around_", featureName, "_in_", chrName,
                "_matrix_bin", binName, "_flank", flankName, ".tab"),
         paste0(matDir, varAcc1classNames[x],
                "_CEN180_consensus_variants_MappedOn_t2t-col.20210610_around_", featureName, "_in_", chrName,
                "_CENranLoc_matrix_bin", binName, "_flank", flankName, ".tab"))
  })
  outDFcolMeans <- lapply(seq_along(varAcc1classNames), function(x) {
    list(paste0(matDir, varAcc1classNames[x],
                "_CEN180_consensus_variants_MappedOn_t2t-col.20210610_around_", featureName, "_in_", chrName,
                "_matrix_bin", binName, "_flank", flankName, "_colMeans.tab"),
         paste0(matDir, varAcc1classNames[x],
                "_CEN180_consensus_variants_MappedOn_t2t-col.20210610_around_", featureName, "_in_", chrName,
                "_CENranLoc_matrix_bin", binName, "_flank", flankName, "_colMeans.tab"))
  })
  
  # Function to create variant frequency matrices for
  # feature loci and random loci (incl. flanking regions)
  # and to calculate mean profiles across all feature loci and random loci
  covMatrix <- function(signal,
                        feature,
                        ranLoc,
                        featureSize,
                        flankSize,
                        winSize,
                        outDF,
                        outDFcolMeans) {
    # feature loci
    set.seed(2840)
    feature_smoothed <- normalizeToMatrix(signal = signal,
                                          target = feature,
                                          value_column = "coverage",
                                          extend = flankSize,
                                          mean_mode = "w0",
                                          w = winSize,
                                          background = 0,
                                          smooth = TRUE,
                                          include_target = TRUE,
                                          target_ratio = featureSize/(featureSize+(flankSize*2)))
    print("feature_smoothed")
    print(feature_smoothed)
    print("feature_smoothed rows = ")
    print(length(feature_smoothed)/round((featureSize/winSize)+((flankSize*2)/winSize)))
    feature_smoothed_DF <- data.frame(feature_smoothed)
    feature_smoothed_DF_colMeans <- as.vector(colMeans(feature_smoothed_DF,
                                                       na.rm = T))
    write.table(feature_smoothed_DF,
                file = outDF[[1]],
                quote = F, sep = "\t", row.names = F, col.names = T)
    write.table(feature_smoothed_DF_colMeans,
                file = outDFcolMeans[[1]],
                quote = F, sep = "\t", row.names = F, col.names = T)
  
    # random loci
    set.seed(8472)
    ranLoc_smoothed <- normalizeToMatrix(signal = signal,
                                         target = ranLoc,
                                         value_column = "coverage",
                                         extend = flankSize,
                                         mean_mode = "w0",
                                         w = winSize,
                                         background = 0,
                                         smooth = TRUE,
                                         include_target = TRUE,
                                         target_ratio = featureSize/(featureSize+(flankSize*2)))
    print("ranLoc_smoothed")
    print(ranLoc_smoothed)
    print("ranLoc_smoothed rows = ")
    print(length(ranLoc_smoothed)/round((featureSize/winSize)+((flankSize*2)/winSize)))
    ranLoc_smoothed_DF <- data.frame(ranLoc_smoothed)
    ranLoc_smoothed_DF_colMeans <- as.vector(colMeans(ranLoc_smoothed_DF,
                                                      na.rm = T))
    write.table(ranLoc_smoothed_DF,
                file = outDF[[2]],
                quote = F, sep = "\t", row.names = F, col.names = T)
    write.table(ranLoc_smoothed_DF_colMeans,
                file = outDFcolMeans[[2]],
                quote = F, sep = "\t", row.names = F, col.names = T)
  }
  
  # Run covMatrix() function on each feature GRanges object to obtain matrices
  # containing normalised feature density values around target and random loci
  mclapply(seq_along(varAcc1classNames), function(x) {
    covMatrix(signal = varAcc1ListGR[[x]],
              feature = featuresGR,
              ranLoc = ranLocGR,
              featureSize = regionBodyLength,
              flankSize = upstream,
              winSize = binSize,
              outDF = outDF[[x]],
              outDFcolMeans = outDFcolMeans[[x]])
    print(paste0(varAcc1classNames[x],
                 " CEN180 consensus variants frequency around ", featureName, " in ", chrName,
                 " profile calculation complete"))
  }, mc.cores = length(varAcc1classNames))
} else {
  # Define matrix and column mean outfiles
  outDF <- lapply(seq_along(varAcc1classNames), function(x) {
    list(paste0(matDir, varAcc1classNames[x],
                "_CEN180_consensus_variants_MappedOn_t2t-col.20210610_around_", featureName, "_in_", chrName,
                "_matrix_bin", binName, "_flank", flankName, ".tab"))
  })
  outDFcolMeans <- lapply(seq_along(varAcc1classNames), function(x) {
    list(paste0(matDir, varAcc1classNames[x],
                "_CEN180_consensus_variants_MappedOn_t2t-col.20210610_around_", featureName, "_in_", chrName,
                "_matrix_bin", binName, "_flank", flankName, "_colMeans.tab"))
  })
  
  # Function to create variant frequency matrices for
  # feature loci and random loci (incl. flanking regions)
  # and to calculate mean profiles across all feature loci and random loci
  covMatrix <- function(signal,
                        feature,
                        featureSize,
                        flankSize,
                        winSize,
                        outDF,
                        outDFcolMeans) {
    # feature loci
    set.seed(2840)
    feature_smoothed <- normalizeToMatrix(signal = signal,
                                          target = feature,
                                          value_column = "coverage",
                                          extend = flankSize,
                                          mean_mode = "w0",
                                          w = winSize,
                                          background = 0,
                                          smooth = TRUE,
                                          include_target = TRUE,
                                          target_ratio = featureSize/(featureSize+(flankSize*2)))
    print("feature_smoothed")
    print(feature_smoothed)
    print("feature_smoothed rows = ")
    print(length(feature_smoothed)/round((featureSize/winSize)+((flankSize*2)/winSize)))
    feature_smoothed_DF <- data.frame(feature_smoothed)
    feature_smoothed_DF_colMeans <- as.vector(colMeans(feature_smoothed_DF,
                                                       na.rm = T))
    write.table(feature_smoothed_DF,
                file = outDF[[1]],
                quote = F, sep = "\t", row.names = F, col.names = T)
    write.table(feature_smoothed_DF_colMeans,
                file = outDFcolMeans[[1]],
                quote = F, sep = "\t", row.names = F, col.names = T)
  }
  
  # Run covMatrix() function on each feature GRanges object to obtain matrices
  # containing normalised feature density values around target and random loci
  mclapply(seq_along(varAcc1classNames), function(x) {
    covMatrix(signal = varAcc1ListGR[[x]],
              feature = featuresGR,
              featureSize = regionBodyLength,
              flankSize = upstream,
              winSize = binSize,
              outDF = outDF[[x]],
              outDFcolMeans = outDFcolMeans[[x]])
    print(paste0(varAcc1classNames[x],
                 " CEN180 consensus variants frequency around ", featureName, " in ", chrName,
                 " profile calculation complete"))
  }, mc.cores = length(varAcc1classNames))
}
