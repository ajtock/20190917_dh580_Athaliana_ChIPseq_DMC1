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

## Genomic definitions
#fai <- read.table(paste0(dirName, "/", refbase, "/", refbase, ".fa.fai"), header = F)
#chrs <- fai$V1[which(gsub(pattern = ".+_", replacement = "", x = fai$V1) %in% chrName)]
#chrLens <- fai$V2[which(gsub(pattern = ".+_", replacement = "", x = fai$V1) %in% chrName)]
## Pericentromeres in TAIR10 as defined in Ziolkowski et al. 2017 Genes Dev. 31: Table S26 
#pericen <- read.table(paste0(dirName, "/", refbase, "/", refbase, ".fa.pericentromeres"), header = T)
#pericen <- pericen[which(gsub(pattern = ".+_", replacement = "", x = pericen$chr) %in% chrName),]
#
## Define genomic regions to be analysed (genomeRegionGR)
#if(genomeRegion == "arm") {
#  genomeRegionGR <- GRanges(seqnames = rep(chrs, 2),
#                            ranges = IRanges(start = c(rep(1, length(chrs)),
#                                                       pericen$end+1),
#                                             end = c(pericen$start-1,
#                                                     chrLens)),
#                            strand = "*")
#  genomeRegionGR <- genomeRegionGR[which(gsub(pattern = ".+_", replacement = "", x = seqnames(genomeRegionGR)@values) %in% chrName)]
#} else if(genomeRegion == "peri") {
#  genomeRegionGR <- GRanges(seqnames = chrs,
#                            ranges = IRanges(start = pericen$start,
#                                             end = pericen$end),
#                            strand = "*")
#  genomeRegionGR <- genomeRegionGR[which(gsub(pattern = ".+_", replacement = "", x = seqnames(genomeRegionGR)@values) %in% chrName)]
#} else if(genomeRegion == "genomewide") {
#  genomeRegionGR <- GRanges(seqnames = chrs,
#                            ranges = IRanges(start = rep(1, length(chrs)),
#                                             end = chrLens),
#                            strand = "*")
#  genomeRegionGR <- genomeRegionGR[which(gsub(pattern = ".+_", replacement = "", x = seqnames(genomeRegionGR)@values) %in% chrName)]
#} else {
#  stop("genomeRegion is not arm, peri or genomewide")
#}
#
## Define genomic regions to be masked from analyses (genomeMaskGR)
#if(genomeRegion == "arm") {
#  genomeMaskGR <- GRanges(seqnames = chrs,
#                          ranges = IRanges(start = pericen$start,
#                                           end = pericen$end),
#                          strand = "*")
#  genomeMaskGR <- genomeMaskGR[which(gsub(pattern = ".+_", replacement = "", x = seqnames(genomeMaskGR)@values) %in% chrName)]
#} else if(genomeRegion == "peri") {
#  genomeMaskGR <- GRanges(seqnames = rep(chrs, 2),
#                          ranges = IRanges(start = c(rep(1, length(chrs)),
#                                                     pericen$end+1),
#                                           end = c(pericen$start-1,
#                                                   chrLens)),
#                          strand = "*")
#  genomeMaskGR <- genomeMaskGR[which(gsub(pattern = ".+_", replacement = "", x = seqnames(genomeMaskGR)@values) %in% chrName)]
#} else if(genomeRegion == "genomewide") {
#  genomeMaskGR <- GRanges()
#  genomeMaskGR <- genomeMaskGR[which(gsub(pattern = ".+_", replacement = "", x = seqnames(genomeMaskGR)@values) %in% chrName)]
#} else {
#  stop("genomeRegion is not arm, peri or genomewide")
#}

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

# Load featuresAcc2_ortho_DF_bed
featuresAcc2_ortho_DF_bed <- read.table(
  paste0("quantiles/", paste0(chrName, collapse = "_"), "/",
         "quantiles_", genomeRegion, "_by_SNP_INDEL_freq_in_bodies/",
         "featuresAcc2_ortho_6quantiles",
         "_", genomeRegion, "_by_SNP_INDEL_freq_in_bodies",
         "_of_Acc2_Chr_genes_in_", refbase, "_",
         paste0(chrName, collapse = "_"), ".bed"),
  header = FALSE, sep = "\t")
colnames(featuresAcc2_ortho_DF_bed) <- c("chr", "start", "end", "ID", "score", "strand")

# Convert into GRanges
# Convert featuresAcc2_ortho_DF_bed into GRanges
# Note addition of 1 to 0-based BED start coordinates
featuresAcc2_orthoGR <- GRanges(seqnames = featuresAcc2_ortho_DF_bed$chr,
                                ranges = IRanges(start = featuresAcc2_ortho_DF_bed$start+1,
                                                 end = featuresAcc2_ortho_DF_bed$end),
                                strand = featuresAcc2_ortho_DF_bed$strand,
                                featureID = featuresAcc2_ortho_DF_bed$ID)

# Load ranLocAcc2_ortho_DF_bed
ranLocAcc2_ortho_DF_bed <- read.table(
  paste0("quantiles/", paste0(chrName, collapse = "_"), "/",
         "quantiles_", genomeRegion, "_by_SNP_INDEL_freq_in_bodies/",
         "featuresAcc2_ortho_6quantiles",
         "_", genomeRegion, "_by_SNP_INDEL_freq_in_bodies",
         "_of_Acc2_Chr_genes_in_", refbase, "_",
         paste0(chrName, collapse = "_"), "_ranLocAcc2.bed"),
  header = FALSE, sep = "\t")
colnames(ranLocAcc2_ortho_DF_bed) <- c("chr", "start", "end", "ID", "score", "strand")

# Convert into GRanges
# Convert ranLocAcc2_ortho_DF_bed into GRanges
# Note addition of 1 to 0-based BED start coordinates
ranLocAcc2_orthoGR <- GRanges(seqnames = ranLocAcc2_ortho_DF_bed$chr,
                              ranges = IRanges(start = ranLocAcc2_ortho_DF_bed$start+1,
                                               end = ranLocAcc2_ortho_DF_bed$end),
                              strand = ranLocAcc2_ortho_DF_bed$strand,
                              featureID = ranLocAcc2_ortho_DF_bed$ID)



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
                     coverage = rep(1, dim(varAcc1)[1]))

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
  GRanges(seqnames = classvarAcc1$ref_chr,
          ranges = IRanges(start = classvarAcc1$ref_start,
                           end = classvarAcc1$ref_end),
          strand = "*",
          coverage = rep(1, dim(classvarAcc1)[1]))
}, mc.cores = length(varAcc1class))
stopifnot(length(varAcc1ListGR[[5]]) ==
          length(varAcc1ListGR[[4]]) + length(varAcc1ListGR[[3]]))

## transition
varAcc1_transition <- varAcc1[(varAcc1$ref_seq == "A" | varAcc1$ref_seq == "G") & (varAcc1$que_seq == "G" | varAcc1$que_seq == "A") |
                              (varAcc1$ref_seq == "C" | varAcc1$ref_seq == "T") & (varAcc1$que_seq == "T" | varAcc1$que_seq == "C"),]
varAcc1_transition_GR <- GRanges(seqnames = varAcc1_transition$ref_chr,
                                 ranges = IRanges(start = varAcc1_transition$ref_start,
                                                  end = varAcc1_transition$ref_end),
                                 strand = "*",
                                 coverage = rep(1, dim(varAcc1_transition)[1]))

## transversion
varAcc1_transversion <- varAcc1[(varAcc1$ref_seq == "A" | varAcc1$ref_seq == "G") & (varAcc1$que_seq == "C" | varAcc1$que_seq == "T") |
                                (varAcc1$ref_seq == "C" | varAcc1$ref_seq == "T") & (varAcc1$que_seq == "A" | varAcc1$que_seq == "G"),]
## Below sanity check doesn't work because due to use of IUPAC ambiguity codes for some alleles
## (e.g., found one in varAcc1[varAcc1$greptype == "A_SNP",][117920,] )
#stopifnot((dim(varAcc1_transition)[1] +
#           dim(varAcc1_transversion)[1]) ==
#           dim(varAcc1[varAcc1$greptype == "A_SNP",])[1])
varAcc1_transversion_GR <- GRanges(seqnames = varAcc1_transversion$ref_chr,
                                   ranges = IRanges(start = varAcc1_transversion$ref_start,
                                                    end = varAcc1_transversion$ref_end),
                                   strand = "*",
                                   coverage = rep(1, dim(varAcc1_transversion)[1]))

# Add transitions and transversions GRanges to varAcc1ListGR
varAcc1ListGR <- c(varAcc1ListGR,
                   varAcc1_transition_GR,
                   varAcc1_transversion_GR)
varAcc1classNames <- c(
                       "Acc1_all",
                       "Acc1_SNP",
                       "Acc1_insertion",
                       "Acc1_deletion",
                       "Acc1_indel",
                       "Acc1_transition",
                       "Acc1_transversion"
                      )


# Load table of varAcc2
varAcc2 <- fread(input = paste0(dirName, "/",
                                substr(x = refbase, start = 22, stop = 24), ".syri.out.gz"),
                  header = FALSE)
varAcc2 <- data.frame(varAcc2)
colnames(varAcc2) <- c("ref_chr", "ref_start", "ref_end", "ref_seq", "que_seq",
                       "que_chr", "que_start", "que_end", "uniqueID", "parentID",
                       "type", "copy_status")
varAcc2$ref_start <- as.numeric(varAcc2$ref_start)
varAcc2$ref_end <- as.numeric(varAcc2$ref_end)
varAcc2$que_start <- as.numeric(varAcc2$que_start)
varAcc2$que_end <- as.numeric(varAcc2$que_end)
varAcc2$ref_chr <- gsub(pattern = "Chr", replacement = "Col_Chr",
                        x = varAcc2$ref_chr)
varAcc2$que_chr <- gsub(pattern = "Chr", replacement = paste0(substr(x = refbase, start = 22, stop = 24), "_Chr"),
                        x = varAcc2$que_chr)

# Remove varAcc2 that are listed more than once relative to the reference genome (Columns 1-5)
varAcc2 <- varAcc2[-which(duplicated(varAcc2[,1:5])),]

# Extract rows corresponding to varAcc2 of type varType
if(varType == "SNP_INDEL") {
  varAcc2 <- varAcc2[which(varAcc2$type %in% c("SNP", "INS", "DEL")),]
} else if(!varType %in% c("all_variants")) {
  varAcc2 <- varAcc2[which(varAcc2$type == varType),]
}

# Flip variant coordinates in wrong orientation for GRanges conversion
varAcc2_end_start <- varAcc2[which(varAcc2$ref_end < varAcc2$ref_start),]
varAcc2 <- varAcc2[-which(varAcc2$ref_end < varAcc2$ref_start),]
varAcc2_start_end <- varAcc2_end_start
varAcc2_start_end$ref_start <- varAcc2_end_start$ref_end
varAcc2_start_end$ref_end <- varAcc2_end_start$ref_start
varAcc2 <- rbind(varAcc2,
                 varAcc2_start_end)

# Order varAcc2 according to chromosome and coordinates
varAcc2 <- varAcc2[ with(varAcc2, order(ref_chr, ref_start, ref_end)), ]

print(paste0("Number of varAcc2 of type ", varType, ":"))
print(nrow(varAcc2))

varAcc2 <- data.frame(varAcc2,
                      greptype = NA)

varAcc2[varAcc2$type == "SNP",]$greptype <- "A_SNP"
varAcc2[varAcc2$type == "INS",]$greptype <- "A_insertion"
varAcc2[varAcc2$type == "DEL",]$greptype <- "A_deletion"

# Convert into GRanges
varAcc2GR <- GRanges(seqnames = varAcc2$ref_chr,
                     ranges = IRanges(start = varAcc2$ref_start,
                                      end = varAcc2$ref_end),
                     strand = "*",
                     coverage = rep(1, dim(varAcc2)[1]))

# Subset varAcc2 by class (grep-ing for "A" within "greptype" field returns all varAcc2)

varAcc2class <- c(
                  "A",
                  "SNP",
                  "insertion",
                  "deletion",
                  "tion"
                 )
varAcc2ListGR <- mclapply(seq_along(varAcc2class), function(x) {
  classvarAcc2 <- varAcc2[grep(varAcc2class[x], varAcc2$greptype),]
  GRanges(seqnames = classvarAcc2$ref_chr,
          ranges = IRanges(start = classvarAcc2$ref_start,
                           end = classvarAcc2$ref_end),
          strand = "*",
          coverage = rep(1, dim(classvarAcc2)[1]))
}, mc.cores = length(varAcc2class))
stopifnot(length(varAcc2ListGR[[5]]) ==
          length(varAcc2ListGR[[4]]) + length(varAcc2ListGR[[3]]))

## transition
varAcc2_transition <- varAcc2[(varAcc2$ref_seq == "A" | varAcc2$ref_seq == "G") & (varAcc2$que_seq == "G" | varAcc2$que_seq == "A") |
                              (varAcc2$ref_seq == "C" | varAcc2$ref_seq == "T") & (varAcc2$que_seq == "T" | varAcc2$que_seq == "C"),]
varAcc2_transition_GR <- GRanges(seqnames = varAcc2_transition$ref_chr,
                                 ranges = IRanges(start = varAcc2_transition$ref_start,
                                                  end = varAcc2_transition$ref_end),
                                 strand = "*",
                                 coverage = rep(1, dim(varAcc2_transition)[1]))

## transversion
varAcc2_transversion <- varAcc2[(varAcc2$ref_seq == "A" | varAcc2$ref_seq == "G") & (varAcc2$que_seq == "C" | varAcc2$que_seq == "T") |
                                (varAcc2$ref_seq == "C" | varAcc2$ref_seq == "T") & (varAcc2$que_seq == "A" | varAcc2$que_seq == "G"),]
## Below sanity check doesn't work because due to use of IUPAC ambiguity codes for some alleles
## (e.g., found one in varAcc2[varAcc2$greptype == "A_SNP",][117920,] )
#stopifnot((dim(varAcc2_transition)[1] +
#           dim(varAcc2_transversion)[1]) ==
#           dim(varAcc2[varAcc2$greptype == "A_SNP",])[1])
varAcc2_transversion_GR <- GRanges(seqnames = varAcc2_transversion$ref_chr,
                                   ranges = IRanges(start = varAcc2_transversion$ref_start,
                                                    end = varAcc2_transversion$ref_end),
                                   strand = "*",
                                   coverage = rep(1, dim(varAcc2_transversion)[1]))

# Add transitions and transversions GRanges to varAcc2ListGR
varAcc2ListGR <- c(varAcc2ListGR,
                   varAcc2_transition_GR,
                   varAcc2_transversion_GR)
varAcc2classNames <- c(
                       "Acc2_all",
                       "Acc2_SNP",
                       "Acc2_insertion",
                       "Acc2_deletion",
                       "Acc2_indel",
                       "Acc2_transition",
                       "Acc2_transversion"
                      )

# Define matrix and column mean outfiles
varAcc1featAcc1outDF <- lapply(seq_along(varAcc1classNames), function(x) {
  list(paste0(matDir, varAcc1classNames[x],
              "_variant_freq_MappedOn_", refbase, "_Acc1_Chr_genes_in_", paste0(chrName, collapse = "_"),
              "_", genomeRegion, "_matrix_bin", binName, "_flank", flankName, ".tab"),
       paste0(matDir, varAcc1classNames[x],
              "_variant_freq_MappedOn_", refbase, "_Acc1_Chr_genes_in_", paste0(chrName, collapse = "_"),
              "_", genomeRegion, "_ranLocAcc1_matrix_bin", binName, "_flank", flankName, ".tab"))
})
varAcc1featAcc1outDFcolMeans <- lapply(seq_along(varAcc1classNames), function(x) {
  list(paste0(matDir, varAcc1classNames[x],
              "_variant_freq_MappedOn_", refbase, "_Acc1_Chr_genes_in_", paste0(chrName, collapse = "_"),
              "_", genomeRegion, "_matrix_bin", binName, "_flank", flankName, "_colMeans.tab"),
       paste0(matDir, varAcc1classNames[x],
              "_variant_freq_MappedOn_", refbase, "_Acc1_Chr_genes_in_", paste0(chrName, collapse = "_"),
              "_", genomeRegion, "_ranLocAcc1_matrix_bin", binName, "_flank", flankName, "_colMeans.tab"))
})

# Define matrix and column mean outfiles
varAcc1featAcc2outDF <- lapply(seq_along(varAcc1classNames), function(x) {
  list(paste0(matDir, varAcc1classNames[x],
              "_variant_freq_MappedOn_", refbase, "_Acc2_Chr_genes_in_", paste0(chrName, collapse = "_"),
              "_", genomeRegion, "_matrix_bin", binName, "_flank", flankName, ".tab"),
       paste0(matDir, varAcc1classNames[x],
              "_variant_freq_MappedOn_", refbase, "_Acc2_Chr_genes_in_", paste0(chrName, collapse = "_"),
              "_", genomeRegion, "_ranLocAcc2_matrix_bin", binName, "_flank", flankName, ".tab"))
})
varAcc1featAcc2outDFcolMeans <- lapply(seq_along(varAcc1classNames), function(x) {
  list(paste0(matDir, varAcc1classNames[x],
              "_variant_freq_MappedOn_", refbase, "_Acc2_Chr_genes_in_", paste0(chrName, collapse = "_"),
              "_", genomeRegion, "_matrix_bin", binName, "_flank", flankName, "_colMeans.tab"),
       paste0(matDir, varAcc1classNames[x],
              "_variant_freq_MappedOn_", refbase, "_Acc2_Chr_genes_in_", paste0(chrName, collapse = "_"),
              "_", genomeRegion, "_ranLocAcc2_matrix_bin", binName, "_flank", flankName, "_colMeans.tab"))
})

# Define matrix and column mean outfiles
varAcc2featAcc2outDF <- lapply(seq_along(varAcc2classNames), function(x) {
  list(paste0(matDir, varAcc2classNames[x],
              "_variant_freq_MappedOn_", refbase, "_Acc2_Chr_genes_in_", paste0(chrName, collapse = "_"),
              "_", genomeRegion, "_matrix_bin", binName, "_flank", flankName, ".tab"),
       paste0(matDir, varAcc2classNames[x],
              "_variant_freq_MappedOn_", refbase, "_Acc2_Chr_genes_in_", paste0(chrName, collapse = "_"),
              "_", genomeRegion, "_ranLocAcc2_matrix_bin", binName, "_flank", flankName, ".tab"))
})
varAcc2featAcc2outDFcolMeans <- lapply(seq_along(varAcc2classNames), function(x) {
  list(paste0(matDir, varAcc2classNames[x],
              "_variant_freq_MappedOn_", refbase, "_Acc2_Chr_genes_in_", paste0(chrName, collapse = "_"),
              "_", genomeRegion, "_matrix_bin", binName, "_flank", flankName, "_colMeans.tab"),
       paste0(matDir, varAcc2classNames[x],
              "_variant_freq_MappedOn_", refbase, "_Acc2_Chr_genes_in_", paste0(chrName, collapse = "_"),
              "_", genomeRegion, "_ranLocAcc2_matrix_bin", binName, "_flank", flankName, "_colMeans.tab"))
})

# Define matrix and column mean outfiles
varAcc2featAcc1outDF <- lapply(seq_along(varAcc2classNames), function(x) {
  list(paste0(matDir, varAcc2classNames[x],
              "_variant_freq_MappedOn_", refbase, "_Acc1_Chr_genes_in_", paste0(chrName, collapse = "_"),
              "_", genomeRegion, "_matrix_bin", binName, "_flank", flankName, ".tab"),
       paste0(matDir, varAcc2classNames[x],
              "_variant_freq_MappedOn_", refbase, "_Acc1_Chr_genes_in_", paste0(chrName, collapse = "_"),
              "_", genomeRegion, "_ranLocAcc1_matrix_bin", binName, "_flank", flankName, ".tab"))
})
varAcc2featAcc1outDFcolMeans <- lapply(seq_along(varAcc2classNames), function(x) {
  list(paste0(matDir, varAcc2classNames[x],
              "_variant_freq_MappedOn_", refbase, "_Acc1_Chr_genes_in_", paste0(chrName, collapse = "_"),
              "_", genomeRegion, "_matrix_bin", binName, "_flank", flankName, "_colMeans.tab"),
       paste0(matDir, varAcc2classNames[x],
              "_variant_freq_MappedOn_", refbase, "_Acc1_Chr_genes_in_", paste0(chrName, collapse = "_"),
              "_", genomeRegion, "_ranLocAcc1_matrix_bin", binName, "_flank", flankName, "_colMeans.tab"))
})

# Function to create variant frequency matrices for
# feature loci and random loci (incl. flanking regions)
# and to calculate mean profiles across all feature loci and random loci
covMatrix <- function(signal,
                      feature,
                      ranLoc,
                      featureSize,
                      flankSize,
                      binSize,
                      outDF,
                      outDFcolMeans) {
  # feature loci
  set.seed(2840)
  feature_smoothed <- normalizeToMatrix(signal = signal,
                                        target = feature,
                                        value_column = "coverage",
                                        extend = flankSize,
                                        mean_mode = "w0",
                                        w = binSize,
                                        background = 0,
                                        smooth = TRUE,
                                        include_target = TRUE,
                                        target_ratio = featureSize/(featureSize+(flankSize*2)))
  print("feature_smoothed")
  print(feature_smoothed)
  print("feature_smoothed rows = ")
  print(length(feature_smoothed)/round((featureSize/binSize)+((flankSize*2)/binSize)))
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
                                       w = binSize,
                                       background = 0,
                                       smooth = TRUE,
                                       include_target = TRUE,
                                       target_ratio = featureSize/(featureSize+(flankSize*2)))
  print("ranLoc_smoothed")
  print(ranLoc_smoothed)
  print("ranLoc_smoothed rows = ")
  print(length(ranLoc_smoothed)/round((featureSize/binSize)+((flankSize*2)/binSize)))
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
# varAcc1featAcc1
mclapply(seq_along(varAcc1classNames), function(x) {
  covMatrix(signal = varAcc1ListGR[[x]],
            feature = featuresAcc1_orthoGR,
            ranLoc = ranLocAcc1_orthoGR,
            featureSize = regionBodyLength,
            flankSize = upstream,
            binSize = binSize,
            outDF = varAcc1featAcc1outDF[[x]],
            outDFcolMeans = varAcc1featAcc1outDFcolMeans[[x]])
  print(paste0(varAcc1classNames[x],
               " Acc1 variant frequency around Acc1 orthologous genes in ", chrName,
               " profile calculation complete"))
}, mc.cores = length(varAcc1classNames))

# varAcc1featAcc2
mclapply(seq_along(varAcc1classNames), function(x) {
  covMatrix(signal = varAcc1ListGR[[x]],
            feature = featuresAcc2_orthoGR,
            ranLoc = ranLocAcc2_orthoGR,
            featureSize = regionBodyLength,
            flankSize = upstream,
            binSize = binSize,
            outDF = varAcc1featAcc2outDF[[x]],
            outDFcolMeans = varAcc1featAcc2outDFcolMeans[[x]])
  print(paste0(varAcc1classNames[x],
               " Acc1 variant frequency around Acc2 orthologous genes in ", chrName,
               " profile calculation complete"))
}, mc.cores = length(varAcc1classNames))

# Run covMatrix() function on each feature GRanges object to obtain matrices
# containing normalised feature density values around target and random loci
# varAcc2featAcc2
mclapply(seq_along(varAcc2classNames), function(x) {
  covMatrix(signal = varAcc2ListGR[[x]],
            feature = featuresAcc2_orthoGR,
            ranLoc = ranLocAcc2_orthoGR,
            featureSize = regionBodyLength,
            flankSize = upstream,
            binSize = binSize,
            outDF = varAcc2featAcc2outDF[[x]],
            outDFcolMeans = varAcc2featAcc2outDFcolMeans[[x]])
  print(paste0(varAcc2classNames[x],
               " Acc2 variant frequency around Acc2 orthologous genes in ", chrName,
               " profile calculation complete"))
}, mc.cores = length(varAcc2classNames))

# varAcc2featAcc1
mclapply(seq_along(varAcc2classNames), function(x) {
  covMatrix(signal = varAcc2ListGR[[x]],
            feature = featuresAcc1_orthoGR,
            ranLoc = ranLocAcc1_orthoGR,
            featureSize = regionBodyLength,
            flankSize = upstream,
            binSize = binSize,
            outDF = varAcc2featAcc1outDF[[x]],
            outDFcolMeans = varAcc2featAcc1outDFcolMeans[[x]])
  print(paste0(varAcc2classNames[x],
               " Acc2 variant frequency around Acc1 orthologous genes in ", chrName,
               " profile calculation complete"))
}, mc.cores = length(varAcc2classNames))
