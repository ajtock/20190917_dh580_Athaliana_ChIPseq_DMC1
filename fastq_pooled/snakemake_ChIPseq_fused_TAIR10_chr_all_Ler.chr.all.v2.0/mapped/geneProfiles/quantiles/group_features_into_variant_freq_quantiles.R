#!/applications/R/R-4.0.0/bin/Rscript

#
# Divide features into quantiles based on varType frequency per kb
# in a given feature region (e.g., promoters).
# Extract and save feature IDs for each quantile for further analyses
# (e.g., GO enrichment and average + 95% CI profile plotting).
# Plot feature quantile varType frequency per kb in a
# box-and-whisker plot or violin plot
#

# Usage:
# /applications/R/R-4.0.0/bin/Rscript group_features_into_variant_freq_quantiles.R SNP_INDEL '/data/public_data/arabidopsis/MPIPZ_Jiao_Schneeberger_2020_NatCommun/Ler' both 'Chr1,Chr2,Chr3,Chr4,Chr5' 2200 2000 2kb '2 kb' 10 10bp genes genomewide 6 fused_TAIR10_chr_all_Ler.chr.all.v2.0

#varType <- "SNP_INDEL"
#dirName <- "/data/public_data/arabidopsis/MPIPZ_Jiao_Schneeberger_2020_NatCommun/Ler"
#align <- "both"
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
#refbase <- "fused_TAIR10_chr_all_Ler.chr.all.v2.0"

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

options(stringsAsFactors = F)
library(EnrichedHeatmap)
library(png)
#library(Cairo)
library(RColorBrewer)
library(circlize)
library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(data.table)
library(parallel)
library(doParallel)
registerDoParallel(cores = detectCores())
print("Currently registered parallel backend name, version and cores")
print(getDoParName())
print(getDoParVersion())
print(getDoParWorkers())

outDir <- paste0(paste0(chrName, collapse = "_"),
                 "/quantiles_", genomeRegion, "_by_", varType,
                 "_freq_in_", orderRegion, "/")
plotDir <- paste0(outDir, "plots/")
system(paste0("[ -d ", outDir, " ] || mkdir -p ", outDir))
system(paste0("[ -d ", plotDir, " ] || mkdir -p ", plotDir))

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

## Load table of representative feature coordinates in BED format
## to be used for extracting ortholog rows containing matching feature parent IDs
#features <- lapply(seq_along(chrName), function(x) {
#  read.table(paste0(dirName, "/", refbase, "/genes/",
#                    "Araport11_representative_mRNA_Col_",
#                    chrName[x], ".bed"),
#             header = F)
#})
#if(length(chrName) > 1) {
#  features <- do.call(rbind, features)
#} else {
#  features <- features[[1]]
#}
## Convert 0-based start coordinates (BED)
## into 1-based start coordinates (for output as TSV below)
#features[,2] <- features[,2]+1
#colnames(features) <- c("chr", "start", "end", "featureID", "score", "strand")
#featuresGR <- GRanges(seqnames = features$chr,
#                      ranges = IRanges(start = features$start,
#                                       end = features$end),
#                      strand = features$strand,
#                      featureID = features$featureID)
#
## Load table of ranLoc coordinates in BED format
#ranLocs <- lapply(seq_along(chrName), function(x) {
#  read.table(paste0(dirName, "/", refbase, "/genes/",
#                    "Araport11_representative_mRNA_",
#                    chrName[x], "_randomLoci.bed"),
#             header = F)
#})
#if(length(chrName) > 1) {
#  ranLocs <- do.call(rbind, ranLocs)
#} else {
#  ranLocs <- ranLocs[[1]]
#}
## Convert 0-based start coordinates (BED)
## into 1-based start coordinates (for output as TSV below)
#ranLocs[,2] <- ranLocs[,2]+1
#colnames(ranLocs) <- c("chr", "start", "end", "ranLocID", "score", "strand")
#ranLocsGR <- GRanges(seqnames = ranLocs$chr,
#                     ranges = IRanges(start = ranLocs$start+1,
#                                      end = ranLocs$end),
#                     strand = ranLocs$strand,
#                     ranLocID = ranLocs$ranLocID)

# Obtain orthologs in the first accession (Acc1) represented in the fused genome (Col)
featuresAcc1 <- readGFF(paste0(dirName, "/",
                               refbase, "/Araport11_GFF3_genes_transposons.201606.gff"))
#featuresAcc1$seqid <- gsub("Chr", "Col_Chr", featuresAcc1$seqid)
featuresAcc1 <- featuresAcc1[featuresAcc1$seqid %in% chrName,]
featuresAcc1 <- featuresAcc1[featuresAcc1$type == "mRNA",]
print(dim(featuresAcc1))
#[1] 52060    25

# Obtain orthologs in the second accession (Acc2) represented in the fused genome (Cvi or Ler)
featuresAcc2 <- readGFF(paste0(dirName, "/",
                               substr(x = refbase, start = 22, stop = 24),
                               ".protein-coding.genes.v2.5.2019-10-09.gff3.gz"))
featuresAcc2$seqid <- gsub("chr", "Chr", featuresAcc2$seqid)
featuresAcc2 <- featuresAcc2[featuresAcc2$seqid %in% chrName,]
featuresAcc2 <- featuresAcc2[featuresAcc2$type == "mRNA",]
print(dim(featuresAcc2))
#[1] 26911    11 

# Load table of orthologs
orthologs <- fread(paste0(substr(x = dirName, start = 1, stop = 68), "/",
                           "orthofinder/ajt200_run_v200721/protein_fasta/OrthoFinder/Results_Jul29/Orthologues/Orthologues_Col/",
                           "Col__v__", substr(x = refbase, start = 22, stop = 24), ".tsv"),
                   header = T)
orthologs <- data.frame(orthologs)
#orthologs <- orthologs[,which(colnames(orthologs) %in% c("Orthogroup", "Col", substr(x = refbase, start = 22, stop = 24)))]

# Use strsplit to convert comma-separated feature ID character strings (one string per accession per ortholog)
# into a list in which each list element is a vector of feature IDs
ortholog_featureIDs_Acc1 <- strsplit(orthologs[,2], split = ", ")
ortholog_featureIDs_Acc2 <- strsplit(orthologs[,3], split = ", ")

## If using representative "features" object instead of "featuresAcc1" object containing all mRNA loci
## Count number of Accesion 1 (Acc1) feature IDs per ortholog that are in the Araport11 feature parent IDs 
#ortholog_featureIDs_Acc1_in_features_sum <- unlist(mclapply(seq_along(ortholog_featureIDs_Acc1), function(x) {
#  sum( gsub(pattern = "\\.\\d+", replacement = "", x = ortholog_featureIDs_Acc1[[x]]) %in%
#       gsub(pattern = "\\.\\d+", replacement = "", x = features$featureID) )
#}, mc.preschedule = TRUE, mc.cores = detectCores()))

# Count number of Accesion 1 (Acc1) feature IDs per ortholog that are in the Araport11 feature IDs 
ortholog_featureIDs_Acc1_in_featuresAcc1_sum <- unlist(mclapply(seq_along(ortholog_featureIDs_Acc1), function(x) {
  sum( ortholog_featureIDs_Acc1[[x]] %in%
       featuresAcc1$ID )
}, mc.preschedule = TRUE, mc.cores = detectCores()))

# Count number of Accession 1 (Acc1) feature IDs per ortholog
ortholog_featureIDs_Acc1_len <- unlist(mclapply(seq_along(ortholog_featureIDs_Acc1), function(x) {
  length( ortholog_featureIDs_Acc1[[x]] )
}, mc.preschedule = TRUE, mc.cores = detectCores()))

# Count number of Accession 2 (Acc2) feature IDs per ortholog
ortholog_featureIDs_Acc2_len <- unlist(mclapply(seq_along(ortholog_featureIDs_Acc2), function(x) {
  length( ortholog_featureIDs_Acc2[[x]] )
}, mc.preschedule = TRUE, mc.cores = detectCores()))

# Get indices for orthologs that contain the
# same number of feature IDs in both accessions, and
# retain those where the number of Accession 1 (Acc1) (and Accession 2; Acc2) features IDs per ortholog is 1, and
# retain those where the number of Accession 1 (Acc1) feature IDs
# per ortholog that are in the Araport11 feature parent IDs is > 0
ortholog_featureIDs_Acc1_in_features_indices <- which(ortholog_featureIDs_Acc1_len == ortholog_featureIDs_Acc2_len &
                                                      ortholog_featureIDs_Acc1_len == 1 &
                                                      ortholog_featureIDs_Acc1_in_featuresAcc1_sum > 0)

orthologs <- orthologs[ortholog_featureIDs_Acc1_in_features_indices,]

# Retain rows in orthologs where both the Acc1 feature ID and the Acc2 feature ID
# is each present in the ID column of featuresAcc1 or featuresAcc2
orthologs <- orthologs[which(orthologs[,2] %in% featuresAcc1$ID &
                             orthologs[,3] %in% featuresAcc2$ID),]

# Extend featuresAcc1 by appending columns from "orthologs" data.frame,
# and order featuresAcc1_ortho by Acc1 feature IDs
featuresAcc1_ortho <- merge(x = featuresAcc1,
                            y = orthologs,
                            by.x = "ID",
                            by.y = colnames(orthologs)[2])
featuresAcc1_ortho <- data.frame(featuresAcc1_ortho[,-dim(featuresAcc1_ortho)[2]],
                                 Col = featuresAcc1_ortho[,1],
                                 placeholder = featuresAcc1_ortho[,dim(featuresAcc1_ortho)[2]])
names(featuresAcc1_ortho)[names(featuresAcc1_ortho) == "placeholder"] <- substr(x = refbase, start = 22, stop = 24)
featuresAcc1_ortho <- featuresAcc1_ortho[ with(featuresAcc1_ortho, order(Col)), ]

# Extend featuresAcc2 by appending columns from "orthologs" data.frame,
# and order featuresAcc2_ortho by Acc1 feature IDs
# Note that fewer
# Acc2 ortholog feature IDs are in featuresAcc2$ID than
# Acc1 ortholog feature IDs are in featuresAcc1$ID
featuresAcc2_ortho <- merge(x = featuresAcc2,
                            y = orthologs,
                            by.x = "ID",
                            by.y = colnames(orthologs)[3])
featuresAcc2_ortho <- data.frame(featuresAcc2_ortho,
                                 placeholder = featuresAcc2_ortho[,1])
names(featuresAcc2_ortho)[names(featuresAcc2_ortho) == "placeholder"] <- substr(x = refbase, start = 22, stop = 24)
featuresAcc2_ortho <- featuresAcc2_ortho[ with(featuresAcc2_ortho, order(Col)), ]

# Ensure row order (by increasing "Col" feature ID column)
# is consistent for featuresAcc1_ortho and featuresAcc2_ortho 
stopifnot(identical(featuresAcc1_ortho$Col, featuresAcc2_ortho$Col))
stopifnot(identical(featuresAcc1_ortho$Col, featuresAcc1_ortho$ID))
stopifnot(identical(featuresAcc2_ortho$Ler, featuresAcc2_ortho$ID))

# Load table of variants
variants <- fread(input = paste0(dirName, "/",
                                 substr(x = refbase, start = 22, stop = 24), ".syri.out.gz"),
                  header = FALSE)
variants <- data.frame(variants)
colnames(variants) <- c("ref_chr", "ref_start", "ref_end", "ref_seq", "que_seq",
                        "que_chr", "que_start", "que_end", "uniqueID", "parentID",
                        "type", "copy_status")
variants$ref_start <- as.numeric(variants$ref_start)
variants$ref_end <- as.numeric(variants$ref_end)
variants$que_start <- as.numeric(variants$que_start)
variants$que_end <- as.numeric(variants$que_end)
variants$ref_chr <- gsub(pattern = "Chr", replacement = "Col_Chr",
                         x = variants$ref_chr)
variants$que_chr <- gsub(pattern = "Chr", replacement = paste0(substr(x = refbase, start = 22, stop = 24), "_Chr"),
                         x = variants$que_chr)

# Remove variants that are listed more than once relative to the reference genome (Columns 1-5)
variants <- variants[-which(duplicated(variants[,1:5])),]

# Extract rows corresponding to variants of type varType
if(varType == "SNP_INDEL") {
  variants <- variants[which(variants$type %in% c("SNP", "INS", "DEL")),]
} else if(!varType %in% c("all_variants")) {
  variants <- variants[which(variants$type == varType),]
}

# Flip variant coordinates in wrong orientation for GRanges conversion
variants_end_start <- variants[which(variants$ref_end < variants$ref_start),]
variants <- variants[-which(variants$ref_end < variants$ref_start),]
variants_start_end <- variants_end_start
variants_start_end$ref_start <- variants_end_start$ref_end
variants_start_end$ref_end <- variants_end_start$ref_start
variants <- rbind(variants,
                  variants_start_end)

# Order variants according to chromosome and coordinates
variants <- variants[ with(variants, order(ref_chr, ref_start, ref_end)), ]

print(paste0("Number of variants of type ", varType, ":"))
print(nrow(variants))

# Convert into GRanges
variantsAcc1GR <- GRanges(seqnames = variants$ref_chr,
                          ranges = IRanges(start = variants$ref_start,
                                           end = variants$ref_end),
                          strand = "*",
                          ref_seq = variants$ref_seq,
                          que_seq = variants$que_seq,
                          uniqueID = variants$uniqueID,
                          type = variants$type)

# Load table of variants
variants <- fread(input = paste0(dirName, "/",
                                 substr(x = refbase, start = 22, stop = 24), ".syri.out.gz"),
                  header = FALSE)
variants <- data.frame(variants)
colnames(variants) <- c("ref_chr", "ref_start", "ref_end", "ref_seq", "que_seq",
                        "que_chr", "que_start", "que_end", "uniqueID", "parentID",
                        "type", "copy_status")
variants$ref_start <- as.numeric(variants$ref_start)
variants$ref_end <- as.numeric(variants$ref_end)
variants$que_start <- as.numeric(variants$que_start)
variants$que_end <- as.numeric(variants$que_end)
variants$ref_chr <- gsub(pattern = "Chr", replacement = "Col_Chr",
                         x = variants$ref_chr)
variants$que_chr <- gsub(pattern = "Chr", replacement = paste0(substr(x = refbase, start = 22, stop = 24), "_Chr"),
                         x = variants$que_chr)

# Remove variants that are listed more than once relative to the reference genome (Columns 4-8)
variants <- variants[-which(duplicated(variants[,4:8])),]

# Extract rows corresponding to variants of type varType
if(varType == "SNP_INDEL") {
  variants <- variants[which(variants$type %in% c("SNP", "INS", "DEL")),]
} else if(!varType %in% c("all_variants")) {
  variants <- variants[which(variants$type == varType),]
}

# Flip variant coordinates in wrong orientation for GRanges conversion
variants_end_start <- variants[which(variants$que_end < variants$que_start),]
variants <- variants[-which(variants$que_end < variants$que_start),]
variants_start_end <- variants_end_start
variants_start_end$que_start <- variants_end_start$que_end
variants_start_end$que_end <- variants_end_start$que_start
variants <- rbind(variants,
                  variants_start_end)

# Order variants according to chromosome and coordinates
variants <- variants[ with(variants, order(que_chr, que_start, que_end)), ]

print(paste0("Number of variants of type ", varType, ":"))
print(nrow(variants))

# Convert into GRanges
variantsAcc2GR <- GRanges(seqnames = variants$que_chr,
                          ranges = IRanges(start = variants$que_start,
                                           end = variants$que_end),
                          strand = "*",
                          ref_seq = variants$ref_seq,
                          que_seq = variants$que_seq,
                          uniqueID = variants$uniqueID,
                          type = variants$type)



# Convert featuresAcc1_ortho into GRanges
featuresAcc1_orthoGR <- GRanges(seqnames = paste0("Col_", featuresAcc1_ortho$seqid),
                                ranges = IRanges(start = featuresAcc1_ortho$start,
                                                 end = featuresAcc1_ortho$end),
                                strand = featuresAcc1_ortho$strand,
                                featureID = featuresAcc1_ortho$ID)
featuresAcc1_orthoGR_ol <- findOverlaps(query = featuresAcc1_orthoGR,
                                        subject = genomeMaskGR,
                                        type = "any", select = "all",
                                        ignore.strand = TRUE)
featuresAcc1_orthoGR <- featuresAcc1_orthoGR[-unique(queryHits(featuresAcc1_orthoGR_ol))]
featuresAcc1_ortho <- featuresAcc1_ortho[which(featuresAcc1_ortho$ID %in% featuresAcc1_orthoGR$featureID),]
stopifnot(all.equal(featuresAcc1_ortho$ID, featuresAcc1_orthoGR$featureID))

# Get ranges corresponding to orderRegion
if(orderRegion == "bodies") {
  featOrReAcc1_orthoGR <- featuresAcc1_orthoGR
} else if(orderRegion == "promoters") {
  # Obtain 1000 bp upstream of start coordinates
  featOrReAcc1_orthoGR <- promoters(featuresAcc1_orthoGR, upstream = 1000, downstream = 0)
} else if(orderRegion == "terminators") {
  # Obtain 1000 bp downstream of end coordinates
  source("/projects/ajt200/Rfunctions/TTSplus.R")
  featOrReAcc1_orthoGR <- TTSplus(featuresAcc1_orthoGR, upstream = -1, downstream = 1000)
} else if(orderRegion == "genes") {
  featOrReAcc1_orthoGR <- GRanges(seqnames = paste0("Col_", featuresAcc1_ortho$seqid),
                                  ranges = IRanges(start = featuresAcc1_ortho$start-1000,
                                                   end = featuresAcc1_ortho$end+1000),
                                  strand = featuresAcc1_ortho$strand,
                                  featureID = featuresAcc1_ortho$ID)
} else {
  stop("orderRegion is none of bodies, promoters, terminators or genes")
}


# Define random loci of the same number and width distribution,
# and in the same per-feature genomeRegionGR as featuresAcc1_orthoGR

# Define function to select randomly positioned loci of the same
# width distribution as featuresAcc1_orthoGR
ranLocStartSelect <- function(coordinates, n) {
  sample(x = coordinates,
         size = n,
         replace = FALSE)
}

# Disable scientific notation (e.g., 59000000 rather than 5.9e+07)
options(scipen = 100)

# Apply ranLocStartSelect() on a per-chromosome basis so that
# ranLocAcc1GR contains the same number of loci per chromosome as featuresAcc1_orthoGR
chrs <- seqlevels(sortSeqlevels(featuresAcc1_orthoGR))
ranLocAcc1GR <- GRanges()
for(i in 1:length(chrs)) {
  featuresAcc1_orthoChrGR <- featuresAcc1_orthoGR[seqnames(featuresAcc1_orthoGR) == chrs[i]]
  genomeRegionChrGR <- genomeRegionGR[seqnames(genomeRegionGR) == chrs[i]]
  # Contract genomeRegionChrGR so that random loci and 2-kb flanking regions
  # do not extend beyond chromosome ends
  end(genomeRegionChrGR) <- end(genomeRegionChrGR)-max(width(featuresAcc1_orthoChrGR))-2000
  start(genomeRegionChrGR) <- start(genomeRegionChrGR)+2000
  # Define seed so that random selections are reproducible
  set.seed(93750174)
  ranLocAcc1ChrStart <- ranLocStartSelect(coordinates = unlist(lapply(seq_along(genomeRegionChrGR), function(x) {
                                                                     start(genomeRegionChrGR[x]) : end(genomeRegionChrGR[x])
                                                                   })),
                                              n = length(featuresAcc1_orthoChrGR))
  ranLocAcc1ChrGR <- GRanges(seqnames = chrs[i],
                             ranges = IRanges(start = ranLocAcc1ChrStart,
                                              width = width(featuresAcc1_orthoChrGR)),
                             strand = strand(featuresAcc1_orthoChrGR),
                             Col_featureID = featuresAcc1_orthoChrGR$featureID)
  ranLocAcc1GR <- append(ranLocAcc1GR, ranLocAcc1ChrGR)
}
ranLocAcc1GR <- sort(ranLocAcc1GR, by = ~ Col_featureID)
stopifnot( identical( as.character(seqnames(featuresAcc1_orthoGR)),
                      as.character(seqnames(ranLocAcc1GR)) ) )
stopifnot( length( findOverlaps(query = ranLocAcc1GR,
                                subject = genomeMaskGR,
                                type = "any", select = "all",
                                ignore.strand = TRUE) ) == 0 )

# Get ranges corresponding to orderRegion
if(orderRegion == "bodies") {
  ranLORAcc1GR <- ranLocAcc1GR
} else if(orderRegion == "promoters") {
  # Obtain 1000 bp upstream of start coordinates
  ranLORAcc1GR <- promoters(ranLocAcc1GR, upstream = 1000, downstream = 0)
} else if(orderRegion == "terminators") {
  # Obtain 1000 bp downstream of end coordinates
  source("/projects/ajt200/Rfunctions/TTSplus.R")
  ranLORAcc1GR <- TTSplus(ranLocAcc1GR, upstream = -1, downstream = 1000)
} else if(orderRegion == "genes") {
  ranLORAcc1GR <- GRanges(seqnames = as.character(seqnames(ranLocAcc1GR)),
                          ranges = IRanges(start = start(ranLocAcc1GR)-1000,
                                           end = end(ranLocAcc1GR)+1000),
                          strand = strand(ranLocAcc1GR))
} else {
  stop("orderRegion is none of bodies, promoters, terminators or genes")
}

# Count variants within orderRegion of each featuresAcc1_ortho and
# divide by orderRegion width in kb
featOrReAcc1_orthoGR_variants <- countOverlaps(query = featOrReAcc1_orthoGR,
                                               subject = variantsAcc1GR,
                                               type = "any",
                                               ignore.strand = TRUE)
featOrReAcc1_orthoGR_variantsPK <- featOrReAcc1_orthoGR_variants/(width(featOrReAcc1_orthoGR)/1000)

# Add variants per kb values to featuresAcc1_ortho dataframe
featuresAcc1_ortho <- data.frame(featuresAcc1_ortho,
                                 orderRegion_variantsPK = featOrReAcc1_orthoGR_variantsPK)

featuresAcc1_ortho$seqid <- paste0("Col_", featuresAcc1_ortho$seqid)

# Group featuresAcc1_ortho into quantiles according to decreasing orderRegion_variantsPK
featuresAcc1_ortho_DF <- data.frame(featuresAcc1_ortho,
                                    percentile = rank(featuresAcc1_ortho[,which(colnames(featuresAcc1_ortho) == "orderRegion_variantsPK")]) /
                                                 length(featuresAcc1_ortho[,which(colnames(featuresAcc1_ortho) == "orderRegion_variantsPK")]),
                                    quantile = as.character(""))
for(x in as.vector(which(sapply(featuresAcc1_ortho_DF, is.list)))) {
  featuresAcc1_ortho_DF[,x] <- sapply(featuresAcc1_ortho_DF[,x], paste, collapse = "; ")
}

# Count variants within orderRegion of each ranLocAcc1GR and
# divide by orderRegion width in kb
ranLORAcc1GR_variants <- countOverlaps(query = ranLORAcc1GR,
                                       subject = variantsAcc1GR,
                                       type = "any",
                                       ignore.strand = TRUE)
ranLORAcc1GR_variantsPK <- ranLORAcc1GR_variants/(width(ranLORAcc1GR)/1000)

# Add variants per kb values to ranLocAcc1 dataframe
ranLocAcc1 <- data.frame(ranLocAcc1GR,
                         orderRegion_variantsPK = ranLORAcc1GR_variantsPK)
ranLocAcc1 <- ranLocAcc1[,-which(names(ranLocAcc1) %in% c("width", "Col_featureID"))]
names(ranLocAcc1)[names(ranLocAcc1) == "seqnames"] <- "seqid"

# Group ranLocAcc1 into quantiles according to decreasing orderRegion_variantsPK
ranLocAcc1_DF <- data.frame(ranLocAcc1,
                            percentile = rank(ranLocAcc1[,which(colnames(ranLocAcc1) == "orderRegion_variantsPK")]) /
                                         length(ranLocAcc1[,which(colnames(ranLocAcc1) == "orderRegion_variantsPK")]),
                            random = as.character(""))

quantilesStats <- data.frame()
for(k in 1:quantiles) {
  # First quantile should span 1 to greater than, e.g., 0.75 proportions of featuresAcc1_ortho
  if(k < quantiles) {
    featuresAcc1_ortho_DF[ !is.na(featuresAcc1_ortho_DF[,which(colnames(featuresAcc1_ortho_DF) == "orderRegion_variantsPK")]) &
                           rank(featuresAcc1_ortho_DF[,which(colnames(featuresAcc1_ortho_DF) == "orderRegion_variantsPK")]) /
                           length(featuresAcc1_ortho_DF[,which(colnames(featuresAcc1_ortho_DF) == "orderRegion_variantsPK")]) <=
                           1-((k-1)/quantiles) &
                           rank(featuresAcc1_ortho_DF[,which(colnames(featuresAcc1_ortho_DF) == "orderRegion_variantsPK")]) /
                           length(featuresAcc1_ortho_DF[,which(colnames(featuresAcc1_ortho_DF) == "orderRegion_variantsPK")]) >
                           1-(k/quantiles), ]$quantile <- paste0("Quantile ", k)
  } else {
  # Final quantile should span 0 to, e.g., 0.25 proportions of featuresAcc1_ortho
    featuresAcc1_ortho_DF[ !is.na(featuresAcc1_ortho_DF[,which(colnames(featuresAcc1_ortho_DF) == "orderRegion_variantsPK")]) &
                           rank(featuresAcc1_ortho_DF[,which(colnames(featuresAcc1_ortho_DF) == "orderRegion_variantsPK")]) /
                           length(featuresAcc1_ortho_DF[,which(colnames(featuresAcc1_ortho_DF) == "orderRegion_variantsPK")]) <=
                           1-((k-1)/quantiles) &
                           rank(featuresAcc1_ortho_DF[,which(colnames(featuresAcc1_ortho_DF) == "orderRegion_variantsPK")]) /
                           length(featuresAcc1_ortho_DF[,which(colnames(featuresAcc1_ortho_DF) == "orderRegion_variantsPK")]) >=
                           1-(k/quantiles), ]$quantile <- paste0("Quantile ", k)
  }
  write.table(featuresAcc1_ortho_DF[featuresAcc1_ortho_DF$quantile == paste0("Quantile ", k),],
              file = paste0(outDir,
                            "quantile", k, "_of_", quantiles,
                            "_", genomeRegion, "_by_", varType, "_freq_in_", orderRegion,
                            "_of_Acc1_Chr_genes_in_", refbase, "_",
                            paste0(chrName, collapse = "_"), ".tsv"),
              quote = FALSE, sep = "\t", row.names = FALSE)
  stats <- data.frame(quantile = as.integer(k),
                      n = as.integer(dim(featuresAcc1_ortho_DF[featuresAcc1_ortho_DF$quantile == paste0("Quantile ", k),])[1]),
                      mean_width = as.integer(round(mean(
                        (featuresAcc1_ortho_DF[featuresAcc1_ortho_DF$quantile == paste0("Quantile ", k),]$end -
                         featuresAcc1_ortho_DF[featuresAcc1_ortho_DF$quantile == paste0("Quantile ", k),]$start) + 1, na.rm = T))),
                      total_width = as.integer(sum(
                        (featuresAcc1_ortho_DF[featuresAcc1_ortho_DF$quantile == paste0("Quantile ", k),]$end -
                         featuresAcc1_ortho_DF[featuresAcc1_ortho_DF$quantile == paste0("Quantile ", k),]$start) + 1, na.rm = T)),
                      mean_orderRegion_variantsPK = as.numeric(mean(featuresAcc1_ortho_DF[featuresAcc1_ortho_DF$quantile == paste0("Quantile ", k),][,which(colnames(featuresAcc1_ortho_DF) == "orderRegion_variantsPK")], na.rm = T)))
  quantilesStats <- rbind(quantilesStats, stats)
}
write.table(quantilesStats,
            file = paste0(outDir,
                          "summary_", quantiles, "quantiles",
                          "_", genomeRegion, "_by_", varType, "_freq_in_", orderRegion,
                          "_of_Acc1_Chr_genes_in_", refbase, "_",
                          paste0(chrName, collapse = "_"), ".tsv"),
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(featuresAcc1_ortho_DF,
            file = paste0(outDir,
                          "featuresAcc1_ortho_", quantiles, "quantiles",
                          "_", genomeRegion, "_by_", varType, "_freq_in_", orderRegion,
                          "_of_Acc1_Chr_genes_in_", refbase, "_",
                          paste0(chrName, collapse = "_"), ".tsv"),
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(featuresAcc1_ortho_DF[,2:10],
            file = paste0(outDir,
                          "featuresAcc1_ortho_", quantiles, "quantiles",
                          "_", genomeRegion, "_by_", varType, "_freq_in_", orderRegion,
                          "_of_Acc1_Chr_genes_in_", refbase, "_",
                          paste0(chrName, collapse = "_"), ".gff3"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
featuresAcc1_ortho_DF_bed <- data.frame(chr = as.character(featuresAcc1_ortho_DF[,2]),
                                        start = as.integer(featuresAcc1_ortho_DF[,5]-1),
                                        end = as.integer(featuresAcc1_ortho_DF[,6]),
                                        name = as.character(featuresAcc1_ortho_DF[,1]),
                                        score = as.numeric(featuresAcc1_ortho_DF[,29]),
                                        strand = as.character(featuresAcc1_ortho_DF[,8]))
write.table(featuresAcc1_ortho_DF_bed,
            file = paste0(outDir,
                          "featuresAcc1_ortho_", quantiles, "quantiles",
                          "_", genomeRegion, "_by_", varType, "_freq_in_", orderRegion,
                          "_of_Acc1_Chr_genes_in_", refbase, "_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


# Divide ranLocAcc1_DF into quantiles based on featuresAcc1_ortho_DF quantile indices
# Get row indices for each feature quantile
quantileIndices <- lapply(1:quantiles, function(k) {
  which(featuresAcc1_ortho_DF$quantile == paste0("Quantile ", k))
})
for(k in 1:quantiles) {
  ranLocAcc1_DF[quantileIndices[[k]],]$random <- paste0("Random ", k)
}
write.table(ranLocAcc1_DF,
            file = paste0(outDir,
                          "featuresAcc1_ortho_", quantiles, "quantiles",
                          "_", genomeRegion, "_by_", varType, "_freq_in_", orderRegion,
                          "_of_Acc1_Chr_genes_in_", refbase, "_",
                          paste0(chrName, collapse = "_"), "_ranLocAcc1.tsv"),
            quote = FALSE, sep = "\t", row.names = FALSE)
ranLocAcc1_DF_bed <- data.frame(chr = as.character(ranLocAcc1_DF[,1]),
                                start = as.integer(ranLocAcc1_DF[,2]-1),
                                end = as.integer(ranLocAcc1_DF[,3]),
                                name = as.integer(1:nrow(ranLocAcc1_DF)),
                                score = as.numeric(ranLocAcc1_DF[,5]),
                                strand = as.character(ranLocAcc1_DF[,4]))
write.table(ranLocAcc1_DF_bed,
            file = paste0(outDir,
                          "featuresAcc1_ortho_", quantiles, "quantiles",
                          "_", genomeRegion, "_by_", varType, "_freq_in_", orderRegion,
                          "_of_Acc1_Chr_genes_in_", refbase, "_",
                          paste0(chrName, collapse = "_"), "_ranLocAcc1.bed"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)



# Convert featuresAcc2_ortho into GRanges
featuresAcc2_ortho <- featuresAcc2_ortho[which(featuresAcc2_ortho$Col %in% featuresAcc1_ortho$Col),]
featuresAcc2_orthoGR <- GRanges(seqnames = paste0(substr(x = refbase, start = 22, stop = 24), "_", featuresAcc2_ortho$seqid),
                                ranges = IRanges(start = featuresAcc2_ortho$start,
                                                 end = featuresAcc2_ortho$end),
                                strand = featuresAcc2_ortho$strand,
                                featureID = featuresAcc2_ortho$ID,
                                Col_featureID = featuresAcc2_ortho$Col)

# Get ranges corresponding to orderRegion
if(orderRegion == "bodies") {
  featOrReAcc2_orthoGR <- featuresAcc2_orthoGR
} else if(orderRegion == "promoters") {
  # Obtain 1000 bp upstream of start coordinates
  featOrReAcc2_orthoGR <- promoters(featuresAcc2_orthoGR, upstream = 1000, downstream = 0)
} else if(orderRegion == "terminators") {
  # Obtain 1000 bp downstream of end coordinates
  source("/projects/ajt200/Rfunctions/TTSplus.R")
  featOrReAcc2_orthoGR <- TTSplus(featuresAcc2_orthoGR, upstream = -1, downstream = 1000)
} else if(orderRegion == "genes") {
  featOrReAcc2_orthoGR <- GRanges(seqnames = paste0(substr(x = refbase, start = 22, stop = 24), "_", featuresAcc2_ortho$seqid),
                                  ranges = IRanges(start = featuresAcc2_ortho$start-1000,
                                                   end = featuresAcc2_ortho$end+1000),
                                  strand = featuresAcc2_ortho$strand,
                                  featureID = featuresAcc2_ortho$ID)
} else {
  stop("orderRegion is none of bodies, promoters, terminators or genes")
}


# Define random loci of the same number and width distribution,
# and in the same per-feature genomeRegionGR as featuresAcc2_orthoGR

# Define function to select randomly positioned loci of the same
# width distribution as featuresAcc2_orthoGR
ranLocStartSelect <- function(coordinates, n) {
  sample(x = coordinates,
         size = n,
         replace = FALSE)
}

# Disable scientific notation (e.g., 59000000 rather than 5.9e+07)
options(scipen = 100)

# Apply ranLocStartSelect() on a per-chromosome basis so that
# ranLocAcc2GR contains the same number of loci per chromosome as featuresAcc2_orthoGR
chrs <- seqlevels(sortSeqlevels(featuresAcc2_orthoGR))
ranLocAcc2GR <- GRanges()
for(i in 1:length(chrs)) {
  featuresAcc2_orthoChrGR <- featuresAcc2_orthoGR[seqnames(featuresAcc2_orthoGR) == chrs[i]]
  genomeRegionChrGR <- genomeRegionGR[seqnames(genomeRegionGR) == chrs[i]]
  # Contract genomeRegionChrGR so that random loci and 2-kb flanking regions
  # do not extend beyond chromosome ends
  end(genomeRegionChrGR) <- end(genomeRegionChrGR)-max(width(featuresAcc2_orthoChrGR))-2000
  start(genomeRegionChrGR) <- start(genomeRegionChrGR)+2000
  # Define seed so that random selections are reproducible
  set.seed(93750174)
  ranLocAcc2ChrStart <- ranLocStartSelect(coordinates = unlist(lapply(seq_along(genomeRegionChrGR), function(x) {
                                                                     start(genomeRegionChrGR[x]) : end(genomeRegionChrGR[x])
                                                                   })),
                                              n = length(featuresAcc2_orthoChrGR))
  ranLocAcc2ChrGR <- GRanges(seqnames = chrs[i],
                             ranges = IRanges(start = ranLocAcc2ChrStart,
                                              width = width(featuresAcc2_orthoChrGR)),
                             strand = strand(featuresAcc2_orthoChrGR),
                             Col_featureID = featuresAcc2_orthoChrGR$Col_featureID)
  ranLocAcc2GR <- append(ranLocAcc2GR, ranLocAcc2ChrGR)
}
ranLocAcc2GR <- sort(ranLocAcc2GR, by = ~ Col_featureID)
stopifnot( identical( as.character(seqnames(featuresAcc2_orthoGR)),
                      as.character(seqnames(ranLocAcc2GR)) ) )
stopifnot( length( findOverlaps(query = ranLocAcc2GR,
                                subject = genomeMaskGR,
                                type = "any", select = "all",
                                ignore.strand = TRUE) ) == 0 )

# Get ranges corresponding to orderRegion
if(orderRegion == "bodies") {
  ranLORAcc2GR <- ranLocAcc2GR
} else if(orderRegion == "promoters") {
  # Obtain 1000 bp upstream of start coordinates
  ranLORAcc2GR <- promoters(ranLocAcc2GR, upstream = 1000, downstream = 0)
} else if(orderRegion == "terminators") {
  # Obtain 1000 bp downstream of end coordinates
  source("/projects/ajt200/Rfunctions/TTSplus.R")
  ranLORAcc2GR <- TTSplus(ranLocAcc2GR, upstream = -1, downstream = 1000)
} else if(orderRegion == "genes") {
  ranLORAcc2GR <- GRanges(seqnames = as.character(seqnames(ranLocAcc2GR)),
                          ranges = IRanges(start = start(ranLocAcc2GR)-1000,
                                           end = end(ranLocAcc2GR)+1000),
                          strand = strand(ranLocAcc2GR))
} else {
  stop("orderRegion is none of bodies, promoters, terminators or genes")
}

# Count variants within orderRegion of each featuresAcc2_ortho and
# divide by orderRegion width in kb
featOrReAcc2_orthoGR_variants <- countOverlaps(query = featOrReAcc2_orthoGR,
                                               subject = variantsAcc2GR,
                                               type = "any",
                                               ignore.strand = TRUE)
featOrReAcc2_orthoGR_variantsPK <- featOrReAcc2_orthoGR_variants/(width(featOrReAcc2_orthoGR)/1000)

# Add variants per kb values to featuresAcc2_ortho dataframe
featuresAcc2_ortho <- data.frame(featuresAcc2_ortho,
                                 orderRegion_variantsPK = featOrReAcc2_orthoGR_variantsPK)

featuresAcc2_ortho$seqid <- paste0(substr(x = refbase, start = 22, stop = 24), "_", featuresAcc2_ortho$seqid)

# Group featuresAcc2_ortho into quantiles according to decreasing orderRegion_variantsPK
featuresAcc2_ortho_DF <- data.frame(featuresAcc2_ortho,
                                    percentile = rank(featuresAcc2_ortho[,which(colnames(featuresAcc2_ortho) == "orderRegion_variantsPK")]) /
                                                 length(featuresAcc2_ortho[,which(colnames(featuresAcc2_ortho) == "orderRegion_variantsPK")]),
                                    quantile = as.character(""))
for(x in as.vector(which(sapply(featuresAcc2_ortho_DF, is.list)))) {
  featuresAcc2_ortho_DF[,x] <- sapply(featuresAcc2_ortho_DF[,x], paste, collapse = "; ")
}

# Count variants within orderRegion of each ranLocAcc2GR and
# divide by orderRegion width in kb
ranLORAcc2GR_variants <- countOverlaps(query = ranLORAcc2GR,
                                       subject = variantsAcc2GR,
                                       type = "any",
                                       ignore.strand = TRUE)
ranLORAcc2GR_variantsPK <- ranLORAcc2GR_variants/(width(ranLORAcc2GR)/1000)

# Add variants per kb values to ranLocAcc2 dataframe
ranLocAcc2 <- data.frame(ranLocAcc2GR,
                         orderRegion_variantsPK = ranLORAcc2GR_variantsPK)
ranLocAcc2 <- ranLocAcc2[,-which(names(ranLocAcc2) %in% c("width", "Col_featureID"))]
names(ranLocAcc2)[names(ranLocAcc2) == "seqnames"] <- "seqid"

# Group ranLocAcc2 into quantiles according to decreasing orderRegion_variantsPK
ranLocAcc2_DF <- data.frame(ranLocAcc2,
                            percentile = rank(ranLocAcc2[,which(colnames(ranLocAcc2) == "orderRegion_variantsPK")]) /
                                         length(ranLocAcc2[,which(colnames(ranLocAcc2) == "orderRegion_variantsPK")]),
                            random = as.character(""))

quantilesStats <- data.frame()
for(k in 1:quantiles) {
  # First quantile should span 1 to greater than, e.g., 0.75 proportions of featuresAcc2_ortho
  if(k < quantiles) {
    featuresAcc2_ortho_DF[ !is.na(featuresAcc2_ortho_DF[,which(colnames(featuresAcc2_ortho_DF) == "orderRegion_variantsPK")]) &
                           rank(featuresAcc2_ortho_DF[,which(colnames(featuresAcc2_ortho_DF) == "orderRegion_variantsPK")]) /
                           length(featuresAcc2_ortho_DF[,which(colnames(featuresAcc2_ortho_DF) == "orderRegion_variantsPK")]) <=
                           1-((k-1)/quantiles) &
                           rank(featuresAcc2_ortho_DF[,which(colnames(featuresAcc2_ortho_DF) == "orderRegion_variantsPK")]) /
                           length(featuresAcc2_ortho_DF[,which(colnames(featuresAcc2_ortho_DF) == "orderRegion_variantsPK")]) >
                           1-(k/quantiles), ]$quantile <- paste0("Quantile ", k)
  } else {
  # Final quantile should span 0 to, e.g., 0.25 proportions of featuresAcc2_ortho
    featuresAcc2_ortho_DF[ !is.na(featuresAcc2_ortho_DF[,which(colnames(featuresAcc2_ortho_DF) == "orderRegion_variantsPK")]) &
                           rank(featuresAcc2_ortho_DF[,which(colnames(featuresAcc2_ortho_DF) == "orderRegion_variantsPK")]) /
                           length(featuresAcc2_ortho_DF[,which(colnames(featuresAcc2_ortho_DF) == "orderRegion_variantsPK")]) <=
                           1-((k-1)/quantiles) &
                           rank(featuresAcc2_ortho_DF[,which(colnames(featuresAcc2_ortho_DF) == "orderRegion_variantsPK")]) /
                           length(featuresAcc2_ortho_DF[,which(colnames(featuresAcc2_ortho_DF) == "orderRegion_variantsPK")]) >=
                           1-(k/quantiles), ]$quantile <- paste0("Quantile ", k)
  }
  write.table(featuresAcc2_ortho_DF[featuresAcc2_ortho_DF$quantile == paste0("Quantile ", k),],
              file = paste0(outDir,
                            "quantile", k, "_of_", quantiles,
                            "_", genomeRegion, "_by_", varType, "_freq_in_", orderRegion,
                            "_of_Acc2_Chr_genes_in_", refbase, "_",
                            paste0(chrName, collapse = "_"), ".tsv"),
              quote = FALSE, sep = "\t", row.names = FALSE)
  stats <- data.frame(quantile = as.integer(k),
                      n = as.integer(dim(featuresAcc2_ortho_DF[featuresAcc2_ortho_DF$quantile == paste0("Quantile ", k),])[1]),
                      mean_width = as.integer(round(mean(
                        (featuresAcc2_ortho_DF[featuresAcc2_ortho_DF$quantile == paste0("Quantile ", k),]$end -
                         featuresAcc2_ortho_DF[featuresAcc2_ortho_DF$quantile == paste0("Quantile ", k),]$start) + 1, na.rm = T))),
                      total_width = as.integer(sum(
                        (featuresAcc2_ortho_DF[featuresAcc2_ortho_DF$quantile == paste0("Quantile ", k),]$end -
                         featuresAcc2_ortho_DF[featuresAcc2_ortho_DF$quantile == paste0("Quantile ", k),]$start) + 1, na.rm = T)),
                      mean_orderRegion_variantsPK = as.numeric(mean(featuresAcc2_ortho_DF[featuresAcc2_ortho_DF$quantile == paste0("Quantile ", k),][,which(colnames(featuresAcc2_ortho_DF) == "orderRegion_variantsPK")], na.rm = T)))
  quantilesStats <- rbind(quantilesStats, stats)
}
write.table(quantilesStats,
            file = paste0(outDir,
                          "summary_", quantiles, "quantiles",
                          "_", genomeRegion, "_by_", varType, "_freq_in_", orderRegion,
                          "_of_Acc2_Chr_genes_in_", refbase, "_",
                          paste0(chrName, collapse = "_"), ".tsv"),
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(featuresAcc2_ortho_DF,
            file = paste0(outDir,
                          "featuresAcc2_ortho_", quantiles, "quantiles",
                          "_", genomeRegion, "_by_", varType, "_freq_in_", orderRegion,
                          "_of_Acc2_Chr_genes_in_", refbase, "_",
                          paste0(chrName, collapse = "_"), ".tsv"),
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(featuresAcc2_ortho_DF[,c(2:9, 1)],
            file = paste0(outDir,
                          "featuresAcc2_ortho_", quantiles, "quantiles",
                          "_", genomeRegion, "_by_", varType, "_freq_in_", orderRegion,
                          "_of_Acc2_Chr_genes_in_", refbase, "_",
                          paste0(chrName, collapse = "_"), ".gff3"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
featuresAcc2_ortho_DF_bed <- data.frame(chr = as.character(featuresAcc2_ortho_DF[,2]),
                                        start = as.integer(featuresAcc2_ortho_DF[,5]-1),
                                        end = as.integer(featuresAcc2_ortho_DF[,6]),
                                        name = as.character(featuresAcc2_ortho_DF[,1]),
                                        score = as.numeric(featuresAcc2_ortho_DF[,15]),
                                        strand = as.character(featuresAcc2_ortho_DF[,8]))
write.table(featuresAcc2_ortho_DF_bed,
            file = paste0(outDir,
                          "featuresAcc2_ortho_", quantiles, "quantiles",
                          "_", genomeRegion, "_by_", varType, "_freq_in_", orderRegion,
                          "_of_Acc2_Chr_genes_in_", refbase, "_",
                          paste0(chrName, collapse = "_"), ".bed"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


# Divide ranLocAcc2_DF into quantiles based on featuresAcc2_ortho_DF quantile indices
# Get row indices for each feature quantile
quantileIndices <- lapply(1:quantiles, function(k) {
  which(featuresAcc2_ortho_DF$quantile == paste0("Quantile ", k))
})
for(k in 1:quantiles) {
  ranLocAcc2_DF[quantileIndices[[k]],]$random <- paste0("Random ", k)
}
write.table(ranLocAcc2_DF,
            file = paste0(outDir,
                          "featuresAcc2_ortho_", quantiles, "quantiles",
                          "_", genomeRegion, "_by_", varType, "_freq_in_", orderRegion,
                          "_of_Acc2_Chr_genes_in_", refbase, "_",
                          paste0(chrName, collapse = "_"), "_ranLocAcc2.tsv"),
            quote = FALSE, sep = "\t", row.names = FALSE)
ranLocAcc2_DF_bed <- data.frame(chr = as.character(ranLocAcc2_DF[,1]),
                                start = as.integer(ranLocAcc2_DF[,2]-1),
                                end = as.integer(ranLocAcc2_DF[,3]),
                                name = as.integer(1:nrow(ranLocAcc2_DF)),
                                score = as.numeric(ranLocAcc2_DF[,5]),
                                strand = as.character(ranLocAcc2_DF[,4]))
write.table(ranLocAcc2_DF_bed,
            file = paste0(outDir,
                          "featuresAcc2_ortho_", quantiles, "quantiles",
                          "_", genomeRegion, "_by_", varType, "_freq_in_", orderRegion,
                          "_of_Acc2_Chr_genes_in_", refbase, "_",
                          paste0(chrName, collapse = "_"), "_ranLocAcc2.bed"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
