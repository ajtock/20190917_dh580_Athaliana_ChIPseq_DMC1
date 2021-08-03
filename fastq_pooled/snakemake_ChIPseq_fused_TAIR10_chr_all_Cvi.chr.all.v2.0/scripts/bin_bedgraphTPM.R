#!/applications/R/R-4.0.0/bin/Rscript

# This R script is called by ../Snakefile
# Apply TPM normalization to counts within
# deepTools bamCoverage-generated per-base bedgraph files

# Usage:
# ./bin_bedgraphTPM.R WT_DMC1_V5_Rep2_ChIP fused_TAIR10_chr_all_Cvi.chr.all.v2.0 both 1

#sampleName <- "WT_DMC1_V5_Rep2_ChIP"
#refbase <- "fused_TAIR10_chr_all_Cvi.chr.all.v2.0"
#align <- "both"
#binSize <- 1

args <- commandArgs(trailingOnly = T)
sampleName <- args[1]
refbase <- args[2]
align <- args[3]
binSize <- as.integer(args[4])

options(stringsAsFactors = F)
library(parallel)
library(plyr)
library(data.table)
library(yaml)
#library(varhandle)
#library(zoo)

# Genomic definitions
fai <- read.table(paste0("data/index/", refbase, ".fa.fai"), header = F)
config <- read_yaml("config.yaml")
ignoreForNormalization <- unlist(strsplit(config$COVERAGE$ignoreForNormalization,
                                          split = " "))
fai <- fai[!(fai$V1 %in% ignoreForNormalization),]
if(!grepl("Chr", fai[,1][1])) {
  chrs <- paste0("Chr", fai[,1])
} else {
  chrs <- fai[,1]
}
chrLens <- fai[,2]

# Make chromosomal coordinates cumulative
# such that the first coordinate of Chr2 is
# equal to the last coordinate of Chr1 + 1
sumchr <- cumsum(c(0, chrLens))
print(sumchr)

# Load bedgraph
sampleProfile <- read.table(paste0("mapped/", align, "/bg/",
                                   sampleName, "_MappedOn_", refbase, "_lowXM_",
                                   align, "_sort_dtnorm.bedgraph"))
sampleProfile <- sampleProfile[!(sampleProfile$V1 %in% ignoreForNormalization),]
if(!grepl("Chr", fai[,1][1])) {
  sampleProfile$V1 <- paste0("Chr", sampleProfile$V1)
}

# Rows where the difference between end and start coordinates is > binSize
sampleProfile_bigWins <- sampleProfile[sampleProfile$V3-sampleProfile$V2 > binSize,]
# Rows where the difference between end and start coordinates is == binSize
sampleProfile <- sampleProfile[sampleProfile$V3-sampleProfile$V2 == binSize,]

# Create a list of big windows, each split into windows of binSize,
# or < binSize if at chromosome end
sampleProfile_bigWinsList <- mclapply(seq_along(1:dim(sampleProfile_bigWins)[1]), function(x) {
  bigWinsSplit <- seq(from = sampleProfile_bigWins[x,]$V2,
                      to = sampleProfile_bigWins[x,]$V3,
                      by = binSize)

  if(bigWinsSplit[length(bigWinsSplit)] < sampleProfile_bigWins[x,]$V3) {
    data.frame(V1 = as.character(sampleProfile_bigWins[x,]$V1),
               V2 = as.integer(c(bigWinsSplit[-length(bigWinsSplit)],
                                 bigWinsSplit[length(bigWinsSplit)])),
               V3 = as.integer(c(bigWinsSplit[-length(bigWinsSplit)]+binSize,
                                 sampleProfile_bigWins[x,]$V3)),
               V4 = as.numeric(sampleProfile_bigWins[x,]$V4))
  } else if (bigWinsSplit[length(bigWinsSplit)] == sampleProfile_bigWins[x,]$V3) {
    data.frame(V1 = as.character(sampleProfile_bigWins[x,]$V1),
               V2 = as.integer(bigWinsSplit[-length(bigWinsSplit)]),
               V3 = as.integer(bigWinsSplit[-length(bigWinsSplit)]+binSize),
               V4 = as.numeric(sampleProfile_bigWins[x,]$V4))
  }
}, mc.cores = detectCores())

sampleProfile_bigWinsDT <- rbindlist(sampleProfile_bigWinsList)
sampleProfile <- rbind.fill(sampleProfile, sampleProfile_bigWinsDT)
sampleProfile <- sampleProfile[order(sampleProfile$V1, sampleProfile$V2),]

chrLenValsList <- mclapply(seq_along(chrs), function (x) {
  chrProfileSample <- sampleProfile[sampleProfile$V1 == chrs[x],]
  if(chrProfileSample[dim(chrProfileSample)[1],]$V3 < chrLens[x]) {
    data.frame(V1 = chrs[x],
               V2 = as.integer(chrProfileSample[dim(chrProfileSample)[1],]$V3),
               V3 = as.integer(chrLens[x]),
               V4 = as.numeric(chrProfileSample[dim(chrProfileSample)[1],]$V4))
  }
}, mc.cores = length(chrs))
sampleProfile_chrLenValsDT <- rbindlist(chrLenValsList)
sampleProfile <- rbind.fill(sampleProfile, sampleProfile_chrLenValsDT)
sampleProfile <- sampleProfile[order(sampleProfile$V1, sampleProfile$V2),]

## Calculate TPM
## See https://rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/
# Divide read counts by the length of each bin in kb (reads per kilobase, RPK)
RPK <- sampleProfile$V4 / (binSize/1e3)
# Sum all the RPK values and divide this number by 1e6
# (the "per million" scaling factor)
scaling_factor <- sum(RPK) / 1e6
# Divide each RPK value by scaling_factor to give TPM
TPM <- RPK / scaling_factor
print("Sum of per-base TPM =")
print(sum(TPM))

bedgraph <- data.frame(chr = as.character(sampleProfile$V1),
                       start0based = as.integer(sampleProfile$V2),
                       end = as.integer(sampleProfile$V3),
                       TPM = as.numeric(TPM))
write.table(bedgraph,
            file = paste0("mapped/", align, "/bg/",
                          sampleName, "_MappedOn_", refbase,
                          "_lowXM_", align, "_sort_TPM.bedgraph"),
            sep = "\t", quote = F, row.names = F, col.names = F)
