#!/applications/R/R-4.0.0/bin/Rscript

# author: Andy Tock
# contact: ajt200@cam.ac.uk
# date: 04.08.2021

# Profile varType variant frequency around genes and random loci

# Usage:
# /applications/R/R-4.0.0/bin/Rscript variant_freq_profiles_around_genes.R SNP_INDEL '/data/public_data/arabidopsis/MPIPZ_Jiao_Schneeberger_2020_NatCommun/Cvi' unique 'Chr1,Chr2,Chr3,Chr4,Chr5' 2200 2000 2kb '2 kb' 10 10bp genes genomewide 6 fused_TAIR10_chr_all_Cvi.chr.all.v2.0 '0.02,0.96'

#varType <- "SNP_INDEL"
#dirName <- "/data/public_data/arabidopsis/MPIPZ_Jiao_Schneeberger_2020_NatCommun/Cvi"
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

args <- commandArgs(trailingOnly = T)
varType <- args[1]
dirName <- args[2]
chrName <- unlist(strsplit(args[3],
                               split = ","))
bodyLength <- as.numeric(args[4])
upstream <- as.numeric(args[6])
downstream <- as.numeric(args[6])
flankName <- args[7]
binSize <- as.numeric(args[9])
orderRegion <- args[11]
genomeRegion <- args[12]
quantiles <- as.numeric(args[13])
refbase <- args[14]

args <- commandArgs(trailingOnly = T)
featureName <- args[1]
region <- args[2]
bodyLength <- as.numeric(args[3])
flankSize <- as.numeric(args[4])
flankName <- as.character(args[5])
binSize <- as.numeric(args[6])

# Source functions to be used in this script
source("/projects/ajt200/Rfunctions/covMatrix_DNAmethMatrix_target_ranLoc_R3.4.0.R")

library(EnrichedHeatmap)

matDir <- paste0("matrices/")
system(paste0("[ -d ", matDir, " ] || mkdir ", matDir))

# Genomic definitions
chrs <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,1])
chrs <- chrs[-length(chrs)]
chrStart <- c(rep(1, times = length(chrs)))
chrLens <- as.vector(read.table("/home/ajt200/analysis/wheat/sRNAseq_meiocyte_Martin_Moore/snakemake_sRNAseq/data/index/wheat_v1.0.fa.sizes")[,2])
chrLens <- chrLens[-length(chrLens)]
centromereStart <- as.vector(read.table(paste0("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/",
                                               "centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11_chr4ALeftmostInterval_chr5ARightTwoIntervals.txt"))[,2])
centromereEnd <- as.vector(read.table(paste0("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/",
                                             "centromeres_outer_CENH3enriched_IWGSC_2018_Science_Table_S11_chr4ALeftmostInterval_chr5ARightTwoIntervals.txt"))[,3])
chrPartitions <- read.table("/home/ajt200/analysis/wheat/wheat_IWGSC_WGA_v1.0_pseudomolecules/chromosome_partitions_IWGSC_2018_Science_Table_S29.txt",
                            header = TRUE)
genomeGR <- GRanges(seqnames = chrs,
                    ranges = IRanges(start = chrStart,
                                     end = chrLens),
                    strand = "*")
genomeGR <- genomeGR[grep(substr(featureName, 10, 10),
                          seqnames(genomeGR))@values]

# Define region to be analysed
if(region == "euchromatin") {
  regionGR <- GRanges(seqnames = rep(chrPartitions$chrom, 2),
                      ranges = IRanges(start = c(rep(1, dim(chrPartitions)[1]),
                                                 chrPartitions$R2b_R3),
                                       end = c(chrPartitions$R1_R2a,
                                               chrLens)),
                      strand = "*")
  regionGR <- regionGR[grep(substr(featureName, 10, 10),
                            seqnames(regionGR))@values]
} else if(region == "interstitial") {
  regionGR <- GRanges(seqnames = rep(chrPartitions$chrom, 2),
                      ranges = IRanges(start = c(chrPartitions$R1_R2a+1,
                                                 chrPartitions$C_R2b),
                                       end = c(chrPartitions$R2a_C,
                                               chrPartitions$R2b_R3-1)),
                      strand = "*")
  regionGR <- regionGR[grep(substr(featureName, 10, 10),
                            seqnames(regionGR))@values]
} else if(region == "proximal") {
  regionGR <- GRanges(seqnames = chrPartitions$chrom,
                      ranges = IRanges(start = chrPartitions$R2a_C+1,
                                       end = chrPartitions$C_R2b-1),
                      strand = "*")
  regionGR <- regionGR[grep(substr(featureName, 10, 10),
                            seqnames(regionGR))@values]
} else if(region == "heterochromatin") {
  regionGR <- GRanges(seqnames = chrPartitions$chrom,
                      ranges = IRanges(start = chrPartitions$R1_R2a+1,
                                       end = chrPartitions$R2b_R3-1),
                      strand = "*")
  regionGR <- regionGR[grep(substr(featureName, 10, 10),
                            seqnames(regionGR))@values]
} else if(region == "centromeres") {
  regionGR <- GRanges(seqnames = chrs,
                      ranges = IRanges(start = centromereStart,
                                       end = centromereEnd),
                      strand = "*")
  regionGR <- regionGR[grep(substr(featureName, 10, 10),
                            seqnames(regionGR))@values]
} else if(region == "genomewide") {
  regionGR <- GRanges(seqnames = chrs,
                      ranges = IRanges(start = rep(1, length(chrs)),
                                       end = chrLens),
                      strand = "*")
  regionGR <- regionGR[grep(substr(featureName, 10, 10),
                            seqnames(regionGR))@values]
} else {
  stop("region is not euchromatin, interstitial, proximal, heterochromatin, centromeres, or genomewide")
}

# Define region to be masked out of analysis
if(region == "euchromatin") {
  maskGR <- GRanges(seqnames = chrPartitions$chrom,
                    ranges = IRanges(start = chrPartitions$R1_R2a+1,
                                     end = chrPartitions$R2b_R3-1),
                    strand = "*")
  maskGR <- maskGR[grep(substr(featureName, 10, 10),
                        seqnames(maskGR))@values]
} else if(region == "interstitial") {
  maskGR <- GRanges(seqnames = rep(chrPartitions$chrom, 3),
                    ranges = IRanges(start = c(rep(1, dim(chrPartitions)[1]),
                                               chrPartitions$R2a_C+1,
                                               chrPartitions$R2b_R3),
                                     end = c(chrPartitions$R1_R2a,
                                             chrPartitions$C_R2b-1,
                                             chrLens)),
                    strand = "*")
  maskGR <- maskGR[grep(substr(featureName, 10, 10),
                        seqnames(maskGR))@values]
} else if(region == "proximal") {
  maskGR <- GRanges(seqnames = rep(chrPartitions$chrom, 2),
                    ranges = IRanges(start = c(rep(1, dim(chrPartitions)[1]),
                                               chrPartitions$C_R2b),
                                     end = c(chrPartitions$R2a_C,
                                             chrLens)),
                    strand = "*")
  maskGR <- maskGR[grep(substr(featureName, 10, 10),
                        seqnames(maskGR))@values]
} else if(region == "heterochromatin") {
  maskGR <- GRanges(seqnames = rep(chrPartitions$chrom, 2),
                    ranges = IRanges(start = c(rep(1, dim(chrPartitions)[1]),
                                               chrPartitions$R2b_R3),
                                     end = c(chrPartitions$R1_R2a,
                                             chrLens)),
                    strand = "*")
  maskGR <- maskGR[grep(substr(featureName, 10, 10),
                        seqnames(maskGR))@values]
} else if(region == "centromeres") {
  maskGR <- GRanges(seqnames = rep(chrPartitions$chrom, 2),
                    ranges = IRanges(start = c(rep(1, dim(chrPartitions)[1]),
                                               centromereEnd+1),
                                     end = c(centromereStart-1,
                                             chrLens)),
                    strand = "*")
  maskGR <- maskGR[grep(substr(featureName, 10, 10),
                        seqnames(maskGR))@values]
} else if(region == "genomewide") {
  maskGR <- GRanges()
  maskGR <- maskGR[grep(substr(featureName, 10, 10),
                        seqnames(maskGR))@values]
} else {
  stop("region is not euchromatin, interstitial, proximal, heterochromatin, centromeres, or genomewide")
}

# Load genes in BED format and convert into GRanges
# Note addition of 1 to 0-based BED start coordinates
genes <- read.table(paste0("/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706_representative_mRNA_in_",
                           substring(featureName, first = 10, last = 16), "_", region, ".bed"),
                    header = F)
colnames(genes) <- c("chr", "start", "end", "name", "score", "strand")
genesGR <- GRanges(seqnames = genes$chr,
                   ranges = IRanges(start = genes$start+1,
                                    end = genes$end),
                   strand = genes$strand,
                   number = genes$name)
genesGR <- genesGR[seqnames(genesGR) != "chrUn"]
# Subset to include only those not overlapping masked region
mask_genes_overlap <- findOverlaps(query = maskGR,
                                   subject = genesGR,
                                   type = "any",
                                   select = "all",
                                   ignore.strand = TRUE)
if(length(mask_genes_overlap) > 0) {
  genesGR <- genesGR[-subjectHits(mask_genes_overlap)]
}
# Load ranLoc in BED format and convert into GRanges
# Note addition of 1 to 0-based BED start coordinates
ranLoc <- read.table(paste0("/home/ajt200/analysis/wheat/annotation/221118_download/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706_representative_mRNA_in_",
                            substring(featureName, first = 10, last = 16), "_", region, "_randomLoci.bed"),
                     header = F)
colnames(ranLoc) <- c("chr", "start", "end", "name", "score", "strand")
ranLocGR <- GRanges(seqnames = ranLoc$chr,
                    ranges = IRanges(start = ranLoc$start+1,
                                     end = ranLoc$end),
                    strand = ranLoc$strand,
                    number = ranLoc$name)
ranLocGR <- ranLocGR[seqnames(ranLocGR) != "chrUn"]
# Subset to include only those not overlapping masked region
mask_ranLoc_overlap <- findOverlaps(query = maskGR,
                                    subject = ranLocGR,
                                    type = "any",
                                    select = "all",
                                    ignore.strand = TRUE)
if(length(mask_ranLoc_overlap) > 0) {
  ranLocGR <- ranLocGR[-subjectHits(mask_ranLoc_overlap)]
}

# Load SNPs
SNPs <- read.table(paste0("/home/ajt200/analysis/wheat/annotation/He_Akhunov_2019_NatGenet_1000exomes_SNPs/",
                          "all.GP08_mm75_het3_publication01142019.ann.vcf"),
                   header = F, skip = 31,
                   colClasses = c(rep(NA, 2),
                                  "NULL",
                                  rep(NA, 2),
                                  rep("NULL", 2),
                                  NA))
# Add this line when working with complete VCF including per-accession genotype columns
#                                  rep("NULL", 812)))
colnames(SNPs) <- c("chr", "pos", "ref", "alt", "info")
## all exome SNPs
SNPsGR <- GRanges(seqnames = SNPs$chr,
                  ranges = IRanges(start = SNPs$pos,
                                   end = SNPs$pos),
                  strand = "*",
                  coverage = rep(1, dim(SNPs)[1]))

# Subset SNPs by class (grep-ing for "A" within VCF "info" field returns all SNPs)
SNPclass <- c(
              "A",
              "upstream_gene_variant",
              "downstream_gene_variant",
              "missense_variant",
              "synonymous_variant",
              "HIGH",
              "MODERATE",
              "LOW",
              "MODIFIER",
              "intron_variant",
              "intergenic"
             )
SNPsListGR <- mclapply(seq_along(SNPclass), function(x) {
  classSNPs <- SNPs[grep(SNPclass[x], SNPs$info),]
  GRanges(seqnames = classSNPs$chr,
          ranges = IRanges(start = classSNPs$pos,
                           end = classSNPs$pos),
          strand = "*",
          coverage = rep(1, dim(classSNPs)[1]))
}, mc.cores = length(SNPclass))

## transition
SNPs_transition <- SNPs[(SNPs$ref == "A" | SNPs$ref == "G") & (SNPs$alt == "G" | SNPs$alt == "A") |
                        (SNPs$ref == "C" | SNPs$ref == "T") & (SNPs$alt == "T" | SNPs$alt == "C"),]
SNPs_transition_GR <- GRanges(seqnames = SNPs_transition$chr,
                              ranges = IRanges(start = SNPs_transition$pos,
                                               end = SNPs_transition$pos),
                              strand = "*",
                              coverage = rep(1, dim(SNPs_transition)[1]))

## transversion
SNPs_transversion <- SNPs[(SNPs$ref == "A" | SNPs$ref == "G") & (SNPs$alt == "C" | SNPs$alt == "T") |
                          (SNPs$ref == "C" | SNPs$ref == "T") & (SNPs$alt == "A" | SNPs$alt == "G"),]
stopifnot((dim(SNPs_transition)[1] +
          dim(SNPs_transversion)[1]) ==
          dim(SNPs)[1])
SNPs_transversion_GR <- GRanges(seqnames = SNPs_transversion$chr,
                                ranges = IRanges(start = SNPs_transversion$pos,
                                                 end = SNPs_transversion$pos),
                                strand = "*",
                                coverage = rep(1, dim(SNPs_transversion)[1]))

# Add transitions and transversions GRanges to SNPsListGR
SNPsListGR <- c(SNPsListGR,
                SNPs_transition_GR,
                SNPs_transversion_GR)
SNPclassNames <- c(
                   "all",
                   "upstream_gene_variant",
                   "downstream_gene_variant",
                   "missense_variant",
                   "synonymous_variant",
                   "HIGH",
                   "MODERATE",
                   "LOW",
                   "MODIFIER",
                   "intron_variant",
                   "intergenic",
                   "transition",
                   "transversion"
                  )

# Define matrix and column mean outfiles
outDF <- lapply(seq_along(SNPclassNames), function(x) {
  list(paste0(matDir, "exome_", SNPclassNames[x],
              "_SNPs_around_", featureName, "_", region,
              "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
       paste0(matDir, "exome_", SNPclassNames[x],
              "_SNPs_around_", featureName, "_", region,
              "_ranLoc_matrix_bin", binSize, "bp_flank", flankName, ".tab"))
})
outDFcolMeans <- lapply(seq_along(SNPclassNames), function(x) {
  list(paste0(matDir, "exome_", SNPclassNames[x],
              "_SNPs_around_", featureName, "_", region,
              "_matrix_bin", binSize, "bp_flank", flankName, "_colMeans.tab"),
       paste0(matDir, "exome_", SNPclassNames[x],
              "_SNPs_around_", featureName, "_", region,
              "_ranLoc_matrix_bin", binSize, "bp_flank", flankName, "_colMeans.tab"))
})

# Run covMatrix() function on each feature GRanges object to obtain matrices
# containing normalised feature density values around target and random loci
mclapply(seq_along(SNPclassNames), function(x) {
  covMatrix(signal = SNPsListGR[[x]],
            feature = genesGR,
            ranLoc = ranLocGR,
            featureSize = bodyLength,
            flankSize = flankSize,
            binSize = binSize,
            outDF = outDF[[x]],
            outDFcolMeans = outDFcolMeans[[x]])
  print(paste0("exome ", SNPclassNames[x],
               " SNP frequency around ", featureName, " ", region,
               " profile calculation complete"))
}, mc.cores = length(SNPclassNames))
