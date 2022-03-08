## Load libraries


suppressPackageStartupMessages({
  library(VariantAnnotation)
  library(GenomicFeatures)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(Biostrings)
  library(ensembldb)
  library(Rsamtools)
  library(GenomicAlignments)
  library(CrispRVariants)
  library(stringr)
})

source("./functions.R")
source("./var_8_3.R")

(v <- readVcf("deepvariant_calls_pass.vcf.gz"))
keep <- geno(v)$DP[,1] > 10 &  mcols(v)$QUAL > 40 # keep only variants with depth > 10
v <- v[keep]

gene<-"PB.749."  #Generalize this!
x <- import("wtc11_corrected.gtf")
txdb <- makeTxDbFromGRanges(x)

new_seqs<-var_8_3(v,x,gene) #for 1 gene
  