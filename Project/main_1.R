## Load libraries
suppressPackageStartupMessages({
  library(VariantAnnotation)
  library(GenomicFeatures)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(GenomicAlignments)
  library(CrispRVariants)
  library(dplyr)
  library(data.table)
})

source("Project/functions_1.R")

x <- import("wtc11_corrected.gtf")
v <- readVcf("deepvariant_calls_pass.vcf.gz")
alns <- "aln_s.bam"

# filter variants
keep <- geno(v)$DP[,1] > 10 &  mcols(v)$QUAL > 40 
v <- v[keep]


# test function
gene <- "PB.22."
new_seqs <- generate_variant_transcripts(v = v,x = x,
                bam_file = alns, gene = gene, verbose = TRUE)


# profile function
library(profvis)
profvis({
gene <- "PB.22."
new_seqs <- generate_variant_transcripts(v = v,x = x,
               bam_file = alns, gene = gene, verbose = TRUE)
})


# #all_genes<-unique(x$gene_id) # to get all genes
# 
#  txdb <- makeTxDbFromGRanges(x)
#  bb<-(genes(txdb))$gene_id # to get all genes
#  X<-sapply(bb,paste0,".",USE.NAMES=F) # to get all genes
#  X<-X[1:10]
# bam_file = alns
# verbose = TRUE
# #gene<-"PB.749."
# oo<-lapply(X, generate_variant_transcripts,v=v,x=x,
#           bam_file = alns, verbose = TRUE)


library(BiocParallel)
genes <- unique(x$gene_id)
system.time({
 # gvts <- lapply(genes, generate_variant_transcripts,
 #                v=v,x=x, bam_file = alns, verbose = TRUE)
 gvts <- bplapply(genes[1:24], generate_variant_transcripts,
                    v = v,x = x, bam_file = alns, verbose = TRUE,
                  BPPARAM = SerialParam())
})


# system.time({
# #  gene <- "PB.749."
#   bam_file = alns
#   verbose = TRUE
#     txdb <- makeTxDbFromGRanges(x)
#     bb<-(genes(txdb))$gene_id # to get all genes
#     X<-sapply(bb,paste0,".",USE.NAMES=F) # to get all genes
#     X<-X[1:10]
#    
#     oo<-lapply(X, generate_variant_transcripts,v=v,x=x,
#                                                 bam_file = alns, verbose = TRUE)
#     y<-do.call(c,do.call(c,oo))
#     writeXStringSet(y,"isoforms_with_variants.fasta")
#   })


# # to help testing
# gene = "PB.749."; bam_file = "aln_s.bam"; verbose = TRUE
# gene = "PB.1."; bam_file = "aln_s.bam"; verbose = TRUE
# # then step through generate_variant_transcripts()

