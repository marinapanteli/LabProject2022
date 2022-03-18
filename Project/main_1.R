## Load libraries
suppressPackageStartupMessages({
  library(VariantAnnotation)
  library(GenomicFeatures)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(GenomicAlignments)
  library(CrispRVariants)
  library(dplyr)
  library(stringr) 
})

source("functions_1.R")

x <- import("wtc11_corrected.gtf")
x <- extract_exons(x)

v <- readVcf("deepvariant_calls_pass.vcf.gz")
keep <- geno(v)$DP[,1] > 10 &  mcols(v)$QUAL > 40 
v <- v[keep]

alns <- "aln_s.bam"

# filter variants


#gene <- "PB.22"  # Select the correct gene
#new_seqs <- generate_variant_transcripts(v = v,x = x,
#                                         bam_file = alns, gene = gene, verbose = TRUE)


(test_variant_seqs_1 <- c("SNP -",test_variant(x, gene="PB.22", which_transcript="PB.22.1" ,ex_seq="ATGTAGATGGGCCCGTC" , v = v,x = x,
                             bam_file = alns, verbose = FALSE), "the ref seq was ATGTAGATGGGCCCGTC and we know: chr1:1319461_C/G"))
  
(test_variant_seqs_2 <- c("SNP +",test_variant(x, gene="PB.749", which_transcript="PB.749.1" ,ex_seq="GGCCCGGATGAGCAGACTCCTGT" , v = v,x = x,
                                               bam_file = alns, verbose = FALSE), "the ref seq was GGCCCGGATGAGCAGACTCCTGT and we know: chr1:95117534_G/C"))

(test_variant_seqs_3 <- c("DEL -",test_variant(x, gene="PB.22", which_transcript="PB.22.32" ,ex_seq="CCTGGCTGCTGGGGAGGAC" , v = v,x = x,
                                                                                              bam_file = alns, verbose = FALSE), "the ref seq was CCTGGCTGCTGGGGAGGAC and we know chr1:1320023_TG/T"))

(test_variant_seqs_4 <- c("DEL +",test_variant(x, gene="PB.749", which_transcript="PB.749.3" ,ex_seq="AAATGAAAAACGTTTGCTAGA" , v = v,x = x,
                                               bam_file = alns, verbose = FALSE), "the ref seq was AAATGAAAAACGTTTGCTAGA and we know: chr1:95194468_AC/A"))

(test_variant_seqs_5 <- c("INS -",test_variant(x, gene="PB.1", which_transcript="PB.1.1" ,ex_seq="CAGAGTGGCCAGCCAC" , v = v,x = x,
                                               bam_file = alns, verbose = FALSE), "the ref seq was CAGAGTGGCCAGCCAC and we know : chr1:15903_G/GC"))

(test_variant_seqs_6 <- c("INS +",test_variant(x, gene="PB.4", which_transcript="PB.4.3" ,ex_seq="TGCACACACGAGCA" , v = v,x = x,
                                               bam_file = alns, verbose = FALSE), "the ref seq was TGCACACACGAGCA and we know: chr1:855316_C/CAT"))









#gene <- unique(x$gene_id)
new_seqs <- generate_variant_transcripts(v = v,x = x,
                                         bam_file = alns, gene = gene, verbose = TRUE)

new_seqs <- do.call(c,new_seqs)
# 
# genes <- unique(x$gene_id)
# X<-genes
# X<-X[1:10]

# profile function
library(profvis)
library(BiocParallel)
profvis({
  # gene <- "PB.22"
  # new_seqs <- generate_variant_transcripts(v = v,x = x,
  #                bam_file = alns, gene = gene, verbose = FALSE)
  genes <- unique(x$gene_id)
  gvts <- bplapply(genes[1:24], generate_variant_transcripts,
                   v = v,x = x, bam_file = alns, verbose = TRUE,
                   BPPARAM = SerialParam())
  gvts <- do.call(c,gvts)
  
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
# oo<-lapply(X, generate_variant_transcripts,v=v,x=x,
#                       bam_file = alns, verbose = TRUE)
# 
library(BiocParallel)
genes <- unique(x$gene_id)
system.time({
  # gvts <- lapply(genes, generate_variant_transcripts,
  #                v=v,x=x, bam_file = alns, verbose = TRUE)
  gvts <- bplapply(genes[1:24], generate_variant_transcripts,
                   v = v,x = x, bam_file = alns, verbose = TRUE,
                   BPPARAM = SerialParam())
  gvts <- do.call(c,gvts)
})


# system.time({
# #  gene <- "PB.749"
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
# gene = "PB.749"; bam_file = "aln_s.bam"; verbose = TRUE
# gene = "PB.1"; bam_file = "aln_s.bam"; verbose = TRUE
# gene = "PB.22"; bam_file = "aln_s.bam"; verbose = TRUE
# # then step through generate_variant_transcripts()
