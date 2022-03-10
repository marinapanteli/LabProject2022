## Load libraries
suppressPackageStartupMessages({
  library(VariantAnnotation)
  library(GenomicFeatures)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(GenomicAlignments)
  library(CrispRVariants)
})

source("functions_1.R")

x <- import("wtc11_corrected.gtf")
v <- readVcf("deepvariant_calls_pass.vcf.gz")
alns <- "aln_s.bam"

# filter variants
keep <- geno(v)$DP[,1] > 10 &  mcols(v)$QUAL > 40 
v <- v[keep]


# library(profvis)
# profvis({
# gene <- "PB.749."  #Generalize this!
# new_seqs <- generate_variant_transcripts(v,x,
#                     bam_file = alns, gene, verbose = FALSE) # for 1 gene
# })


#all_genes<-unique(x$gene_id) # to get all genes

system.time({
#  gene <- "PB.749."
    txdb <- makeTxDbFromGRanges(x)
    bb<-(genes(txdb))$gene_id # to get all genes
    X<-sapply(bb,paste0,".",USE.NAMES=F) # to get all genes
    oo<-lapply(X, generate_variant_transcripts,v=v,x=x,
                                                bam_file = alns, verbose = TRUE)
#    new_seqs <- generate_variant_transcripts(v,x,
#                                           bam_file = alns, gene, verbose = TRUE) # for 1 gene
  #old_seqs <- generate_old_transcripts(v,x,
   #                                        bam_file = alns, gene, verbose = FALSE) # for 1 gene
  
  })


# # to help testing
# gene = "PB.749."; bam_file = "aln_s.bam"; verbose = TRUE
# # then step through generate_variant_transcripts()

