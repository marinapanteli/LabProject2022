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


# TESTS: Remember to change gene and run the 4 following lines of code before test!
gene <- "PB.749"  # Select the correct gene
new_seqs <- generate_variant_transcripts(v = v,x = x,
                                         bam_file = alns, gene = gene, verbose = TRUE)
x_ex <- x[x$gene_id==gene]
x_exs <- split(x_ex, x_ex$transcript_id)




####### 


# Careful, in case an insertion or a deletion has happened on the left of this area of variation. That is why there might be a shift.
# ### SNP - 
a <- new_seqs$PB.22.1.a
b <- new_seqs$PB.22.1.b

tr_gr <- x_ex[x_ex$transcript_id == "PB.22.1"]
gr <- GRanges(seqnames = seqnames(tr_gr),
              IRanges(ranges(tr_gr)), strand = "+")

dss_ref <- getSeq(Hsapiens, gr,
                  as.character = TRUE)
ref_seq <- DNAStringSet(paste(dss_ref, collapse = ""))

alleles_of_interest<-DNAStringSet(c(as.character(a), as.character(b)))
loc1<-str_locate(as.character(ref_seq),"ATGTAGATGGGCCCGTC")[1]
loc2<-str_locate(as.character(ref_seq),"ATGTAGATGGGCCCGTC")[2]

substring(as.character(reverseComplement(alleles_of_interest)),loc1,loc2)
# ### SNP +
# Run for PB.749.3 (which is i=3). GGCCCGGATGAGCAGACTCCTGT is in the ref_seq
#substring(as.character((new_seqs[[i]])),str_locate(as.character(ref_seq),"GGCCCGGATGAGCAGACTCCTGT")[1],str_locate(as.character(ref_seq),"GGCCCGGATGAGCAGACTCCTGT")[2])
a <- new_seqs$PB.749.3.a
b <- new_seqs$PB.749.3.b

tr_gr <- x_ex[x_ex$transcript_id == "PB.749.3"]
gr <- GRanges(seqnames = seqnames(tr_gr),
              IRanges(ranges(tr_gr)), strand = "+")

dss_ref <- getSeq(Hsapiens, gr,
                  as.character = TRUE)
ref_seq <- DNAStringSet(paste(dss_ref, collapse = ""))

alleles_of_interest<-DNAStringSet(c(as.character(a), as.character(b)))
loc1<-str_locate(as.character(ref_seq),"GGCCCGGATGAGCAGACTCCTGT")[1]
loc2<-str_locate(as.character(ref_seq),"GGCCCGGATGAGCAGACTCCTGT")[2]

substring(as.character((alleles_of_interest)),loc1,loc2)

# ### DEL -  
a <- new_seqs$PB.22.32.a
b <- new_seqs$PB.22.32.b

tr_gr <- x_ex[x_ex$transcript_id == "PB.22.32"]
gr <- GRanges(seqnames = seqnames(tr_gr),
              IRanges(ranges(tr_gr)), strand = "+")

dss_ref <- getSeq(Hsapiens, gr,
                  as.character = TRUE)
ref_seq <- DNAStringSet(paste(dss_ref, collapse = ""))

alleles_of_interest<-DNAStringSet(c(as.character(a), as.character(b)))
loc1<-str_locate(as.character(ref_seq),"CCTGGCTGCTGGGGAGGAC")[1]
loc2<-str_locate(as.character(ref_seq),"CCTGGCTGCTGGGGAGGAC")[2]

substring(as.character(reverseComplement(alleles_of_interest)),loc1,loc2)
# 
# 
# ### DEL +
a <- new_seqs$PB.749.3.a
b <- new_seqs$PB.749.3.b

tr_gr <- x_ex[x_ex$transcript_id == "PB.749.3"]
gr <- GRanges(seqnames = seqnames(tr_gr),
              IRanges(ranges(tr_gr)), strand = "+")

dss_ref <- getSeq(Hsapiens, gr,
                  as.character = TRUE)
ref_seq <- DNAStringSet(paste(dss_ref, collapse = ""))

alleles_of_interest<-DNAStringSet(c(as.character(a), as.character(b)))
loc1<-str_locate(as.character(ref_seq),"AAATGAAAAACGTTTGCTAGA")[1]
loc2<-str_locate(as.character(ref_seq),"AAATGAAAAACGTTTGCTAGA")[2]

substring(as.character((alleles_of_interest)),loc1,loc2)
# 
# ### INS -  
a <- new_seqs$PB.1.1.a
b <- new_seqs$PB.1.1.b

tr_gr <- x_ex[x_ex$transcript_id == "PB.1.1"]
gr <- GRanges(seqnames = seqnames(tr_gr),
              IRanges(ranges(tr_gr)), strand = "+")

dss_ref <- getSeq(Hsapiens, gr,
                  as.character = TRUE)
ref_seq <- DNAStringSet(paste(dss_ref, collapse = ""))

alleles_of_interest<-DNAStringSet(c(as.character(a), as.character(b)))
loc1<-str_locate(as.character(ref_seq),"CAGAGTGGCCAGCCAC")[1]
loc2<-str_locate(as.character(ref_seq),"CAGAGTGGCCAGCCAC")[2]

substring(as.character(reverseComplement(alleles_of_interest)),loc1,loc2)

# 
# 
# ### INS +
a <- new_seqs$PB.4.3.a
b <- new_seqs$PB.4.3.b

tr_gr <- x_ex[x_ex$transcript_id == "PB.4.3"]
gr <- GRanges(seqnames = seqnames(tr_gr),
              IRanges(ranges(tr_gr)), strand = "+")

dss_ref <- getSeq(Hsapiens, gr,
                  as.character = TRUE)
ref_seq <- DNAStringSet(paste(dss_ref, collapse = ""))

alleles_of_interest<-DNAStringSet(c(as.character(a), as.character(b)))
loc1<-str_locate(as.character(ref_seq),"TGCACACACGAGCA")[1]
loc2<-str_locate(as.character(ref_seq),"TGCACACACGAGCA")[2]

substring(as.character((alleles_of_interest)),loc1,loc2)
















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
