var_8_3<-function(v,x,gene){

reads <- "aln_s.bam"    #test this here, not in main?
x_ex<-extract_exons(x)  #A
chromosome<-as.character(seqnames(x_ex))[1]

gr <- GRanges(seqnames(x_ex)[1],IRanges(start = min(start(x_ex)),end = max(end(x_ex))),strand = strand(x_ex)[1])

(fop <- findOverlapPairs(rowRanges(v), x_ex))

sbp <- ScanBamParam(which=gr, what=scanBamWhat())
bamWhat(sbp) <- "seq"
ga <- readGAlignments(reads, param=sbp, use.names = TRUE)
(variation_regs <- unique(first(fop)))

rngs <- get_read_ranges(ga) #B

x_exs <- split(x_ex, x_ex$transcript_id)

# intersect all reads with all transcripts, form matrix
m <- matrix(0, ncol=length(variation_regs), nrow = length(rngs), dimnames = list(names(rngs),names(variation_regs)))
rngs_u <- unlist(rngs)
fo <- findOverlaps(rngs_u, ranges(variation_regs))
pop <- split(subjectHits(fo), names(rngs_u)[queryHits(fo)]) # which variations in a certain alignment read, and below notate with 0/1 if absent/present
for(i in 1:length(pop)) {
  m[names(pop)[i],unique(pop[[i]])] <- 1
}

new_seqs <- vector()
for (k in second(fop)$transcript_id) {
  
  f <- first(fop)[second(fop)$transcript_id==k]
  (nm <- names(f))
  
  # only take reads that overlap all variation of this transcript
  w <- rowSums(m[,nm,drop=FALSE])==length(nm) 
  (rds <- rownames(m)[w]) 
  
  # extract transcript seq from genome + positions where we must insert variation (transcript coord.)
  dss_ref <- getSeq(Hsapiens, x_ex[x_ex$transcript_id==k])
  ref_seq <- DNAStringSet(paste(as.character(dss_ref), collapse=""))
  names(ref_seq) <- k
  (mtt <- mapToTranscripts(variation_regs[nm], x_exs[k]))
  st <- start(ranges(variation_regs[nm]))
  
  df <- narrowal(st,rds,ga,chromosome) #C
  to_insert <- get_alleles(df) #D
  new_seqs <- seqs_with_var(to_insert,ref_seq,mtt,x_ex,k,new_seqs) #E
}
new_seqs<-do.call(c, new_seqs)
}
