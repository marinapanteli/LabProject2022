get_read_ranges <- function(ga) {
  rngs <- cigarRangesAlongReferenceSpace(cigar(ga), 
                                         pos = start(ga),
                                         drop.empty.ranges = FALSE,
                                         N.regions.removed = FALSE,
                                         with.ops = TRUE,
                                         reduce.ranges = FALSE)
  names(rngs) <- names(ga)
  # remove 'N' ranges
  rngs <- lapply(rngs, function(u) {
    nm <- names(u)
    names(u) <- NULL
    u[nm != "N"]
  })
  IRangesList(rngs)
}


extract_exons <- function(x, gene) {
  x_ex <- x[x$type=="exon" & seqnames(x) != "chrM"]
  g <- grepl(gene, x_ex$transcript_id, 
             fixed = TRUE)
  x_ex <- x_ex[g]
  w <- which(seqlevels(x_ex)=="chrM")
  seqlevels(x_ex) <- seqlevels(x_ex)[-w]
  x_ex
}

narrowal <- function(variation, ga){
  
  starts <- start(variation)
  names(starts) <- names(variation)
  
  lapply(starts, function(u) {
    ir <- IRanges(start = u, width = 2)
    gr <- GRanges(seqnames(ga)[1],
                  ir, strand = "+")
    na <- narrowAlignments(ga, gr)
    #as.character(mcols(na)$seq)
    mcols(na)
  })
  
}


seqs_with_var <- function(to_insert,ref_seq,mtt,x_ex,tr){
  seqs <- list()
  to_insert <- to_insert[lengths(to_insert)!=0]
  for(i in 1:length(to_insert)) { # loop through alleles
    seqs[[i]] <- ref_seq
    for(j in length(mtt):1) {     # loop through variations in 3'->5' direction
      
      if(start(mtt)[j]==width(ref_seq)){
        subseq(seqs[[i]], start=start(mtt)[j], 
               end=start(mtt)[j]) <- to_insert[[i]][j]
      } else {
        subseq(seqs[[i]], start=start(mtt)[j], 
             end=start(mtt)[j]+1) <- to_insert[[i]][j]
      }
    }     
    # remove '-'
    seqs[[i]] <- DNAStringSet(gsub("-","",seqs[[i]]))

    names(seqs[[i]]) <- paste0(names(seqs[[i]]), ".", letters[i])
  }
  seqs <- unlist(DNAStringSetList(seqs))
  
  if (as.character(strand(x_ex)[1]) == "-")
    seqs <- reverseComplement(seqs)
  
  seqs
}

get_alleles <- function(df){
  z <- matrix(unlist(df), ncol=length(df))
  alleles <- as.data.frame(z) %>% group_by_all %>% count
  alleles <- arrange(alleles)
  
  alleles <- alleles[seq_len(min(2,dim(alleles)[1])),]
  
  strspl = data.table::transpose(as.list(alleles[-dim(alleles)[2]]))
}

generate_variant_transcripts <- function(v, x, 
                                         gene, 
                                         bam_file = "aln_s.bam", verbose = TRUE) {
  
  if(verbose)
    message(paste0("working on '", gene, "'"))
  
  # extract all exons for given gene  
  x_ex <- extract_exons(x, gene) ##
  
  # find overlaps b/w variation and this gene
  fop <- findOverlapPairs(rowRanges(v), x_ex) 

  #check that there IS an overlap
  if (!isEmpty(pintersect(fop))) {
  
    # read alignments for this region
    gr <- range(x_ex)
    sbp <- ScanBamParam(which = gr, what = scanBamWhat())
    bamWhat(sbp) <- "seq"
    ga <- readGAlignments(bam_file, 
                          param = sbp, use.names = TRUE)
    
    # full set of variants for this region
    variation_regs <- unique(S4Vectors::first(fop))
    
    # ranges for each read
    rngs <- get_read_ranges(ga) ##
    
    # split into transcripts
    x_exs <- split(x_ex, x_ex$transcript_id)
    
    # intersect all reads w/ all transcripts, form matrix
    m <- matrix(0,
                ncol = length(variation_regs),
                nrow = length(rngs),
                dimnames = list(names(rngs), names(variation_regs))
    )
    rngs_u <- unlist(rngs)
    fo <- findOverlaps(rngs_u, ranges(variation_regs))
    
    # which variations in a certain alignment read, and 
    # below notate with 0/1 if absent/present
    pop <- split(subjectHits(fo),
                 names(rngs_u)[queryHits(fo)]) 
    for (i in 1:length(pop)) {
      m[names(pop)[i], unique(pop[[i]])] <- 1
    }
    
    transcripts <- unique(S4Vectors::second(fop)$transcript_id)
    new_seqs <- vector("list", length(transcripts))
    
    # remove reads that end right at the variation
    ga <- ga[ !(end(ga) %in% start(variation_regs)) ]
    
    # narrow all alignments around variation
    df_tot <- narrowal(variation_regs,ga)

    # loop through affected transcripts    
    for (i in seq_len(length(transcripts))) {
      tr <- transcripts[i]
      if(verbose)
        message(tr)

      f <- S4Vectors::first(fop)[S4Vectors::second(fop)$transcript_id == tr]
      nm <- names(f)
      m1 <- as.data.frame(m)[nm]

      # only take reads that overlap all variation pos. of this transcript
      w <- rowSums(m1[, nm, drop = FALSE]) == length(nm)
      rds <- rownames(m1)[w]
      m2 <- subset(m1, rownames(m1) %in% rds)

      # extract transcript seq from genome + 
      # positions where we must insert variation (transcript coord.)
      dss_ref <- getSeq(Hsapiens, x_ex[x_ex$transcript_id == tr])
      ref_seq <- DNAStringSet(paste(as.character(dss_ref),
                                    collapse = ""))
      names(ref_seq) <- tr
      mtt <- mapToTranscripts(variation_regs[nm], x_exs[tr])

      df <- df_tot[nm]
      
      df1 <- lapply(df, function(u) {
        keep <- rownames(u) %in% rownames(m2)
        as.character(u$seq[keep])
      })
      
      # get alleles and insert them into reference sequence
      to_insert <- get_alleles(df1)
      new_seqs[[i]] <- seqs_with_var(to_insert, ref_seq, 
                                     mtt, x_ex, tr)
  
    }
  }
  

  transcripts_unmod <- setdiff(unique(x_ex$transcript_id), transcripts)
  seqs_unmod <- vector("list", length(transcripts_unmod))

  # check for transcripts with NO overlap
  if(!isEmpty(transcripts_unmod)) {

    for (j in seq_len(length(transcripts_unmod))) {
      tr_unmod <- transcripts_unmod[j]
      if(verbose)
        message(tr_unmod)

      dss_ref_unmod <- getSeq(Hsapiens, x_ex[x_ex$transcript_id == tr_unmod])
      ref_seq_unmod <- DNAStringSet(paste(as.character(dss_ref_unmod),
                                        collapse = ""))
      names(ref_seq_unmod) <- tr_unmod

      if(as.character(strand(x_ex[x_ex$transcript_id == tr_unmod])[1])=='-'){
        ref_seq_unmod <- reverseComplement(ref_seq_unmod)
      }

      seqs_unmod[[j]] <- ref_seq_unmod
    }
  }
  
  append(new_seqs, seqs_unmod)

}
