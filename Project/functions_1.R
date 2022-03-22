# get_read_ranges <- function(ga) {
#   nm_new <- names(ga)
#   rngs_new <- cigarRangesAlongReferenceSpace(cigar(ga),
#                                          pos = start(ga),
#                                          drop.empty.ranges = FALSE,
#                                          N.regions.removed = FALSE,
#                                          with.ops = TRUE,
#                                          reduce.ranges = FALSE) %>%
#     setNames(nm_new)
# 
#   nms_new <- rep(nm_new, lengths(rngs_new))
#   rngs_u_new <- unlist(rngs_new, use.names = FALSE)
# 
#   keep_new <- names(rngs_u_new) != "N"
#   rngs_new <- split(rngs_u_new[keep_new], nms_new[keep_new])[nm_new]
# 
#   dokimh2  <- lapply(rngs_new, function(u) {
#     nm_new <- names(u)
#     names(u) <- NULL
#     u[]
#   })
# 
#   dokimh2
# 
# }
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


extract_exons <- function(x, gene = NULL) {
  x_ex <- x[x$type=="exon" & seqnames(x) != "chrM"]
  if(!is.null(gene))
    x_ex <- x_ex[x_ex$gene_id == gene]
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
    mcols(na)
  })
  
}


seqs_with_var <- function(to_insert, ref_seq, mtt, ref_base, alt_base){
  to_insert <- to_insert[lengths(to_insert)!=0]
  seqs <- vector("list", length(to_insert))
  
  start <- vector("list", length(mtt))
  
  for(i in 1:length(to_insert)) { # loop through alleles
    seqs[[i]] <- ref_seq
    for(j in length(mtt):1) {     # loop through variations in 3'->5' direction
      
      #start[j] = start(mtt)[j]
      st <- start(mtt)[j]
      en <- st+1
      
      
      # Check if the variation is on the very last base. Deletion can deal with it normally but the other types of variation are concerning.
      if(start(mtt)[j]== width(ref_seq) && !nchar(ref_base[[j]])>nchar(alt_base[[j]])){
        # subseq(seqs[[i]], start=start[[j]], 
        #        end=start[[j]]) <- to_insert[[i]][j]
        en <- st
      }
        
      subseq(seqs[[i]], start=st, 
               end=en) <- to_insert[[i]][j]
      
    } 
    names(seqs[[i]]) <- paste0(names(seqs[[i]]), ".", letters[i])  
  }  
  
  seqs <- unlist(DNAStringSetList(seqs))
  
}



get_alleles <- function(df){
  z <- matrix(unlist(df), ncol=length(df), 
              dimnames = list(NULL,names(df)))
  alleles <- as.data.frame(z) %>% group_by_all %>% count %>%
    arrange(desc(n)) 
  # take top 2 alleles, transpose, return list
  alleles[seq_len(min(2,nrow(alleles))),] %>%
    select(-n) %>% as.data.frame %>% 
    t %>% as.data.frame %>% as.list
}

generate_variant_transcripts <- function(v, x,
                                         gene, 
                                         bam_file = "aln_s.bam", verbose = FALSE) {
  
  if(verbose)
    message(paste0("working on '", gene, "'"))
  
  # extract exons for given gene  
  x_ex <- x[x$gene_id==gene] ##
  # split into transcripts
  x_exs <- split(x_ex, x_ex$transcript_id)
  
  # find overlaps b/w variation and this gene
  fop <- findOverlapPairs(rowRanges(v), x_ex) 
  transcripts <- unique(S4Vectors::second(fop)$transcript_id)
  
  new_seqs <- vector("list", length(transcripts))
  
  #check that there IS an overlap
  if (!isEmpty(pintersect(fop))) {
    
    # read alignments for this region
    gr <- range(x_ex)
    sbp <- ScanBamParam(which = gr, what = scanBamWhat())
    bamWhat(sbp) <- "seq"
    ga <- readGAlignments(bam_file, 
                          param = sbp, use.names = TRUE)
    
    if(verbose)
      message(paste0(length(ga), " total reads."))
    
    # full set of variants for this region
    variation_regs <- unique(S4Vectors::first(fop))
    
    # ranges for each read
    # rngs_old <- get_read_ranges_old(ga) ##
    rngs <- get_read_ranges(ga) ##
    
    
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
    ##m_1 OK
    # remove reads that never overlap
    keep <- rowSums(m)>0
    m <- m[keep,,drop=FALSE]
    #m_2 NOT OK
    ga <- ga[keep]
    
    # remove reads that end right at the variation
    ga <- ga[ !(end(ga) %in% start(variation_regs)) ]
    
    # narrow all alignments around variation
    df_tot <- narrowal(variation_regs,ga)
    
    ###########
    

    for (p in 1:length(df_tot)) {
      
      if(width(variation_regs$REF[p])>length(variation_regs$ALT[p])){ # Deletions
        
        for(q in 1:length(df_tot[p][[1]]$seq)){
          
          if(df_tot[p][[1]]$seq[q]!=variation_regs$REF[p] && df_tot[p][[1]]$seq[q]!=variation_regs$REF[p] ){
            
            read_name <- rownames(df_tot[p][[1]])
            m[read_name,p]<-0
            
          }
          
        }
        
      }else{ # SNPs, Insertions
        
        for(q in 1:length(df_tot[p][[1]]$seq)){
          
          if(substring(df_tot[p][[1]]$seq[q],1,width(df_tot[p][[1]]$seq[q])-1)!=variation_regs$REF[p][[1]] && substring(df_tot[p][[1]]$seq[q],1,width(df_tot[p][[1]]$seq[q])-1)!=variation_regs$ALT[p][[1]]){
            
            read_name <- rownames(df_tot[p][[1]])
            m[read_name,p]<-0
            
          }
          
        }
        
      }
      
    }
    
    
    
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
      tr_gr <- x_ex[x_ex$transcript_id == tr]
      # 
      # gr <- GRanges(seqnames = seqnames(tr_gr),
      #                 IRanges(ranges(tr_gr)), strand = "+")
      # 
      # 
      # 
      # 
      # names(ref_seq) <- tr
      # 
      # mtt <- mapToTranscripts(variation_regs[nm], x_exs[tr])
      # 
      
      strand(tr_gr) <- "+"
      grl <- GRangesList(tr_gr)
      names(grl) <- tr
      mtt <- mapToTranscripts(variation_regs[nm], grl)
      dss_ref <- getSeq(Hsapiens, grl,
                        as.character = TRUE)
      ref_seq <- DNAStringSet(paste(dss_ref, collapse = ""))
      
      df <- df_tot[nm]
      
      df1 <- lapply(df, function(u) {
        keep <- rownames(u) %in% rownames(m2)
        as.character(u$seq[keep])
      })
      
      ref_base <- variation_regs[nm]$REF
      alt_base <- variation_regs[nm]$ALT
      
      #get alleles and insert them into reference sequence
      
      to_insert <- get_alleles(df1) 
      new_seqs[[i]] <- seqs_with_var(to_insert, ref_seq, mtt,ref_base, alt_base)
      
      
      # Recording change in transcript
      new_info_line <- paste(c(tr, paste(nm, collapse=";"), paste(lapply(as.data.frame(variation_regs[nm])[,2], function(u){paste(c(u,u+1), collapse="-")}),collapse=";"), paste(lapply(as.array(start(mtt)), function(u){paste(c(u,u+1), collapse="-")}),collapse=";"), paste(lapply(to_insert, function(u){paste(as.array(u),collapse=";")}), collapse = "    ")), collapse="    ")
      
      write(new_info_line, file = "my_file.txt", append = TRUE)
      
    }
  }
  
  
  transcripts_unmod <- setdiff(names(x_exs), transcripts)
  seqs_unmod <- vector("list", length(transcripts_unmod))
  
  # check for transcripts with NO overlap
  
  if(!isEmpty(transcripts_unmod)) {
    
    for (j in seq_len(length(transcripts_unmod))) {
      tr_unmod <- transcripts_unmod[j]
      
      tr_gr_unmod <- x_ex[x_ex$transcript_id == tr_unmod]
      
      if(verbose)
        message(tr_unmod)
      
      
      strand(tr_gr_unmod) <- "+"
      grl_unmod <- GRangesList(tr_gr_unmod)
      names(grl_unmod) <- tr_unmod
      dss_ref_unmod <- getSeq(Hsapiens, grl_unmod,
                              as.character = TRUE)
      ref_seq_unmod <- DNAStringSet(paste(dss_ref_unmod, collapse = ""))

      names(ref_seq_unmod) <- tr_unmod

      seqs_unmod[[j]] <- ref_seq_unmod
    }
  }
  
  
  
  if(as.character(strand(x_ex))[1]=='-'){
    seqs_unmod <- lapply(seqs_unmod, function(u) {
      rev_unmod<-reverseComplement(u)
    }) ############
    new_seqs <- lapply(new_seqs, function(u) {
      rev_unmod<-reverseComplement(u)
    })
  }
  
  
  new_seqs <- append(new_seqs, seqs_unmod)
  unlist(DNAStringSetList(new_seqs))  
}





test_variant <- function(new_seqs, gene, which_transcript,ex_seq,v = v,x =x,
                         bam_file = alns, verbose = FALSE ){
  x_ex <- x[x$gene_id==gene]
  x_exs <- split(x_ex, x_ex$transcript_id)
  
  new_seqs <- generate_variant_transcripts(v = v,x = x,
                                           bam_file = alns, gene = gene, verbose = FALSE)
  
  
  a <- new_seqs[[paste0(which_transcript,".a")]]
  b <- new_seqs[[paste0(which_transcript,".b")]]
  
  tr_gr <- x_ex[x_ex$transcript_id == which_transcript]
  gr <- GRanges(seqnames = seqnames(tr_gr),
                IRanges(ranges(tr_gr)), strand = "+")
  
  dss_ref <- getSeq(Hsapiens, gr,
                    as.character = TRUE)
  ref_seq <- DNAStringSet(paste(dss_ref, collapse = ""))
  
  alleles_of_interest<-DNAStringSet(c(as.character(a), as.character(b)))
  loc1<-str_locate(as.character(ref_seq),ex_seq)[1]
  loc2<-str_locate(as.character(ref_seq),ex_seq)[2]
  
  if(as.character(strand(x_ex))[1]=='-'){
    ret_substring <- substring(as.character(reverseComplement(alleles_of_interest)),loc1,loc2)
  } else {
    ret_substring <- substring(as.character((alleles_of_interest)),loc1,loc2)
  }
  ret_substring
}

test_NO_variant <- function(new_seqs, gene, which_transcript,ex_seq,v = v,x =x,
                         bam_file = alns, verbose = FALSE ){
  x_ex <- x[x$gene_id==gene]
  x_exs <- split(x_ex, x_ex$transcript_id)
  
  new_seqs <- generate_variant_transcripts(v = v,x = x,
                                           bam_file = alns, gene = gene, verbose = FALSE)
  
  
  a <- new_seqs[[which_transcript]]
  tr_gr <- x_ex[x_ex$transcript_id == which_transcript]
  gr <- GRanges(seqnames = seqnames(tr_gr),
                IRanges(ranges(tr_gr)), strand = "+")
  
  dss_ref <- getSeq(Hsapiens, gr,
                    as.character = TRUE)
  ref_seq <- DNAStringSet(paste(dss_ref, collapse = ""))
  
  alleles_of_interest<-DNAStringSet(c(as.character(a)))
  loc1<-str_locate(as.character(ref_seq),ex_seq)[1]
  loc2<-str_locate(as.character(ref_seq),ex_seq)[2]
  
  if(as.character(strand(x_ex))[1]=='-'){
    ret_substring <- substring(as.character(reverseComplement(alleles_of_interest)),loc1,loc2)
  } else {
    ret_substring <- substring(as.character((alleles_of_interest)),loc1,loc2)
  }
  ret_substring
}




