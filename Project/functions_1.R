get_read_ranges <- function(ga) {
  nm <- names(ga)
  rngs <- cigarRangesAlongReferenceSpace(cigar(ga), 
                                         pos = start(ga),
                                         drop.empty.ranges = FALSE,
                                         N.regions.removed = FALSE,
                                         with.ops = TRUE,
                                         reduce.ranges = FALSE) %>% 
    setNames(nm)
  
  nms <- rep(nm, lengths(rngs))
  rngs_u <- unlist(rngs, use.names = FALSE)
  
  keep <- names(rngs_u) != "N"
  split(unname(rngs_u[keep]), nms[keep])[nm]
}


extract_exons <- function(x, gene = NULL) {
  
  x_ex <- x[x$type == "exon" & seqnames(x) != "chrM"]
  
  if (!is.null(gene))
    x_ex <- x_ex[x_ex$gene_id == gene]
  
  w <- which(seqlevels(x_ex) == "chrM")
  seqlevels(x_ex) <- seqlevels(x_ex)[-w]
  x_ex
}


narrowal <- function(variation, ga) {
  
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


seqs_with_var <- function(to_insert, ref_seq, mtt, ref_base, alt_base) {
  
  to_insert <- to_insert[lengths(to_insert) != 0]
  seqs <- vector("list", length(to_insert))
  
  start <- vector("list", length(mtt))
  
  for (i in 1:length(to_insert)) { # loop through alleles
    seqs[[i]] <- ref_seq
    for (j in length(mtt):1) {    # loop through variations in 3'->5' direction
      st <- start(mtt)[j]
      
      if (nchar(ref_base[[j]])>2) {
        en <- st + nchar(ref_base[[j]])-1
      } else {
        en <- st+1
        }
      
      # Check if the variation is on the very last base. 
      #Deletion can deal with it normally but the other types cannot.
      if (start(mtt)[j] ==  width(ref_seq) && 
          !nchar(ref_base[[j]])>nchar(alt_base[[j]])) {
        en <- st
      }
      
      subseq(seqs[[i]], start = st, 
               end = en) <- to_insert[[i]][j]
    } 
    names(seqs[[i]]) <- paste0(names(seqs[[i]]), ".", letters[i])  
  }  
  seqs <- unlist(DNAStringSetList(seqs))
}


get_alleles <- function(df) {
  
  z <- matrix(unlist(df), ncol = length(df), 
              dimnames = list(NULL,names(df)))
  alleles <- as.data.frame(z) %>% group_by_all %>% count %>%
    arrange(desc(n)) 
  # take top 2 alleles, transpose, return list
  alle <- alleles[seq_len(min(2,nrow(alleles))), ] %>%
    select(-n) %>% as.data.frame %>% 
    t %>% as.data.frame %>% as.list
  alle2_length <- vector()
  alle1_length <- length(which(apply(z, 1, function(x) 
    identical(as.character(x), alle[[1]]))))
  
  if (length(alle) == 2) {
    alle2_length <- length(which(apply(z, 1, function(x) 
      identical(as.character(x), alle[[2]]))))}
  alle_lengths <- c(alle1_length, alle2_length)
  
  list(alle, alle_lengths)
}


generate_variant_transcripts <- function(v, x, gene, bam_file = "aln_s.bam",
                                         verbose = FALSE) {
  
  if (verbose)
    message(paste0("working on '", gene, "'"))
  
  # extract exons for given gene  
  x_ex <- x[x$gene_id == gene] 
  # split into transcripts
  x_exs <- split(x_ex, x_ex$transcript_id)
  
  # find overlaps b/w variation and this gene
  fop <- findOverlapPairs(rowRanges(v), x_ex) 
  transcripts <- unique(S4Vectors::second(fop)$transcript_id)
  
  new_seqs <- vector("list", length(transcripts))

  for_next_loop <- vector()
  
  #check that there IS an overlap
  if (!isEmpty(pintersect(fop))) {
    # read alignments for this region
    gr <- range(x_ex)
    sbp <- ScanBamParam(which = gr, what = scanBamWhat())
    bamWhat(sbp) <- "seq"
    ga <- readGAlignments(bam_file, 
                          param = sbp, use.names = TRUE)
    
    if (verbose)
      message(paste0(length(ga), " total reads."))
    
    # full set of variants for this region
    variation_regs <- unique(S4Vectors::first(fop))
    
    # ranges for each read
    rngs <- get_read_ranges(ga)
    
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

    # remove reads that never overlap
    keep <- rowSums(m)>0
    m <- m[keep,,drop = FALSE]
    ga <- ga[keep]
    
    # remove reads that end right at the variation
    ga <- ga[!(end(ga) %in% start(variation_regs))]
    
    # narrow all alignments around variation
    df_tot <- narrowal(variation_regs,ga)
    
    #apoth<-m

    # loop through affected transcripts    
    for (i in seq_len(length(transcripts))) {
      tr <- transcripts[i]
      if (verbose)
        message(tr)
      
      if (all(m == 0)) {
        for_next_loop <- c(for_next_loop,tr)
        next
      }
      
      f <- S4Vectors::first(fop)[S4Vectors::second(fop)$transcript_id == tr]
      nm <- names(f)
      m1 <- as.data.frame(m)[nm]
      
      if (all(m1 == 0)) {
        for_next_loop <- c(for_next_loop,tr)
        next
      }
      
      # only take reads that overlap all variation pos. of this transcript
      w <- rowSums(m1[ , nm, drop = FALSE]) == length(nm)
      rds <- rownames(m1)[w]
      m2 <- subset(m1, rownames(m1) %in% rds)
      
      if (all(sum(isEmpty(m2)) == length(m1))) {
        for_next_loop <- c(for_next_loop,tr)
        next
      }
      
      
      variation_spec <- variation_regs[names(variation_regs)%in%colnames(m2)]
      
      df_tot_spec <- df_tot[names(df_tot)%in%colnames(m2)]
      df_tot_specs <- lapply(df_tot_spec, function(u){
        as.data.frame(u)[rownames(as.data.frame(u)) %in% rownames(m2), ]})
      
      #apoth2<-m2

      for (p in 1:length(m2)) {

          if (width(variation_spec$REF[p])>
              length(variation_spec$ALT[p])) { # Deletions

            read_nuc <- df_tot_specs[[p]]

          } else {

            read_nuc <- substring(df_tot_specs[[p]],
                                  1,nchar(df_tot_specs[[p]])-1)
          }

        w <- read_nuc != as.character(variation_spec$REF[p]) &
        read_nuc != as.character(variation_spec$ALT[[p]])
        
        m2[w,p] <- 0
      }
    
      w2 <- rowSums(m2[ , nm, drop = FALSE]) == length(m2)
      rds2 <- rownames(m2)[w2]
      m2 <- subset(m2, rownames(m2) %in% rds2)
     
      tr_gr <- x_ex[x_ex$transcript_id == tr]
      
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
      
      allele_info <- get_alleles(df1) 
      to_insert <- allele_info[1][[1]]
      new_seqs[[i]] <- seqs_with_var(to_insert, ref_seq, mtt,ref_base, alt_base)
      
      # Recording change in transcript
      ref_alt <- paste(nm, collapse = ";")
      
      repl_range <- vector()
      
      for (t in 1:length(ref_base)) {
        if (width(ref_base[t]) == 1 || width(ref_base[t]) == 2) {
          repl_range <- c(repl_range, 1)
        } else {
          repl_range <- c(repl_range, width(ref_base[t])-1)
        }
      }
      
      gen_loc_list <- list(as.data.frame(variation_regs[nm])[ ,2], repl_range)
      tr_loc_list <- list(as.array(start(mtt)), repl_range)
      
      loc_fun <- function(el1, el2){
        paste(c(el1,el1+el2), collapse = "-")
      }
    
     gen_loc <- paste(mapply(loc_fun,gen_loc_list[[1]],gen_loc_list[[2]]),
                   collapse = ";")
      
     tr_loc <- paste(mapply(loc_fun,tr_loc_list[[1]],tr_loc_list[[2]]),
                     collapse = ";")
      
     alleles_obt <- paste(lapply(to_insert, 
                          function(u){paste(as.array(u),collapse = ";")}), 
                   collapse = "    ")
      
     new_info_line <- paste(c(tr, ref_alt, gen_loc, tr_loc, alleles_obt, 
                               allele_info[[2]]),  
                   collapse = "    ")
     
     write(new_info_line, file = "my_file.txt", append = TRUE)
      
    }
  }
  
  
  transcripts_unmod <- setdiff(names(x_exs), transcripts)
  seqs_unmod <- vector("list", length(transcripts_unmod))
  
  transcripts_unmod <- c(for_next_loop, transcripts_unmod)
  # check for transcripts with NO overlap
  
  if (!isEmpty(transcripts_unmod)) {
    
    for (j in seq_len(length(transcripts_unmod))) {
      tr_unmod <- transcripts_unmod[j]
      tr_gr_unmod <- x_ex[x_ex$transcript_id == tr_unmod]
      
      if (verbose)
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
  
  
  
  if (as.character(strand(x_ex))[1] == '-') {
    seqs_unmod <- lapply(seqs_unmod, function(u) {
      rev_unmod <- reverseComplement(u)
    }) 
    new_seqs <- lapply(new_seqs, function(u) {
      rev_unmod <- reverseComplement(u)
    })
  }
  
  new_seqs <- append(new_seqs, seqs_unmod)
  unlist(DNAStringSetList(new_seqs))  
}


test_variant <- function(new_seqs, gene, which_transcript,ex_seq,v = v,x = x,
                         bam_file = alns, verbose = FALSE) {
  
  x_ex <- x[x$gene_id == gene]
  x_exs <- split(x_ex, x_ex$transcript_id)
  
  new_seqs <- generate_variant_transcripts(v = v,x = x, bam_file = alns, 
                                           gene = gene, verbose = FALSE)
  
  
  a <- new_seqs[[paste0(which_transcript,".a")]] # allele a
  b <- new_seqs[[paste0(which_transcript,".b")]] # allele b
  
  tr_gr <- x_ex[x_ex$transcript_id == which_transcript]
  gr <- GRanges(seqnames = seqnames(tr_gr),
                IRanges(ranges(tr_gr)), strand = "+")
  
  dss_ref <- getSeq(Hsapiens, gr,
                    as.character = TRUE)
  ref_seq <- DNAStringSet(paste(dss_ref, collapse = ""))
  
  alleles_of_interest <- DNAStringSet(c(as.character(a), as.character(b)))
  loc1 <- str_locate(as.character(ref_seq),ex_seq)[1]
  loc2 <- str_locate(as.character(ref_seq),ex_seq)[2]
  
  if (as.character(strand(x_ex))[1] == '-') {
    rev_alleles_of_interest <- reverseComplement(alleles_of_interest)
    ret_substring <- substring(as.character(rev_alleles_of_interest),loc1,loc2)
  } else {
    ret_substring <- substring(as.character((alleles_of_interest)),loc1,loc2)
  }
  ret_substring
}


test_NO_variant <- function(new_seqs, gene, which_transcript,ex_seq,v = v,x = x,
                         bam_file = alns, verbose = FALSE) {
  
  x_ex <- x[x$gene_id == gene]
  x_exs <- split(x_ex, x_ex$transcript_id)
  
  new_seqs <- generate_variant_transcripts(v = v,x = x, bam_file = alns, 
                                           gene = gene, verbose = FALSE)
  
  
  a <- new_seqs[[which_transcript]] # allele a
  tr_gr <- x_ex[x_ex$transcript_id == which_transcript]
  gr <- GRanges(seqnames = seqnames(tr_gr),
                IRanges(ranges(tr_gr)), strand = "+")
  
  dss_ref <- getSeq(Hsapiens, gr,
                    as.character = TRUE)
  ref_seq <- DNAStringSet(paste(dss_ref, collapse = ""))
  
  alleles_of_interest <- DNAStringSet(c(as.character(a)))
  loc1 <- str_locate(as.character(ref_seq),ex_seq)[1]
  loc2 <- str_locate(as.character(ref_seq),ex_seq)[2]
  
  if (as.character(strand(x_ex))[1] == '-') {
    rev_alleles_of_interest <- reverseComplement(alleles_of_interest)
    ret_substring <- substring(as.character(rev_alleles_of_interest),loc1,loc2)
  } else {
    ret_substring <- substring(as.character((alleles_of_interest)),loc1,loc2)
  }
  ret_substring
}




