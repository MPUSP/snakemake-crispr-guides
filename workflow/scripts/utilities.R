# add off target candidates with biostrings
addOffTargetSequences <- function(
    guideset, guide_aligner, max_cores, txdb, genome_dna, genome_index) {
  if (guide_aligner == "biostrings") {
    guideset_chunks <- split(
      guideset,
      ceiling(
        seq_along(guideset) /
          (length(guideset) /
            max_cores)
      )
    )

    guideset_chunks <- mclapply(
      mc.cores = max_cores,
      X = guideset_chunks,
      FUN = function(guide_chunk) {
        addSpacerAlignments(
          guide_chunk,
          aligner = "biostrings",
          txObject = txdb,
          custom_seq = genome_dna,
          n_mismatches = 3,
          n_max_alignments = 1000,
          addSummary = TRUE,
          all_alignments = TRUE
        )
      }
    )
    # merge all chunks to single guideSet again
    guideset <- Reduce("c", guideset_chunks)
  } else if (guide_aligner == "bowtie") {
    guideset <- addSpacerAlignments(
      guideset,
      aligner = "bowtie",
      aligner_index = paste0(genome_index, "/index"),
      bsgenome = bsgenome,
      addSummary = TRUE,
      n_mismatches = 3,
      custom_seq = genome_dna,
    )
  }
  return(guideset)
}


# add ON target scores based on spacer and protospacer sequence
addAllOnTargetScores <- function(
    guideset, score_methods, spacer_length, tss_window = NULL) {
  guideset <- addOnTargetScores(
    guideset,
    methods = setdiff(score_methods, c("tssdist", "genrich"))
  )

  # rescale scores to a range between 0 and 1 if they exceed this range
  rescale_score <- function(x) {
    (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  }
  for (score in setdiff(score_methods, c("tssdist", "genrich"))) {
    score_val <- guideset@elementMetadata[[paste0("score_", score)]]
    if (min(score_val, na.rm = TRUE) < 0 || max(score_val, na.rm = TRUE) > 1) {
      guideset@elementMetadata[[paste0("score_", score)]] <- score_val %>%
        rescale_score()
    }
  }

  # add score based on G enrichment in seed region (pos -4 to -14 from PAM),
  # see Miao & Jahn et al., The Plant Cell, 2023
  if ("genrich" %in% score_methods) {
    guideset$score_genrich <- guideset$protospacer %>%
      subseq(start = spacer_length - 13, end = spacer_length - 3) %>%
      letterFrequency("G", as.prob = TRUE) %>%
      as.numeric() %>%
      rescale_score()
  }

  # add score based on distance to TSS; lower dist = higher score
  if ("tssdist" %in% score_methods) {
    if (is.null(tss_window)) {
      guideset$score_tssdist <- 0
    } else {
      guideset$score_tssdist <- 1 - (abs(guideset$dist_to_tss) / max(abs(tss_window)))
    }
  }

  return(guideset)
}

# add information on possible bad seed (seed region = 12 bp upstream of PAM)
addBadSeeds <- function(guideset, bad_seeds) {
  for (seed_pattern in bad_seeds) {
    guideset@elementMetadata[[paste0("seed_", seed_pattern)]] <- vmatchPattern(
      pattern = seed_pattern,
      subject = subseq(guideset$protospacer, start = spacer_length - 11, end = spacer_length)
    ) %>%
      sapply(FUN = function(x) length(width(x))) %>%
      as.logical()
  }
  return(guideset)
}

# reduce overlap: from each set of overlapping guides,
# select the worst one by lowest mean score and remove it
filter_overlaps <- function(guideset, guide_length = spacer_length) {
  index_overlap <- findOverlaps(guideset,
    maxgap = spacer_length,
    drop.self = TRUE, drop.redundant = TRUE
  ) %>%
    as.list() %>%
    lapply(function(x) {
      x[which.min(guideset$score_all[x])]
    }) %>%
    unlist() %>%
    {
      .[!duplicated(.)]
    }
  index_overlap
}

# function to export guide tables as GFF file
export_gff <- function(
    df, filename, feature = "sgRNA", source = "snakemake-crispr-guides") {
  gr <- makeGRangesFromDataFrame(df,
    seqnames.field = "seqnames",
    keep.extra.columns = TRUE
  )
  gr$type <- feature
  gr$ID <- gr$group_name
  gr$Parent <- NA
  rtracklayer::export(
    gr, filename,
    version = "3",
    source = source
  )
}
