# LOAD PACKAGES
# ------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(Biostrings)
  library(GenomeInfoDbData)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(crisprBase)
  library(crisprDesign)
  library(crisprBwa)
  library(crisprScore)
  library(parallel)
})


# CONFIGURATION
# ------------------------------
sm_params <- snakemake@params[[1]]
max_cores <- snakemake@params[[2]]
for (param in names(sm_params)) {
  assign(param, sm_params[[param]], envir = .GlobalEnv)
}
output_top <- snakemake@output[["guideRNAs_top"]]
output_fail <- snakemake@output[["guideRNAs_fail"]]
output_ntc <- snakemake@output[["guideRNAs_ntc"]]


# STAGE 1 : FILE PREPARATION
# ------------------------------
# import genome sequence
messages <- c("importing genome sequence and annotation")
genome_fasta <- snakemake@input[["fasta"]]
genome_gff <- snakemake@input[["gff"]]
genome_index <- snakemake@input[["bowtie_index"]]
genome_dna <- Biostrings::readDNAStringSet(genome_fasta)
load(snakemake@input[["seqinfo"]])
seqinfo(genome_dna) <- seqinfo_genome
load(snakemake@input[["bsgenome"]])

# import genome annotation with chromosome metadata
txdb <- makeTxDbFromGFF(
  file = genome_gff,
  organism = unname(genome(seqinfo_genome)[1]),
  chrominfo = seqinfo_genome
)

# optionally select only specific chromosome/region
if (!is.null(target_region)) {
  genome_dna <- genome_dna[names(genome_dna) %in% target_region]
}


# STAGE 2 : FINDING GUIDES
# ------------------------------
# predict all possible guides for target regions
messages <- append(
  messages,
  paste0("predicting guide RNAs for enzyme: ", crispr_enzyme)
)
data(list = c(crispr_enzyme), package = "crisprBase")
list_pred_guides <- findSpacers(
  x = genome_dna,
  canonical = canonical,
  both_strands = TRUE,
  spacer_len = spacer_length,
  crisprNuclease = get(crispr_enzyme)
)

# extract (pseudo) transcription start sites (TSS)
list_tx <- transcripts(txdb)
list_tss <- resize(list_tx, width = 1)
list_tss$ID <- list_tss$tx_name

# map guides to TSS windows
list_pred_guides <- addTssAnnotation(
  list_pred_guides,
  list_tss,
  tss_window = tss_window
)

# remove all guides that do not target within a TSS window
index_non_target <- unlist(mclapply(
  mc.cores = max_cores,
  X = list_pred_guides$tssAnnotation,
  FUN = nrow
))
list_pred_guides <- list_pred_guides[index_non_target != 0]

# clip TSS windows where the 5'-UTR extends into another transcript
# then remove all guides that are beyond the 3'- or 5' end of the target gene
list_intergenic <- {start(list_tx) - c(1, end(list_tx[-length(list_tx)]))} %>%
  c(., unname(tail(seqlengths(txdb), 1) - tail(end(list_tx), 1))) %>%
  replace(., . > abs(tss_window[1]), abs(tss_window[1])) %>%
  replace(., . < tss_window[1] / 2, tss_window[1] / 2)
max_5UTR_index <- seq_along(list_tx)
max_5UTR_index[as.logical(strand(list_tx) == "-")] <-
  max_5UTR_index[as.logical(strand(list_tx) == "-")] + 1
max_5UTR <- list_intergenic[max_5UTR_index]
list_tx_window <- flank(list_tx, max_5UTR)
list_tx_window <- punion(list_tx, list_tx_window)
list_tss_window <- punion(
  flank(list_tss, abs(tss_window[1])),
  resize(list_tss, width = tss_window[2])
)

# select only guides whose PAM overlaps with the TSS window
list_tss_window <- pintersect(list_tss_window, list_tx_window)
genome(list_tss_window) <- "custom"
list_targets <- as.list(findOverlaps(
  list_pred_guides, list_tss_window,
  ignore.strand = TRUE
))
names(list_targets) <- names(list_pred_guides)
list_targets_keep <- unlist(lapply(list_targets, length) == 1)
list_pred_guides <- list_pred_guides[list_targets_keep]
list_targets <- list_targets[list_targets_keep]

# trim original TSS annotation to one entry per guide
df_tss_annot <- list_pred_guides$tssAnnotation %>%
  as_tibble() %>%
  right_join(
    by = c("group_name", "tx_id"),
    enframe(
      unlist(list_targets),
      name = "group_name",
      value = "tx_id"
    )
  ) %>%
  dplyr::select(-group, -chr, -strand) %>%
  mutate(tx_width = setNames(width(list_tx), list_tx$tx_id)[tx_id]) %>%
  arrange(as.numeric(str_sub(group_name, 8L)))

# add modified TSS information to metadata
list_pred_guides@elementMetadata <- cbind(
  list_pred_guides@elementMetadata,
  df_tss_annot
)

# remove guides where a transcript was mapped but not the right one
list_pred_guides <- list_pred_guides[!is.na(list_pred_guides$tss_id)]

# remove guides targeting RNAs (optional)
if (filter_rna) {
  list_rna_targets <- grepl("^rna-", list_pred_guides$tx_name)
  list_pred_guides <- list_pred_guides[!list_rna_targets]
}

# add spacer and PAM sequence features
list_pred_guides <- addSequenceFeatures(list_pred_guides)

# add off target candidates with biostrings
if (guide_aligner == "biostrings") {
  list_pred_guides_chunks <- split(
    list_pred_guides,
    ceiling(
      seq_along(list_pred_guides) /
        (length(list_pred_guides) /
          max_cores)
    )
  )

  list_pred_guides_chunks <- mclapply(
    mc.cores = max_cores,
    X = list_pred_guides_chunks,
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
  list_pred_guides <- Reduce("c", list_pred_guides_chunks)
} else if (guide_aligner == "bowtie") {
  list_pred_guides <- addSpacerAlignments(
    list_pred_guides,
    aligner = "bowtie",
    aligner_index = paste0(genome_index, "/index"),
    bsgenome = bsgenome,
    addSummary = TRUE,
    n_mismatches = 3,
    custom_seq = genome_dna,
  )
}

# add OFF target scores based on alignment
list_pred_guides <- addOffTargetScores(list_pred_guides)

# add ON target scores based on spacer and protospacer sequence
list_pred_guides <- addOnTargetScores(
  list_pred_guides,
  methods = setdiff(score_methods, c("tssdist", "genrich"))
)

# rescale scores to a range between 0 and 1 if they exceed this range
rescale_score <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}
for (score in setdiff(score_methods, c("tssdist", "genrich"))) {
  score_val <- list_pred_guides@elementMetadata[[paste0("score_", score)]]
  if (min(score_val, na.rm = TRUE) < 0 || max(score_val, na.rm = TRUE) > 1) {
    list_pred_guides@elementMetadata[[paste0("score_", score)]] <- score_val %>%
      rescale_score()
  }
}

# add score based on distance to TSS; lower dist = higher score
list_pred_guides$score_tssdist <- 1 - (abs(list_pred_guides$dist_to_tss) / max(abs(tss_window)))

# add score based on G enrichment in seed region (pos -4 to -14 from PAM),
# see Miao & Jahn et al., The Plant Cell, 2023
list_pred_guides$score_genrich <- list_pred_guides$protospacer %>%
  subseq(start = spacer_length - 13, end = spacer_length - 3) %>%
  letterFrequency("G", as.prob = TRUE) %>%
  as.numeric() %>%
  rescale_score()

# add information on possible restriction sites
if (!is.null(restriction_sites)) {
  list_pred_guides <- addRestrictionEnzymes(
    list_pred_guides,
    patterns = restriction_sites,
    includeDefault = FALSE
  )
}

# add information on possible bad seed (seed region = 12 bp upstream of PAM)
for (seed_pattern in bad_seeds) {
  list_pred_guides@elementMetadata[[paste0("seed_", seed_pattern)]] <- vmatchPattern(
    pattern = seed_pattern,
    subject = subseq(list_pred_guides$protospacer, start = spacer_length - 11, end = spacer_length)
  ) %>%
    sapply(FUN = function(x) length(width(x))) %>%
    as.logical()
}


# STAGE 3 : FILTER GUIDES
# ------------------------------
# guides are filtered on A) hard boundaries
# and B) ranked by off-/on-target scores

# filter by strand
n_guides_pre_filter <- length(list_pred_guides)
if (strands != "both") {
  if (strands == "coding") {
    filter_strand <- strand(list_pred_guides) != list_pred_guides$tss_strand
  } else if (strands == "template") {
    filter_strand <- strand(list_pred_guides) == list_pred_guides$tss_strand
  } else {
    stop("parameter 'strands' must be one of 'coding', 'template', or 'both'")
  }
  messages <- append(messages, paste0(
    "Removed ", sum(!filter_strand),
    " guide RNAs not targeting the ", strands, " strand"
  ))
} else {
  filter_strand <- rep(TRUE, length(list_pred_guides))
}

# filter by GC content
filter_gc_low <- list_pred_guides$percentGC >= gc_content_range[1]
filter_gc_high <- list_pred_guides$percentGC <= gc_content_range[2]

messages <- append(messages, paste0(
  "Removed ", sum(!filter_gc_low),
  " guide RNAs with GC content below ", gc_content_range[1], "%"
))
messages <- append(messages, paste0(
  "Removed ", sum(!filter_gc_high),
  " guide RNAs with GC content above ", gc_content_range[2], "%"
))

# filter by poly-T stretches
filter_polyt <- !list_pred_guides$polyT
messages <- append(messages, paste0(
  "Removed ", sum(!filter_polyt),
  " guide RNAs with poly-T stretches"
))

# filter by starting-G stretches
filter_startg <- !list_pred_guides$startingGGGGG
messages <- append(messages, paste0(
  "Removed ", sum(!filter_startg),
  " guide RNAs with 5 starting G stretches"
))

# filter by off-target hits
if (filter_multi_targets) {
  filter_offtargets <- with(list_pred_guides, n0 == 1 & n1 == 0 & n2 == 0 & n3 == 0)
} else {
  filter_offtargets <- with(list_pred_guides, n0 >= 1 & n1 == 0 & n2 == 0 & n3 == 0)
}
messages <- append(messages, paste0(
  "Removed ", sum(!filter_offtargets),
  " guide RNAs with off targets (1 to 4 nt mismatches allowed)"
))

# filter by off-target scores (only if multi-targets are not allowed)
if (filter_multi_targets) {
  filter_offtarget_scores <- unname(list_pred_guides$score_cfd == 1 & list_pred_guides$score_mit == 1)
} else {
  filter_offtarget_scores <- rep_along(list_pred_guides, TRUE)
}
messages <- append(messages, paste0(
  "Removed ", sum(!filter_offtarget_scores),
  " guide RNAs with CFD or MIT off-target score < 1"
))

# filter by bad seed sequences
filter_badseeds <- list_pred_guides@elementMetadata[paste0("seed_", bad_seeds)] %>%
  as.data.frame() %>%
  apply(1, function(x) !any(x))
messages <- append(messages, paste0(
  "Removed ", sum(!filter_badseeds),
  " guide RNAs with a bad seed sequence"
))

# apply filters
filter_by_all <- filter_strand &
  filter_gc_low &
  filter_gc_high &
  filter_polyt &
  filter_startg &
  filter_offtargets &
  filter_offtarget_scores &
  filter_badseeds
list_pred_guides <- list_pred_guides[filter_by_all]

# make composite score (mean of all scores)
list_pred_guides$score_all <- list_pred_guides@elementMetadata[paste0("score_", score_methods)] %>%
  apply(1, weighted.mean, w = score_weights, na.rm = TRUE)

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

# this loop filters guides iteratively until no overlapping guides remain
guides_filtered <- filter_overlaps(list_pred_guides)
while (length(guides_filtered) > 0) {
  messages <- append(messages, paste0(
    "Removed ", length(guides_filtered),
    " guide RNAs overlapping with better ones"
  ))
  list_pred_guides <- list_pred_guides[-guides_filtered]
  guides_filtered <- filter_overlaps(list_pred_guides)
}

# filter duplicated guides (only when "filter_multi_targets == FALSE")
if (!filter_multi_targets) {
  filter_duplicated <- !duplicated(list_pred_guides$protospacer)
  filter_duplicated_revcom <- sapply(list_pred_guides$protospacer, function(x) {
    !(x %in% reverseComplement(list_pred_guides$protospacer))
  })
  list_pred_guides <- list_pred_guides[filter_duplicated & filter_duplicated_revcom]
  list_pred_guides$multi_target <- ifelse(list_pred_guides$n0 > 1, TRUE, FALSE)
  messages <- append(messages, paste0(
    "Removed ", sum(!(filter_duplicated & filter_duplicated_revcom)),
    " guide RNAs with duplicated protospacer sequence"
  ))
}

# apply filtering by score
if (!is.null(filter_top_n)) {
  index_top_n <- tapply(
    list_pred_guides$score_all %>% setNames(names(list_pred_guides)),
    list_pred_guides$tx_name,
    function(x) names(x)[order(x, decreasing = TRUE)][1:filter_top_n]
  ) %>%
    unlist() %>%
    unname() %>%
    na.omit()
  messages <- append(messages, paste0(
    "Removed ", length(list_pred_guides) - length(index_top_n),
    " guide RNAs with low rank for on-target scores"
  ))
  list_pred_guides <- list_pred_guides[names(list_pred_guides) %in% index_top_n]
} else if (!is.null(filter_score_threshold)) {
  index_score <- list_pred_guides$score_all >= filter_score_threshold
  list_pred_guides <- list_pred_guides[index_score]
  messages <- append(messages, paste0(
    "Removed ", length(list_pred_guides) - sum(index_score),
    " guide RNAs falling below the on-target score threshold of ",
    filter_score_threshold
  ))
}

# print final numbers
n_guides_removed <- n_guides_pre_filter - length(list_pred_guides)
messages <- append(messages, "----------------------------------")
messages <- append(messages, paste0(
  "Removed a total of ", n_guides_removed, " guide RNAs (",
  round(n_guides_removed / n_guides_pre_filter * 100, 1), "%)"
))
messages <- append(messages, paste0(
  "A final list of ", length(list_pred_guides), " guide RNAs was exported"
))

if (!is.null(target_region)) {
  list_tx_total <- list_tx[seqnames(list_tx) %in% target_region]$tx_name
} else {
  list_tx_total <- list_tx$tx_name
}
if (filter_rna) {
  list_rna_targets <- grepl("^rna-", list_tx_total)
  list_tx_total <- list_tx_total[!list_rna_targets]
}

list_no_guides <- setdiff(list_tx_total, unique(list_pred_guides$tx_name))
if (length(list_no_guides) >= 1) {
  messages <- append(messages, paste0(
    "guide RNAs were not found for ", length(list_no_guides),
    " targets: ", paste(list_no_guides, collapse = ", ")
  ))
}

# optionally create a set of random control guides
# we do a very reduced check of these guides (only off-targets)
# since all on-target scores don't apply and would produce errors
if (no_target_controls > 0) {
  list_ntc_guides <- sapply(1:no_target_controls, FUN = function(x) {
    ntc <- sample(x = c("A", "C", "T", "G"), size = spacer_length, replace = TRUE)
    ntc <- paste(ntc, collapse = "")
    names(ntc) <- paste0("NTC_", str_pad(x, width = nchar(max(no_target_controls)), pad = "0"))
    ntc
  })

  list_ntc_guides <- GuideSet(
    ids = names(list_ntc_guides),
    protospacers = unname(list_ntc_guides),
    seqnames = rep(seqnames(bsgenome)[1], no_target_controls),
    bsgenome = bsgenome,
    pams = rep(colnames(get(crispr_enzyme)@cutSites)[1], no_target_controls),
    CrisprNuclease = get(crispr_enzyme),
    targetOrigin = "customSequences",
    customSequences = genome_dna
  )

  list_ntc_guides <- addSpacerAlignments(
    list_ntc_guides,
    aligner = "biostrings",
    txObject = txdb,
    custom_seq = genome_dna,
    n_mismatches = 3,
    n_max_alignments = 1000,
    addSummary = TRUE,
    all_alignments = TRUE
  )

  df_ntc_guides <- list_ntc_guides %>%
    as.data.frame() %>%
    as_tibble()
}

# STAGE 4 : EXPORT RESULTS
# ------------------------------
# prepare result table
df_pred_guides <- list_pred_guides %>%
  as.data.frame() %>%
  as_tibble() %>%
  dplyr::select(-tssAnnotation) %>%
  mutate(
    width = spacer_length,
    start = ifelse(strand == "-", start + 1, start - 20),
    end = ifelse(strand == "-", end + 20, end - 1)
  )

# add potential linkers on 5' and/or 3' end
if (!is.null(fiveprime_linker) || !is.null(threeprime_linker)) {
  df_pred_guides <- df_pred_guides %>%
    mutate(
      fiveprime_linker = fiveprime_linker,
      threeprime_linker = threeprime_linker,
      seq_with_linkers = paste0(fiveprime_linker, protospacer, threeprime_linker)
    )
  if (no_target_controls > 0) {
    df_ntc_guides <- df_ntc_guides %>%
      mutate(
        fiveprime_linker = fiveprime_linker,
        threeprime_linker = threeprime_linker,
        seq_with_linkers = paste0(fiveprime_linker, protospacer, threeprime_linker)
      )
  }
}

# export table with NTCs
if (no_target_controls > 0) {
  write_csv(df_ntc_guides, output_ntc)
}

# export results as csv table
write_csv(df_pred_guides, output_top)

# export table with transcripts where no guide is available
df_no_guides <- list_tx[list_tx$tx_name %in% list_no_guides] %>%
  as.data.frame()
write_csv(df_no_guides, output_fail)

# export log
write_lines(
  file = snakemake@log[["path"]],
  x = paste0("DESIGN_GUIDES: ", messages)
)
