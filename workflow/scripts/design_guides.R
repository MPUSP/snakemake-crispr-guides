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
source("workflow/scripts/utilities.R")
sm_params <- snakemake@params[[1]]
max_cores <- snakemake@params[[2]]
for (param in names(sm_params)) {
  assign(param, sm_params[[param]], envir = .GlobalEnv)
}
output_files <- paste0("guideRNAs_", target_type, ".RData")
output_tx <- snakemake@output$list_tx


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
list_guides <- findSpacers(
  x = genome_dna,
  canonical = canonical,
  both_strands = TRUE,
  spacer_len = spacer_length,
  crisprNuclease = get(crispr_enzyme)
)

# pcocess portion of guides matching supplied TSS windows
if ("target" %in% target_type) {
  # import genome annotation with chromosome metadata
  txdb <- makeTxDbFromGFF(
    file = genome_gff,
    organism = unname(genome(seqinfo_genome)[1]),
    chrominfo = seqinfo_genome
  )

  # extract (pseudo) transcription start sites (TSS)
  list_tx <- transcripts(txdb)
  list_tss <- resize(list_tx, width = 1)
  list_tss$ID <- list_tss$tx_name

  # map guides to TSS windows
  list_guides <- addTssAnnotation(
    list_guides,
    list_tss,
    tss_window = tss_window
  )

  # for targets: remove all guides that do not target within a TSS window
  index_non_target <- unlist(mclapply(
    mc.cores = max_cores,
    X = list_guides$tssAnnotation,
    FUN = nrow
  ))
  list_target <- list_guides[index_non_target != 0]

  # clip TSS windows where the 5'-UTR extends into another transcript
  # then remove all guides that are beyond the 3'- or 5' end of the target gene
  list_tg_intergenic <-
    {
      start(list_tx) - c(1, end(list_tx[-length(list_tx)]))
    } %>%
    c(., unname(tail(seqlengths(txdb), 1) - tail(end(list_tx), 1))) %>%
    replace(., . > abs(tss_window[1]), abs(tss_window[1])) %>%
    replace(., . < tss_window[1] / 2, tss_window[1] / 2)
  max_5UTR_index <- seq_along(list_tx)
  max_5UTR_index[as.logical(strand(list_tx) == "-")] <-
    max_5UTR_index[as.logical(strand(list_tx) == "-")] + 1
  max_5UTR <- list_tg_intergenic[max_5UTR_index]
  list_tx_window <- flank(list_tx, max_5UTR)
  list_tx_window <- punion(list_tx, list_tx_window)
  list_tss_window <- punion(
    flank(list_tss, abs(tss_window[1])),
    resize(list_tss, width = tss_window[2])
  )
  genome(list_tss_window) <- "custom"
  genome(list_tx_window) <- "custom"

  # for intergenic regions: remove all guides that target within transcripts
  index_intergenic <- findOverlaps(list_guides, list_tx_window, ignore.strand = TRUE)
  list_intergenic <- list_guides[-index_intergenic@from]

  # select only guides whose PAM overlaps with the TSS window
  list_tss_window <- pintersect(list_tss_window, list_tx_window)
  list_overlaps <- as.list(findOverlaps(
    list_target, list_tss_window,
    ignore.strand = TRUE
  ))
  names(list_overlaps) <- names(list_target)
  list_targets_keep <- unlist(lapply(list_overlaps, length) == 1)
  list_target <- list_target[list_targets_keep]
  list_overlaps <- list_overlaps[list_targets_keep]

  # trim original TSS annotation to one entry per guide
  df_tss_annot <- list_target$tssAnnotation %>%
    as_tibble() %>%
    right_join(
      by = c("group_name", "tx_id"),
      enframe(
        unlist(list_overlaps),
        name = "group_name",
        value = "tx_id"
      )
    ) %>%
    dplyr::select(-group, -chr, -strand) %>%
    mutate(tx_width = setNames(width(list_tx), list_tx$tx_id)[tx_id]) %>%
    arrange(as.numeric(str_sub(group_name, 8L)))

  # add modified TSS information to metadata
  list_target@elementMetadata <- cbind(
    list_target@elementMetadata,
    df_tss_annot
  )

  # remove guides where a transcript was mapped but not the right one
  list_target <- list_target[!is.na(list_target$tss_id)]

  # add spacer and PAM sequence features
  list_target <- addSequenceFeatures(list_target)

  # add off target candidates with biostrings
  list_target <- addOffTargetSequences(list_target, guide_aligner, max_cores, txdb, genome_dna, genome_index)

  # add OFF target scores based on alignment
  list_target <- addOffTargetScores(list_target)

  # add ON target scores based on spacer and protospacer sequence
  list_target <- addAllOnTargetScores(list_target, score_methods, spacer_length, tss_window)

  # add information on possible restriction sites
  if (!is.null(restriction_sites)) {
    list_target <- addRestrictionEnzymes(list_target, patterns = restriction_sites, includeDefault = FALSE)
  }

  # add information on possible bad seed (seed region = 12 bp upstream of PAM)
  list_target <- addBadSeeds(list_target, bad_seeds)
} else if ("intergenic" %in% target_type) {
  list_intergenic <- list_guides
}

if ("intergenic" %in% target_type) {
  # a similar procedure for benchmarking intergenic guides is applied as
  # as for targets (e.g. genes)
  # first, to filter intergenic regions by score later on, we apply
  # a tiling approach: guides are grouped into bins of 1000 nt width
  # for later filtering
  list_intergenic$tx_window <- as_tibble(list_intergenic) %>%
    group_by(seqnames) %>%
    mutate(tx_window = cut_interval(start, length = tiling_window, labels = FALSE)) %>%
    pull(tx_window)
  list_intergenic$tx_name <- with(
    list_intergenic,
    paste0(
      seqnames, "_",
      formatC(tx_window, width = nchar(as.character(max(tx_window))), flag = "0")
    )
  )
  list_intergenic <- addSequenceFeatures(list_intergenic)
  list_intergenic <- addOffTargetSequences(list_intergenic, guide_aligner, max_cores, txdb, genome_dna, genome_index)
  list_intergenic <- addOffTargetScores(list_intergenic)
  list_intergenic <- addAllOnTargetScores(list_intergenic, score_methods, spacer_length, tss_window = NULL)
  if (!is.null(restriction_sites)) {
    list_intergenic <- addRestrictionEnzymes(list_intergenic, patterns = restriction_sites, includeDefault = FALSE)
  }
  list_intergenic <- addBadSeeds(list_intergenic, bad_seeds)
}

# optionally create a set of random control guides
# we do a very reduced check of these guides (sequences and off-targets)
# since all on-target scores don't apply and would produce errors
if ("ntc" %in% target_type & no_target_controls > 0) {
  list_ntc <- sapply(1:no_target_controls, FUN = function(x) {
    ntc <- sample(x = c("A", "C", "T", "G"), size = spacer_length, replace = TRUE)
    ntc <- paste(ntc, collapse = "")
    names(ntc) <- paste0("NTC_", str_pad(x, width = nchar(max(no_target_controls)), pad = "0"))
    ntc
  })

  list_ntc <- GuideSet(
    ids = names(list_ntc),
    protospacers = unname(list_ntc),
    seqnames = rep(seqnames(bsgenome)[1], no_target_controls),
    bsgenome = bsgenome,
    pams = rep(colnames(get(crispr_enzyme)@cutSites)[1], no_target_controls),
    CrisprNuclease = get(crispr_enzyme),
    targetOrigin = "customSequences",
    customSequences = genome_dna
  )

  list_ntc <- addSequenceFeatures(list_ntc)
  list_ntc <- addOffTargetSequences(list_ntc, guide_aligner, max_cores, txdb, genome_dna, genome_index)
  list_ntc <- suppressWarnings(addOffTargetScores(list_ntc))
  list_ntc$score_cfd <- with(list_ntc, replace(score_cfd, is.na(score_cfd), 1))
  list_ntc$score_mit <- with(list_ntc, replace(score_mit, is.na(score_mit), 1))
  if (!is.null(restriction_sites)) {
    list_ntc <- addRestrictionEnzymes(list_ntc, patterns = restriction_sites, includeDefault = FALSE)
  }
  list_ntc <- addBadSeeds(list_ntc, bad_seeds)
}

# export final annotated guide sets for next module
output_dir <- "results/design_guides/"
for (target in target_type) {
  save(
    list = c(paste0("list_", target)),
    file = paste0(output_dir, "guideRNAs_", target, ".RData")
  )
}

# export list of transcripts/target genes for filtering module
if ("target" %in% target_type) {
  df_tx <- as_tibble(list_tx)
  if (!is.null(target_region)) {
    df_tx <- filter(df_tx, seqnames %in% target_region)
  }
  write_csv(df_tx, output_tx)
} else {
  write_csv(data.frame(x = "dummy"), output_tx)
}

# export log
write_lines(
  file = snakemake@log[["path"]],
  x = paste0("DESIGN_GUIDES: ", messages)
)
