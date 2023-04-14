# LOAD PACKAGES
# ------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(GenomeInfoDbData)
  library(GenomicFeatures)
  library(crisprBase)
  library(crisprDesign)
  library(crisprBwa)
  library(crisprScore)
  library(ORFik)
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


# STAGE 1 : FILE PREPARATION
# ------------------------------
# import genome sequence
messages <- c("importing genome sequence and annotation")
genome_fasta <- snakemake@input[["fasta"]]
genome_gff <- snakemake@input[["gff"]]
genome_dna <- Biostrings::readDNAStringSet(genome_fasta)
genome_name <- str_remove_all(names(genome_dna)[1], "^NC\\_[0-9]*\\.[0-9]* |\\,.*")
seqinfo_genome <- seqinfo(genome_dna)
isCircular(seqinfo_genome) <- rep_along(seqlevels(seqinfo_genome), FALSE)
genome(seqinfo_genome) <- genome_name

# import genome annotation
txdb <- makeTxDbFromGFF(
  file = genome_gff
)

# check if sequence annotation is identical for sequence and annotation
if (any(seqlevels(seqinfo_genome) != seqlevels(txdb))) {
  seqlevels(seqinfo_genome) <- seqlevels(txdb)
  seqinfo(genome_dna) <- seqinfo_genome
}

# re-import genome annotation with chromosome metadata
txdb <- makeTxDbFromGFF(
  file = genome_gff,
  organism = genome_name,
  chrominfo = seqinfo_genome
)

# extract all transcripts
list_cds <- ORFik::loadRegion(txdb, "cds")
if (!is.null(gene_subset)) {
  list_cds <- list_cds[1:gene_subset]
}

# optionally extract region surrounding each start codon
# to target probable promoters
list_cds_leader <- extendLeaders(list_cds, extension = five_utr)


# STAGE 2 : FINDING GUIDES
# ------------------------------
# extract sequences for all CDS start regions
list_cds_seq <- GenomicFeatures::extractTranscriptSeqs(
  genome_dna,
  list_cds_leader
)

# predict all possible guides for target regions
messages <- append(
  messages,
  paste0("predicting guide RNAs for enzyme: ", crispr_enzyme)
)
data(list = c(crispr_enzyme), package = "crisprBase")
list_pred_guides <- crisprDesign::findSpacers(
  x = list_cds_seq,
  canonical = canonical,
  both_strands = both_strands,
  spacer_len = spacer_length,
  crisprNuclease = get(crispr_enzyme)
)

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
      crisprDesign::addSpacerAlignments(
        guide_chunk,
        aligner = guide_aligner,
        txObject = txdb,
        custom_seq = genome_dna,
        n_mismatches = 4,
        n_max_alignments = 1000,
        addSummary = TRUE,
        all_alignments = TRUE
      )
    }
  )

  # merge all chunks to single guideSet again
  list_pred_guides <- Reduce("c", list_pred_guides_chunks)
} else {
  stop("Only 'biostrings' currently supported for off target search;
    methods 'bowtie' and 'bwa' not yet implemented.")
}

# add OFF target scores based on alignment
list_pred_guides <- crisprDesign::addOffTargetScores(list_pred_guides)

# add ON target scores based on spacer and protospacer sequence
list_pred_guides <- crisprDesign::addOnTargetScores(
  list_pred_guides,
  methods = score_methods
)

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

# filter by GC content
n_guides_pre_filter <- length(list_pred_guides)
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
filter_offtargets <- with(list_pred_guides, n0 == 1 & n1 == 0 & n2 == 0 & n3 == 0 & n4 == 0)
messages <- append(messages, paste0(
  "Removed ", sum(!filter_offtargets),
  " guide RNAs with off targets (1 to 4 nt mismatches allowed)"
))

# filter by off-target scores
filter_offtarget_scores <- unname(list_pred_guides$score_cfd == 1 & list_pred_guides$score_cfd == 1)
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
filter_by_all <- filter_gc_low &
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

# apply filtering by score
if (!is.null(filter_top_n)) {
  index_top_n <- tapply(
    list_pred_guides$score_all %>% setNames(names(list_pred_guides)),
    list_pred_guides$region,
    function(x) names(x)[order(x, decreasing = TRUE)][1:10]
  ) %>%
    unlist() %>%
    unname()
  messages <- append(messages, paste0(
    "Removed ", length(list_pred_guides) - length(index_top_n),
    " guide RNAs with low rank for on-target scores"
  ))
  list_pred_guides <- list_pred_guides[names(list_pred_guides) %in% index_top_n]
} else if (!is.null(filter_score_threshold)) {
  messages <- append(messages, "Filtering by on-target score threshold not yet implemented")
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

# print final numbers
n_guides_removed <- n_guides_pre_filter - length(list_pred_guides)
messages <- append(messages, "----------------------------------")
messages <- append(messages, paste0(
  "Removed a total of ", n_guides_removed, " guide RNAs (",
  round(n_guides_removed / n_guides_pre_filter * 100, 1), "%)"
))
messages <- append(messages, paste0(
  "A final list of ", length(list_pred_guides), " guide RNAs is exported"
))


# STAGE 4 : EXPORT RESULTS
# ------------------------------
# export results as GuideSet and data.frame / csv
df_pred_guides <- list_pred_guides %>%
  as.data.frame() %>%
  as_tibble()
write_csv(df_pred_guides, output_top)

write_lines(
  file = snakemake@log[["path"]],
  x = paste0("DESIGN_GUIDES: ", messages)
)
