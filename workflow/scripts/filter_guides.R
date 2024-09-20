# LOAD PACKAGES
# ------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(Biostrings)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(crisprBase)
  library(crisprDesign)
})


# CONFIGURATION
# ------------------------------
save(snakemake, file = "sm.Rdata")
source(paste0(snakemake@scriptdir, "/utilities.R"))
sm_params <- snakemake@params[[1]]
max_cores <- snakemake@params[[2]]
for (param in names(sm_params)) {
  assign(param, sm_params[[param]], envir = .GlobalEnv)
}
target_region <- snakemake@config$design_guides$target_region
spacer_length <- snakemake@config$design_guides$spacer_length
strands <- snakemake@config$design_guides$strands
bad_seeds <- snakemake@config$design_guides$bad_seeds
score_methods <- snakemake@config$design_guides$score_methods
score_weights <- snakemake@config$design_guides$score_weights
restriction_sites <- snakemake@config$design_guides$restriction_sites
input_file <- snakemake@input$guideset
input_tx <- snakemake@input$list_tx
output_csv <- snakemake@output$guideset_csv
output_gff <- snakemake@output$guideset_gff


# STAGE 1 : FILE PREPARATION
# ------------------------------
# import guideset
messages <- c(paste0("Importing guide set: ", input_file))
target_type <- str_remove_all(tail(str_split_1(input_file, "\\/"), 1), "guideRNAs_|.RData")
load(input_file)
list_guides <- get(paste0("list_", target_type))
list_tx <- read_csv(input_tx, show_col_types = FALSE)


# STAGE 2 : FILTER GUIDES BY PROPERTY
# -----------------------------------
# guides are filtered on A) hard boundaries

# remove guides targeting RNAs (optional)
if (filter_rna & "tx_name" %in% colnames(list_guides@elementMetadata)) {
  list_rna_targets <- grepl("^rna-", list_guides$tx_name)
  list_guides <- list_guides[!list_rna_targets]
}

# filter by strand; this is gene-orientation specific for targets
n_guides_pre_filter <- length(list_guides)
if (strands != "both") {
  if (strands == "coding") {
    if (target_type == "target") {
      filter_strand <- strand(list_guides) != list_guides$tss_strand
    } else if (target_type == "intergenic") {
      filter_strand <- strand(list_guides) == "+"
    } else if (target_type == "ntc") {
      filter_strand <- rep(TRUE, length(list_guides))
    }
  } else if (strands == "template") {
    if (target_type == "target") {
      filter_strand <- strand(list_guides) == list_guides$tss_strand
    } else if (target_type == "intergenic") {
      filter_strand <- strand(list_guides) == "-"
    } else if (target_type == "ntc") {
      filter_strand <- rep(TRUE, length(list_guides))
    }
  } else {
    stop("parameter 'strands' must be one of 'coding', 'template', or 'both'")
  }
  messages <- append(messages, paste0(
    "Removed ", sum(!filter_strand),
    " guide RNAs not targeting the ", strands, " strand"
  ))
} else {
  filter_strand <- rep(TRUE, length(list_guides))
}

# filter by GC content
filter_gc_low <- list_guides$percentGC >= gc_content_range[1]
filter_gc_high <- list_guides$percentGC <= gc_content_range[2]
messages <- append(messages, paste0(
  "Removed ", sum(!filter_gc_low),
  " guide RNAs with GC content below ", gc_content_range[1], "%"
))
messages <- append(messages, paste0(
  "Removed ", sum(!filter_gc_high),
  " guide RNAs with GC content above ", gc_content_range[2], "%"
))

# filter by poly-T stretches
filter_polyt <- !list_guides$polyT
messages <- append(messages, paste0(
  "Removed ", sum(!filter_polyt),
  " guide RNAs with poly-T stretches"
))

# filter by starting-G stretches
filter_startg <- !list_guides$startingGGGGG
messages <- append(messages, paste0(
  "Removed ", sum(!filter_startg),
  " guide RNAs with 5 starting G stretches"
))

# filter by off-target hits
if (target_type == "ntc") {
  filter_offtargets <- with(list_guides, n0 == 0 & n1 == 0 & n2 == 0 & n3 == 0)
} else {
  if (filter_multi_targets) {
    filter_offtargets <- with(list_guides, n0 == 1 & n1 == 0 & n2 == 0 & n3 == 0)
  } else {
    filter_offtargets <- with(list_guides, n0 >= 1 & n1 == 0 & n2 == 0 & n3 == 0)
  }
}
messages <- append(messages, paste0(
  "Removed ", sum(!filter_offtargets),
  " guide RNAs with off targets (matching with up to 4 nt mismatches)"
))

# filter by off-target scores (only if multi-targets are not allowed)
if (filter_multi_targets) {
  filter_offtarget_scores <- unname(list_guides$score_cfd == 1 & list_guides$score_mit == 1)
} else {
  filter_offtarget_scores <- rep_along(list_guides, TRUE)
}
messages <- append(messages, paste0(
  "Removed ", sum(!filter_offtarget_scores),
  " guide RNAs with CFD or MIT off-target score < 1"
))

# filter by bad seed sequences
filter_badseeds <- list_guides@elementMetadata[paste0("seed_", bad_seeds)] %>%
  as.data.frame() %>%
  apply(1, function(x) !any(x))
messages <- append(messages, paste0(
  "Removed ", sum(!filter_badseeds),
  " guide RNAs with a bad seed sequence"
))

# filter by resctriction sites
if (!is.null(restriction_sites)) {
  filter_restriction <- list_guides@elementMetadata$enzymeAnnotation %>%
    as.data.frame() %>%
    dplyr::select(-c(1, 2)) %>%
    apply(1, function(x) !any(x))
  messages <- append(messages, paste0(
    "Removed ", sum(!filter_restriction),
    " guide RNAs matching a restriction site"
  ))
} else {
  filter_restriction <- rep(TRUE, length(list_guides))
}

# apply filters
filter_by_all <- filter_strand &
  filter_gc_low &
  filter_gc_high &
  filter_polyt &
  filter_startg &
  filter_offtargets &
  filter_offtarget_scores &
  filter_badseeds &
  filter_restriction
list_guides <- list_guides[filter_by_all]


# STAGE 3 : FILTER GUIDES BY SCORE
# -----------------------------------
# guides are filtered on B) on-target score/ranking
# this does not apply for NTC controls
if (target_type != "ntc") {
  # make composite score (mean of all scores)
  list_guides$score_all <- list_guides@elementMetadata[paste0("score_", score_methods)] %>%
    apply(1, weighted.mean, w = score_weights, na.rm = TRUE)

  # this loop filters guides iteratively until no overlapping guides remain
  guides_filtered <- filter_overlaps(list_guides)
  round_counter <- 1
  while (length(guides_filtered) > 0) {
    messages <- append(messages, paste0(
      "Removed ", length(guides_filtered),
      " guide RNAs overlapping with better ones",
      " (round ", round_counter, ")"
    ))
    list_guides <- list_guides[-guides_filtered]
    guides_filtered <- filter_overlaps(list_guides)
    round_counter <- round_counter + 1
  }

  # filter duplicated guides (only when "filter_multi_targets == FALSE")
  if (!filter_multi_targets) {
    filter_duplicated <- !duplicated(list_guides$protospacer)
    filter_duplicated_revcom <- sapply(list_guides$protospacer, function(x) {
      !(x %in% reverseComplement(list_guides$protospacer))
    })
    list_guides <- list_guides[filter_duplicated & filter_duplicated_revcom]
    list_guides$multi_target <- ifelse(list_guides$n0 > 1, TRUE, FALSE)
    messages <- append(messages, paste0(
      "Removed ", sum(!(filter_duplicated & filter_duplicated_revcom)),
      " guide RNAs with duplicated protospacer sequence"
    ))
  }

  # apply filtering by score
  # option 1 is top N guides per gene or per kb for intergenic regions
  filter_top_n <- ifelse(
    target_type == "target",
    filter_best_per_gene,
    filter_best_per_tile
  )
  if (!is.null(filter_top_n)) {
    index_top_n <- tapply(
      list_guides$score_all %>% setNames(names(list_guides)),
      list_guides$tx_name,
      function(x) names(x)[order(x, decreasing = TRUE)][1:filter_top_n]
    ) %>%
      unlist() %>%
      unname() %>%
      na.omit()
    messages <- append(messages, paste0(
      "Removed ", length(list_guides) - length(index_top_n),
      " guide RNAs with low rank for on-target scores"
    ))
    list_guides <- list_guides[names(list_guides) %in% index_top_n]
  } else if (!is.null(filter_score_threshold)) {
    # option 2 is by score threshold
    index_score <- list_guides$score_all >= filter_score_threshold
    list_guides <- list_guides[index_score]
    messages <- append(messages, paste0(
      "Removed ", length(list_guides) - sum(index_score),
      " guide RNAs falling below the on-target score threshold of ",
      filter_score_threshold
    ))
  } else {
    warn_no_filter <- paste0(
      "Warning: Removed no guide RNAs by score! ",
      "Neither 'filter_best_per_gene/tile' nor 'filter_score_threshold' are specified."
    )
    messages <- append(messages, warn_no_filter)
    message(warn_no_filter)
  }
}

# print final numbers
n_guides_removed <- n_guides_pre_filter - length(list_guides)
messages <- append(messages, "----------------------------------")
messages <- append(messages, paste0(
  "Removed a total of ", n_guides_removed, " guide RNAs (",
  round(n_guides_removed / n_guides_pre_filter * 100, 1), "%)"
))
messages <- append(messages, paste0(
  "A final list of ", length(list_guides), " guide RNAs was exported"
))

# STAGE 4 : DETERMINE FAILURES
# ------------------------------
#
if (target_type == "target") {
  if (!is.null(target_region)) {
    list_tx_total <- filter(list_tx, seqnames %in% target_region)$tx_name
  } else {
    list_tx_total <- list_tx$tx_name
  }
  if (filter_rna) {
    list_rna_targets <- grepl("^rna-", list_tx_total)
    list_tx_total <- list_tx_total[!list_rna_targets]
  }

  list_no_guides <- setdiff(list_tx_total, unique(list_guides$tx_name))
  if (length(list_no_guides) >= 1) {
    messages <- append(messages, paste0(
      "guide RNAs were not found for ", length(list_no_guides),
      " targets: ", paste(list_no_guides, collapse = ", ")
    ))
  }
}

# STAGE 4 : EXPORT RESULTS
# ------------------------------
# prepare result table
df_guides <- list_guides %>%
  as.data.frame() %>%
  as_tibble() %>%
  dplyr::select(!any_of(c("tssAnnotation"))) %>%
  mutate(
    enzymeAnnotation = FALSE,
    width = spacer_length,
    start = ifelse(strand == "-", start + 1, start - 20),
    end = ifelse(strand == "-", end + 20, end - 1)
  )

# add potential linkers on 5' and/or 3' end
if (!is.null(fiveprime_linker) || !is.null(threeprime_linker)) {
  df_guides <- df_guides %>%
    mutate(
      fiveprime_linker = fiveprime_linker,
      threeprime_linker = threeprime_linker,
      seq_with_linkers = paste0(fiveprime_linker, protospacer, threeprime_linker)
    )
}

# check if guides remain, otherwise throw warning
if (!nrow(df_guides)) {
  warn_no_guides <- "the final list of guide RNAs after filtering is 0, omitting export."
  messages <- append(messages, warn_no_guides)
  warning(warn_no_guides)
  export_as_gff <- FALSE
}

# export results as CSV table and GFF file (optional)
write_csv(df_guides, output_csv)
if (export_as_gff) {
  export_gff(
    df = df_guides,
    filename = str_replace(output_gff, ".csv$", ".gff"),
    feature = paste0("sgRNA_", target_type)
  )
}

# export table with transcripts where no guide is available
if (target_type == "target") {
  df_no_guides <- filter(list_tx, tx_name %in% list_no_guides) %>%
    as.data.frame()
  write_csv(df_no_guides, str_replace(output_csv, "target", "target_failed"))
}

# export log
write_lines(
  file = snakemake@log[["path"]],
  x = paste0("DESIGN_GUIDES: ", messages)
)
