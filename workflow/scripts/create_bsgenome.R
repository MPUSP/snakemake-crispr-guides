# LOAD PACKAGES
# ------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(devtools)
  library(Biostrings)
  library(GenomicFeatures)
  library(BSgenome)
})

# IMPORT GENOME
# ------------------------------
output <- snakemake@output[["path"]]
messages <- c("importing genome sequence")
genome_fasta <- snakemake@input[["fasta"]]
genome_gff <- snakemake@input[["gff"]]
genome_dna <- Biostrings::readDNAStringSet(genome_fasta)
genome_seqlevels <- str_split_i(names(genome_dna), " ", i = 1)
seqinfo_genome <- seqinfo(genome_dna)
seqlevels(seqinfo_genome) <- genome_seqlevels
isCircular(seqinfo_genome) <- rep_along(seqlevels(seqinfo_genome), FALSE)

df_taxonomy <- loadTaxonomyDb()
header_gff <- read_lines(genome_gff, n_max = 10)
header_species <- na.omit(str_match(header_gff, ";species=.*"))[1]

if (!is.na(header_species)) {
  messages <- c(messages, "trying to guess taxonomy ID from GFF file")
  taxon_id <- header_species %>%
    str_extract("wwwtax.cgi(%3F|\\?)id(%3D|=)[0-9]+") %>%
    str_extract("[0-9]+$") %>%
    as.numeric()
  messages <- c(messages, paste("extracted taxon ID:", taxon_id))
  tax_entry <- filter(df_taxonomy, tax_id == taxon_id)
  if (!is.na(tax_entry[1, "genus"])) {
    genome_name <- paste(tax_entry[1, "genus"], tax_entry[1, "species"])
    messages <- c(messages, paste("found genome name:", genome_name)) 
  } else {
    genome_name <- NA
    messages <- c(messages, paste0("the extracted Taxonomy ID '", taxon_id, "' is not valid"))
  }
} else {
  genome_name <- NA
}

if (is.na(genome_name)) {
  messages <- c(messages, "trying to guess genome name from fasta file")
  genome_name <- str_remove_all(
    names(genome_dna)[1],
    paste0(genome_seqlevels[1], " |chromosome|\\,.*")
  )
  messages <- c(messages, paste("extracted genome name:", genome_name))
  df_taxonomy$name <- paste(df_taxonomy$genus, df_taxonomy$species)
  tax_entry <- filter(df_taxonomy, name == genome_name)
  taxon_id <- as.numeric(tax_entry[1, "tax_id"])
  if (!is.na(taxon_id)) {
    messages <- c(messages, paste("found taxon ID:", taxon_id))
  } else {
    genome_name <- NA
  }
}

if (is.na(genome_name)) {
  messages <- c(messages, "taxonomy guessing failed: falling back to arbitrary taxon ID")
  genome_name <- "root"
  taxon_id <- 1
}

genome_common <- str_flatten(str_split_1(genome_name, " ")[1:2], " ", na.rm = TRUE)
genome(seqinfo_genome) <- genome_name
genome_build_tag <- na.omit(str_match(header_gff, "#!genome-build [a-zA-Z0-9]+"))[1]
if (!is.na(genome_build_tag)) {
  genome_build <- str_extract(genome_build_tag, "[a-zA-Z0-9]+$")
} else {
  genome_build <- ""
}

# import genome annotation
quiet_txdb <- quietly(makeTxDbFromGFF)
txdb_msg <- quiet_txdb(file = genome_gff)
txdb <- txdb_msg$result
messages <- append(messages, paste0("warning: ", txdb_msg$warnings))
messages <- append(messages, txdb_msg$messages)

# check if sequence annotation is identical for sequence and annotation
if (!all(seqlevels(txdb) %in% seqlevels(seqinfo_genome))) {
  stop(
    paste0(
      "the following chromosome(s) are annotated in GFF file but not in FASTA sequence: ",
      setdiff(seqlevels(seqinfo_genome), seqlevels(txdb))
    )
  )
} else {
  seqinfo(genome_dna) <- seqinfo_genome
}

# re-import genome annotation with chromosome metadata
txdb <- quiet_txdb(
  file = genome_gff,
  organism = genome_name,
  chrominfo = seqinfo_genome,
  taxonomyId = taxon_id
)$result

# export sequence file
messages <- append(messages, "exporting sequence files")
dir.create(paste0(output, "/seqs_srcdir"))
txdb_metadata <- metadata(txdb) %>% deframe()
for (seqs in seqlevels(genome_dna)) {
  genome_dna[seqs] %>%
    writeXStringSet(paste0(output, "/seqs_srcdir/", seqs, ".fa"), format = "fasta")
}

messages <- append(messages, "forging *.2bit file for BSgenome preparation")
forgeSeqFiles(
  provider = "NCBI",
  genome = genome_build,
  seqnames = txdb$user_seqlevels,
  suffix = ".fa",
  seqs_srcdir = paste0(output, "/seqs_srcdir"),
  seqs_destdir = paste0(output, "/seqs_srcdir"),
  ondisk_seq_format = "2bit",
  verbose = FALSE
)

# create the BSgenome object, based on GH issue:
# https://github.com/Bioconductor/BSgenome/issues/3
messages <- append(messages, "creating BSgenome object")
bsgenome <- BSgenome:::BSgenome(
  organism = genome_common,
  common_name = genome_common,
  provider = "NCBI",
  provider_version = genome_build,
  release_date = txdb_metadata["Creation time"],
  release_name = txdb_metadata["Data source"],
  source_url = txdb_metadata["Data source"],
  seqnames = seqlevels(genome_dna),
  circ_seqs = NA,
  mseqnames = NULL,
  seqs_pkgname = NA_character_,
  seqs_dirpath = paste0(output, "/seqs_srcdir")
)

messages <- append(messages, paste0("exporting BSgenome object to ", output, "/bsgenome.RData"))
save(bsgenome, file = paste0(output, "/bsgenome.RData"))
messages <- append(messages, paste0("exporting sequence metadata to ", output, "/seqinfo.RData"))
save(seqinfo_genome, file = paste0(output, "/seqinfo.RData"))

# export log file
write_lines(
  file = snakemake@log[["path"]],
  x = paste0("BSGENOME: ", messages)
)
