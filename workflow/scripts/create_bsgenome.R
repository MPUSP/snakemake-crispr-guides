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
genome_seqlevels <- str_extract(names(genome_dna), "^NC\\_[0-9]*\\.[0-9]*")
genome_name <- str_remove_all(names(genome_dna)[1], "^NC\\_[0-9]*\\.[0-9]* |\\,.*")
genome_common <- str_flatten(str_split_1(genome_name, " ")[1:2], " ")
seqinfo_genome <- seqinfo(genome_dna)
seqlevels(seqinfo_genome) <- genome_seqlevels
isCircular(seqinfo_genome) <- rep_along(seqlevels(seqinfo_genome), FALSE)
genome(seqinfo_genome) <- genome_name
genome_build <- str_remove(read_lines(genome_gff, n_max = 5)[4], "#!genome-build ")

# import genome annotation
txdb <- makeTxDbFromGFF(
  file = genome_gff
)

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
txdb <- makeTxDbFromGFF(
  file = genome_gff,
  organism = genome_name,
  chrominfo = seqinfo_genome
)

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
  verbose = TRUE
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
