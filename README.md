# snakemake-crispr-guides

![Platform](https://img.shields.io/badge/platform-all-green)
[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/MPUSP/snakemake-crispr-guides/workflows/Tests/badge.svg?branch=main)](https://github.com/MPUSP/snakemake-crispr-guides/actions?query=branch%3Amain+workflow%3ATests)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)

---

A Snakemake workflow for the design of small guide RNAs (sgRNAs) for CRISPR applications.

- [snakemake-crispr-guides](#snakemake-crispr-guides)
  - [Usage](#usage)
  - [Workflow overview](#workflow-overview)
  - [Installation](#installation)
    - [Snakemake](#snakemake)
    - [Additional tools](#additional-tools)
  - [Running the workflow](#running-the-workflow)
    - [Input data](#input-data)
    - [Execution](#execution)
    - [Parameters](#parameters)
    - [Off-target scores](#off-target-scores)
    - [On-target scores](#on-target-scores)
  - [Output](#output)
  - [Authors](#authors)
  - [References](#references)


## Usage

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=MPUSP%2Fsnakemake-crispr-guides).

If you use this workflow in a paper, don't forget to give credits to the author(s) by citing the URL of this (original) repository and its DOI (see above).

## Workflow overview

<!-- include logo-->
<!-- <img src="docs/images/logo.png" align="right" /> -->

---

This workflow is a best-practice workflow for the automated generation of guide RNAs for CRISPR applications. It's main purpose is to provide a simple, efficient and easy-to-use framework to design thousands of guides simultaneously for CRISPR libraries from as little input as an organism's name/genome ID. For the manual design of single guides, users are instead referred to even simpler web resources such as [Chop-Chop](http://chopchop.cbu.uib.no/), [CRISPick](https://portals.broadinstitute.org/gppx/crispick/public), or [Cas-OFFinder/Cas-Designer](http://www.rgenome.net/cas-designer/).

This workflow relies to a large degree on the underlying [Bioconductor package ecosystem `crisprVerse`](http://bioconductor.org/packages/release/bioc/html/crisprVerse.html), published in 2022 by:

> Hoberecht, L., Perampalam, P., Lun, A. et al. _A comprehensive Bioconductor ecosystem for the design of CRISPR guide RNAs across nucleases and technologies_. Nat Commun 13, 6568 (**2022**). https://doi.org/10.1038/s41467-022-34320-7.

The workflow is built using [snakemake](https://snakemake.readthedocs.io/en/stable/) and consists of the following steps:

1. Obtain genome database in `fasta` and `gff` format (`python`, [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/))
   1. Using automatic download from NCBI with a `RefSeq` ID
   2. Using user-supplied files
2. Find all possible guide RNAs for the given sequence, with many options for customization (`R`, `crisprVerse`)
3. Collect on-target and off-target scores (`R`, `crisprVerse`, `Bowtie2`)
4. Filter and rank guide RNAs based on scores and return final list (`R`, `crisprVerse`)
5. Generate report with overview figures and statistics (`R markdown`)

If you want to contribute, report issues, or suggest features, please get in touch on [github](https://github.com/MPUSP/snakemake-crispr-guides).

## Installation

### Snakemake

Step 1: Install snakemake with `conda`, `mamba`, `micromamba` (or any another `conda` flavor). This step generates a new conda environment called `snakemake-crispr-guides`, which will be used for all further installations.

```
conda create -c conda-forge -c bioconda -n snakemake-crispr-guides snakemake
```

Step 2: Activate conda environment with snakemake

```
source /path/to/conda/bin/activate
conda activate snakemake-crispr-guides
```

Alternatively, install `snakemake` using pip:

```
pip install snakemake
```

Or install `snakemake` globally from linux archives:

```
sudo apt install snakemake
```

### Additional tools

**Important note:**

All other dependencies for the workflow are **automatically pulled as `conda` environments** by snakemake, when running the workflow with the `--use-conda` parameter (recommended).

In case the workflow should be executed **without automatically built `conda` environments**, the packages need to be installed manually.


**crisprVerse (Bioconductor)**

This will install R and all required dependencies of R- and Bioconductor-packages.

```bash
conda install -c bioconda bioconductor-crisprverse
```

In case a package is missing, update packages and (re-) install conventional R packages from within R:

```r
update.packages()
install.packages("BiocManager")
```

For Bioconductor packages, use for example:

```r
BiocManager::install("GenomeInfoDbData")
```

## Running the workflow

### Input data

The workflow requires the following input:

1. An NCBI Refseq ID, e.g. `GCF_000006945.2`. Find your genome assembly and ID on [NCBI genomes](https://www.ncbi.nlm.nih.gov/data-hub/genome/)
2. OR use a custom pair of `*.fasta` file and `*.gff` file that describe the genome of choice


### Execution

To run the workflow from command line, change the working directory.

```
cd /path/to/snakemake-crispr-guides
```

Adjust the global and module-specific options in the default config file `config/config.yml`.
Before running the entire workflow, you can perform a dry run using:

```
snakemake --dry-run
```

To run the complete workflow with test files, execute the following command. The definition of the number of compute cores is mandatory.

```
snakemake --cores 10 --use-conda
```

To supply a custom config file and/or use options that override the defaults, run the workflow like this:

```
snakemake --cores 10 --use-conda \
  --configfile 'config/my_config.yml' \
  --config \
  option='input'
```

### Parameters

This table lists all parameters that can be used to run the workflow.

| parameter              | type      | details                                   | default                                                           |
| ---------------------- | --------- | ----------------------------------------- | ----------------------------------------------------------------- |
| GET_GENOME             |           |                                           |                                                                   |
| database               | character | one of `ncbi`, `manual`                   | `ncbi`                                                            |
| assembly               | character | RefSeq ID                                 | `GCF_000006945.2`                                                 |
| fasta                  | path      | optional input                            | `Null`                                                            |
| gff                    | path      | optional input                            | `Null`                                                            |
| DESIGN_GUIDES          |           |                                           |                                                                   |
| target_region          | numeric   | use subset of regions for testing         | `["NC_003277.2"]`                                                 |
| tss_window             | numeric   | upstream/downstream window around TSS     | `[0, 500]`                                                        |
| circular               | logical   | is the genome circular?                   | `False`                                                           |
| canonical              | logical   | only canonical PAM sites are included     | `True`                                                            |
| both_strands           | logical   | find guides on "+"/"-" or only "+" strand | `True`                                                            |
| spacer_length          | numeric   | desired length of guides                  | `20`                                                              |
| guide_aligner          | character | one of `biostrings`, `bowtie`             | `biostrings`                                                      |
| crispr_enzyme          | character | CRISPR enzyme ID                          | `SpCas9`                                                          |
| gc_content_range       | numeric   | range of allowed GC content               | `[30, 70]`                                                        |
| score_methods          | character | see _crisprScore_ package                 | `["ruleset1", "ruleset3", "crisprater", "crisprscan", "tssdist"]` |
| score_weights          | numeric   | opt. weights when calculating mean score  | `[1, 1, 1, 1, 1]`                                                 |
| restriction_sites      | character | sequences to omit in entire guide         | `Null`                                                            |
| bad_seeds              | character | sequences to omit in seed region          | `["ACCCA", "ATACT", "TGGAA"]`                                     |
| filter_top_n           | numeric   | max number of guides to return            | `10`                                                              |
| filter_score_threshold | numeric   | mean score to use as lower limit          | `Null`                                                            |

### Off-target scores

The pipeline maps each guide RNA to the target genome and -by default- counts the number of alternative alignments with 1, 2, 3, or 4 mismatches. All guide RNAs that map to any other position including up to 4 allowed mismatches is removed.

### On-target scores

The list of available on-target scores in the [R crisprScore package](https://github.com/crisprVerse/crisprScore) is larger than the different scores included by default. It is important to note that the computation of some scores does not necessarily make sense for the design of every CRISRP library. For example, several scores were obtained from analysis of Cas9 cutting efficiency in human cell lines. For such scores it is questionable if they are uesful for the design of a very differenr library, for example a dCas9 CRISPR inhibition library in bacteria.

Another good reason to exclude some scores is the computational resources they require. Particularly deep learning-derived scores are calculated by machine learning models that require both a lot of extra resources in terms of disk space (downloaded and installed via `basilisk` and `conda` environments) and processing power (orders of magnitude longer computation time).

Users can look up all available scores on the [R crisprScore github page](https://github.com/crisprVerse/crisprScore) and decide which ones should be included. In addition, the default behavior of the pipeline is to compute an average score and select the top N guides based on it. The average score is the weighted mean of all single scores and the weights can be defined in the `config/config.yml` file. If a score should be excluded from the ranking, it's weight can simply be set to zero.

## Output

The workflow generates the following output from its modules:

<details markdown="1">
<summary>get_genome</summary>

- `genome.fasta`: Supplied or downloaded fasta file
- `genome.gff`: Supplied or downloaded gff file
- `log.txt`: Log file for this module

</details>

<details markdown="1">
<summary>design_guides</summary>

- `guideRNAs_top.csv`: Table with top N guide RNAs per gene remaining after filtering
- `log.txt`: Log file for this module

</details>

## Authors

- The custom `snakemake`, `R`, `R markdown`, and `python` scripts were written by Michael Jahn, PhD
- Affiliation: [Max-Planck-Unit for the Science of Pathogens](https://www.mpusp.mpg.de/) (MPUSP), Berlin, Germany
- Visit the MPUSP github page at https://github.com/MPUSP for info on this workflow and other projects
- Visit the author's github page at https://github.com/m-jahn for info on other projects

## References

- Essential tools are linked in the top section of this document
- The core of this workflow is the [Bioconductor package `crisprVerse`](http://bioconductor.org/packages/release/bioc/html/crisprVerse.html):

> Hoberecht, L., Perampalam, P., Lun, A. et al. _A comprehensive Bioconductor ecosystem for the design of CRISPR guide RNAs across nucleases and technologies_. Nat Commun 13, 6568 (**2022**). https://doi.org/10.1038/s41467-022-34320-7.
