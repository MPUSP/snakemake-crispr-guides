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

This workflow is a best-practice workflow for the automated generation of guide RNAs for CRISPR applications. It's main purpose is to provide a simple, efficient and easy-to-use framework to design thousands of guides simultaneously for CRISPR libraries from as little input as an organism's name/genome ID. For the manual design of single guides, users are instead referred to the even simpler web resources, such as [Chop-Chop](http://chopchop.cbu.uib.no/), [CRISPick](https://portals.broadinstitute.org/gppx/crispick/public), or [Cas-OFFinder/Cas-Designer](http://www.rgenome.net/cas-designer/).

This workflow relies to a large degree on the underlying [Bioconductor package ecosystem `crisprVerse`](http://bioconductor.org/packages/release/bioc/html/crisprVerse.html), published in 2022 by:

> Hoberecht, L., Perampalam, P., Lun, A. et al. _A comprehensive Bioconductor ecosystem for the design of CRISPR guide RNAs across nucleases and technologies_. Nat Commun 13, 6568 (**2022**). https://doi.org/10.1038/s41467-022-34320-7.

The workflow is built using [snakemake](https://snakemake.readthedocs.io/en/stable/) and consists of the following steps:

1. TODO: Obtain genome database from NCBI or use user-supplied fasta file (`python`, [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/))
2. Step 2
3. Step 3

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

The workflow requires the following input files:

1. a genome
2. an (organism) database in `*.fasta` format _OR_ a NCBI Refseq ID. Decoys (`rev_` prefix) will be added if necessary
3.

### Execution

To run the workflow from command line, change the working directory.

```
cd /path/to/snakemake-crispr-guides
```

Adjust the global and module-specific options in the config file under `config/config.yml`.
Before running the entire workflow, you can perform a dry run using:

```
snakemake --configfile 'config/config.yml' --dry-run
```

To run the complete workflow with test files, execute the following command. The definition of the number of compute cores is mandatory.

```
snakemake --configfile 'config/config.yml' --cores 10 --use-conda
```

To supply options that override the defaults, run the workflow like this:

```
snakemake --cores 10 --use-conda \
  --configfile 'config/config.yml' \
  --config \
  option='input'
```

### Parameters

This table lists all parameters that can be used to run the workflow.

| parameter | type                     | details         | example                                         |
| --------- | ------------------------ | --------------- | ----------------------------------------------- |
| genome    | `*.genbank` OR refseq ID | plain text      | `test/input/genome/test.gbk`, `GCF_000009045.1` |
| output    | path                     | valid directory | `test/output/`                                  |

## Output

The workflow generates the following output from its modules:

<details markdown="1">
<summary>output</summary>

- `file.txt`: An example file
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
