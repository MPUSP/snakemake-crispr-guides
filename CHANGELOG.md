# Changelog

## 1.0.0 (2023-07-13)


### Features

* add major module to find, score and filter guide RNAs ([4f93545](https://github.com/MPUSP/snakemake-crispr-guides/commit/4f935454f1496fef2872fe26771fcf959c4fe5ff))
* added basic R markdown report module ([3950bcd](https://github.com/MPUSP/snakemake-crispr-guides/commit/3950bcde6a23ff600ef0c60d2301567087e2c6cf))
* added CI test dir ([bdfe96c](https://github.com/MPUSP/snakemake-crispr-guides/commit/bdfe96c8400cf13edd80e91b1888328d8cc7db38))
* added dev branch to GH actions ([303de61](https://github.com/MPUSP/snakemake-crispr-guides/commit/303de61ca37ad63d3a39de119096df10206ea94f))
* added function to exclude RNAs from guide design ([63f0b0e](https://github.com/MPUSP/snakemake-crispr-guides/commit/63f0b0ec2ce995c8d23ebf2489516cae63685425))
* added functionality to design random control guides ([a32185b](https://github.com/MPUSP/snakemake-crispr-guides/commit/a32185ba426b6de48c4c5c249620ef1e018a1f8f))
* added functionality to handle and visualize guides with multiple targets. Closes [#11](https://github.com/MPUSP/snakemake-crispr-guides/issues/11) ([ac4cd79](https://github.com/MPUSP/snakemake-crispr-guides/commit/ac4cd79088a461be873be80d59ddba017115ce5c))
* added log entry with number of non-targeted genes ([1b1653f](https://github.com/MPUSP/snakemake-crispr-guides/commit/1b1653f2c464258823851487ec32b7b0b13e476c))
* added module to convert html report to pdf ([37f30a4](https://github.com/MPUSP/snakemake-crispr-guides/commit/37f30a4a35a55875d9bef04b308b7f0c111fc359))
* added module to generate BSgenome object on the fly, for use with bowtie, closes [#4](https://github.com/MPUSP/snakemake-crispr-guides/issues/4) ([c344d0d](https://github.com/MPUSP/snakemake-crispr-guides/commit/c344d0d70ea42b147b76df41d8dbe99c91526ab1))
* added module to make bowtie2 index ([cf11e73](https://github.com/MPUSP/snakemake-crispr-guides/commit/cf11e73821fc56a4d30f12e65f43d521686458a1))
* added more TSS features to result table, and new viz figures ([8ffa6a2](https://github.com/MPUSP/snakemake-crispr-guides/commit/8ffa6a24fa74ba593a40f36222850a7e0f994b48))
* added option to add arbitrary linkers to guide RNA ends ([d931cb7](https://github.com/MPUSP/snakemake-crispr-guides/commit/d931cb722c00cb7e75f0a673f46a400cf69dbbb7))
* added optional weights to calculate mean score ([5678c38](https://github.com/MPUSP/snakemake-crispr-guides/commit/5678c38689dcd374aa483074f2f191d3dcea22e2))
* added output of failed guide RNA design to report ([4019861](https://github.com/MPUSP/snakemake-crispr-guides/commit/4019861fecd5a3b5b9ea9dccb953261ce6a7e109))
* added percentages to report plots ([da8b5d9](https://github.com/MPUSP/snakemake-crispr-guides/commit/da8b5d9fdb8ef8e8b5185297bece70b6e246881b))
* added score for G eenrichment, closes [#5](https://github.com/MPUSP/snakemake-crispr-guides/issues/5) ([a1b35a1](https://github.com/MPUSP/snakemake-crispr-guides/commit/a1b35a136ffc79fcf2c979e8ec1f784d00602e5d))
* changed license to GPL ([9ce508a](https://github.com/MPUSP/snakemake-crispr-guides/commit/9ce508a30f8b9ba3ee00fcd5d36a298fe9ab58d2))
* export table of targets where guide design failed ([e377f4b](https://github.com/MPUSP/snakemake-crispr-guides/commit/e377f4b8a29ca88803c216c3d8a57fbf0d34b997))
* filter for guides that fall within TSS window _and_ transcript region + visualization ([39e9cea](https://github.com/MPUSP/snakemake-crispr-guides/commit/39e9ceadebb3211c879a2c7ccfeac664ea3befb3))
* finished module for obtaining genomes from NCBI or user ([5c1fa2b](https://github.com/MPUSP/snakemake-crispr-guides/commit/5c1fa2b97c92e9aa362d05a33d117f686ae50e54))
* first draft of module to generate bsgenome object on the fly ([a5f8b77](https://github.com/MPUSP/snakemake-crispr-guides/commit/a5f8b77633d6d0ff0dc4751563a2eb65e614eb6e))
* fixed small bugs to improve usage of custom fasta/gff files ([ab76a25](https://github.com/MPUSP/snakemake-crispr-guides/commit/ab76a25b153d70415faae4113e8bf1044a3662af))
* improved prediction around TSS and added TSS dist score, closes [#3](https://github.com/MPUSP/snakemake-crispr-guides/issues/3) ([290b675](https://github.com/MPUSP/snakemake-crispr-guides/commit/290b67569d1702da1068bca6cdeea5396b95f98c))
* included visualization of guide RNA location, and option for strand preference ([d345099](https://github.com/MPUSP/snakemake-crispr-guides/commit/d3450996b001d579078bf10885adad5a41aa9705))
* simplified handling of result files ([ddcb925](https://github.com/MPUSP/snakemake-crispr-guides/commit/ddcb925fe87a4ab2d1dd6f9c911543a14d44b44f))


### Bug Fixes

* added missing documentation for paramter filter_multi_targets ([32476e5](https://github.com/MPUSP/snakemake-crispr-guides/commit/32476e57be7fba97112bd29aa946adfff8520395))
* added r-base environment to conda list ([e2d7cbf](https://github.com/MPUSP/snakemake-crispr-guides/commit/e2d7cbfa92e6e6d884696afec41addf06fc88140))
* added variable TSS windows, closes [#9](https://github.com/MPUSP/snakemake-crispr-guides/issues/9) ([014e6ec](https://github.com/MPUSP/snakemake-crispr-guides/commit/014e6ec7ca24071b415ad1f8412d081ebfe0a447))
* changed syntax for env defs, use dry run for test ([32f5e82](https://github.com/MPUSP/snakemake-crispr-guides/commit/32f5e82525e514e455f66aef948b26e63d405db2))
* conda env error nothing provides libgcc-ng ([ffc9cec](https://github.com/MPUSP/snakemake-crispr-guides/commit/ffc9cec6129459b40b8672a05bd48f1f26a9f0e6))
* fixed small syntax error for env ([bd1775b](https://github.com/MPUSP/snakemake-crispr-guides/commit/bd1775baec6880a7ee31b371a0612d7b853b1ba9))
* improvements to score correlation plot ([7afc9cb](https://github.com/MPUSP/snakemake-crispr-guides/commit/7afc9cb488f3b053766b7739ee4d69dfeef71927))
* removed typo in dry snakemake test command ([8e3e5c0](https://github.com/MPUSP/snakemake-crispr-guides/commit/8e3e5c0b94769a139710ad11171a7448d440d67b))
* replaced bowtie2 with bowtie1 as a req for crisprDesign module ([f8cd5fb](https://github.com/MPUSP/snakemake-crispr-guides/commit/f8cd5fbd985ab2a4ed80c73c66fdbf3157c20567))
* small improvements to report and docs ([d24c36c](https://github.com/MPUSP/snakemake-crispr-guides/commit/d24c36c7889efc4dba083a75c8ea859525d40a29))
* small update to guide RNA drawing ([18eac19](https://github.com/MPUSP/snakemake-crispr-guides/commit/18eac19aa440b3c8b6bf89ac9bf018d2a892c16f))
* some corrections to pipeline overview ([2e519ed](https://github.com/MPUSP/snakemake-crispr-guides/commit/2e519ed8a1b7b98b4b45cbdef4890128dbeb0ae2))