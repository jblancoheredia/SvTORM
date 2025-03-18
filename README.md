# mskcc/svtorm

[![GitHub Actions CI Status](https://github.com/mskcc/svtorm/actions/workflows/ci.yml/badge.svg)](https://github.com/mskcc/svtorm/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/mskcc/svtorm/actions/workflows/linting.yml/badge.svg)](https://github.com/mskcc/svtorm/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/mskcc/svtorm)

## Pipeline Overview

![Pipeline Steps](assets/SVtorm.svg)

## Introduction

**mskcc/svtorm** is a bioinformatics pipeline that calls structural variants from target sequenced data.

![Pipeline Steps](assets/SVtorm.svg)

0. Starts with a pre-made BAM file.
1. Calling SVs
   - ([`Delly`](https://github.com/dellytools/delly))
   - ([`Svaba`](https://github.com/walaj/svaba))
   - ([`Manta`](https://github.com/Illumina/manta))
   - ([`Gridss`](https://github.com/PapenfussLab/gridss))
2. Merging Calls ([`SURVIVOR`](https://github.com/fritzsedlazeck/SURVIVOR))
3. Bed to Interval list ([`GATK`](https://gatk.broadinstitute.org/hc/en-us/articles/360035531852-Intervals-and-interval-lists))
4. ReCalling ([`Gridss`](https://github.com/PapenfussLab/gridss))
5. Filtering Calls ([`SURVIVOR`](https://github.com/fritzsedlazeck/SURVIVOR))
6. Annotate SVs ([`iAnnotateSV`](https://github.com/mskcc/iAnnotateSV))
7. Draw SVs ([`DrawSV`](https://github.com/jblancoheredia/DrawSV))
8. Check for expected SVs in Controls 
9. Present QC for raw reads ([`MultiQC`](http://multiqc.info/)) 

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

<!-- TODO nf-core: Describe the minimum required steps to execute the pipeline, e.g. how to prepare samplesheets.
     Explain what rows and columns represent. For instance (please edit as appropriate):

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
```

Each row represents a fastq file (single-end) or a pair of fastq files (paired end).

-->

Now, you can run the pipeline using:

<!-- TODO nf-core: update the following command to include all required parameters for a minimal example -->

```bash
nextflow run mskcc/svtorm \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

## Credits

mskcc/svtorm was originally written by blancoj@mskcc.org.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use mskcc/svtorm for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).