# ClinVar, reimagined

## motivation

During the creation of [Talos](https://www.github.com/populationgenomics/automated-interpretation-pipeline), a tool for identifying clinically relevant variants in large cohorts, we use [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) ratings as a contributing factor in determining pathogenicity. During development of this tool we determined that the default summaries generated in ClinVar were highly conservative; see the [table here](https://www.ncbi.nlm.nih.gov/clinvar/docs/clinsig/#agg_germline) describing the aggregate classification logic.

## content

This repository contains an alternative algorithm ([described here](docs/algorithm.md)) for re-aggregating the individual ClinVar submissions, generating decisions which favour clear assignment of pathogenic/benign ratings instead of defaulting to 'conflicting'. These ratings are not intended as a replacement of ClinVar's own decisions, but may provide value by showing that that though conflicting submissions exist, there is a clear bias towards either benign or pathogenic ratings.

We aim to re-run this process monthly, and publish the resulting files as Releases. You can also run this locally (see `Usage`)

## outputs

* JSON file & Hail Table of all revised decisions
* Sites-only VCF of all revised decisions; this can be used as a custom annotation source in VEP or other annotators
* Sites-only VCF of all Pathogenic-SNVs, for feeding into the ACMG criteria PM5 analysis
* Annotated VCF of Pathogenic SNVs using `bcftools csq`
* JSON file & Hail Table of all Pathogenic missense changes, indexed on Transcript and Codon. This is usable as a PM5 annotation resource.


# Usage

## Download

Before cloning this repository, install [git-lfs](https://github.com/git-lfs/git-lfs) to automatically download the GFF3 files used by BCFTools. These are described in more detail in the relevant [README](bcftools_data/README.txt).

If git-lfs is not available, the relevant GFF3 files can be easily downloaded from the Ensembl FTP site as  documented in the README.

## Non-CPG usage

Example scripts outline the steps needed to run this locally. There are two versions, suffixed with `_docker` for running inside a container, or `_local` if you have the package and bcftools installed locally.

1. [optional] Build a Docker image from this repository
2. Create a local directory called `data`
3. Copy a reference genome for GRCh37 or 38 into `data`, call it `ref.fa`
4. If running in a container, use `docker run -v "$(pwd)/data":/data <image_name> bash example_script_docker.sh`
5. If running natively, `bash example_script_local.sh`

All outputs will be generated in the `/data` directory

## CPG-Flow

Internally at CPG, this workflow is run using [CPG-Flow](https://github.com/populationgenomics/cpg-flow), an in-house Hail Batch based workflow executor. The folder `src/clinvarbitration/cpg_flow` contains elements relating to that workflow:

* an [example config file](src/clinvarbitration/config_template.toml), with enough entries populated that a standard CPG user could dry-run the workflow locally
* a [workflow runner script](src/clinvarbitration/run_workflow.py)
* a definition of all [workflow stages](src/clinvarbitration/stages.py)

The intention is that once the Dockerfile within this repository is used, this workflow can be triggered like so:

```bash
analysis-runner \
    --skip-repo-checkout \
    --image <URI of the docker image> \
    --config <path to a config file> \
    --dataset seqr \
    --description 'resummarise_clinvar' \
    -o resummarise_clinvar \
    --access-level standard \
    run_workflow
```

A config file is required containing a few entries, some relating to this workflow specifically, some relating to cpg-flow setup:

* `workflow.driver_image`: populated by analysis-runner, points to _this_ docker image
* `site_blacklist`: list of ClinVar submitters to ignore. Useful in removing noise, or blinding to _self_ submissions
* `ref_fasta`: required to run bcftools csq. Must match the `genome_build`
* `genome_build`: used to decide whether ClinVar/Annotation is sourced using GRCh37 or GRCh38 (default)


## acknowledgements

* [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar), for providing the data which this process is based on
