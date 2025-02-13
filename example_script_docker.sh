#!/usr/bin/env bash

set -ex

# run these commands in advance of running this script:

# create a docker image from this repository
#docker build --platform linux/arm64/v8 -t clinvarbitration:example .

# copy a reference genome for your intended assembly into data
# for the purposes of this script I'm assuming it is present in data and called ref.fa

# run this script inside the docker image
#docker run -v "$(pwd)/data":/data clinvarbitration:example bash example_script_docker.sh

wget -O /data/variant_summary.txt.gz https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
wget -O /data/submission_summary.txt.gz https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/submission_summary.txt.gz

# re-summarise the downloaded clinvar data
# use the `--assembly GRCh37` flag if the data is from GRCh37, defaults to GRCh38
resummary -v "/data/variant_summary.txt.gz" -s "/data/submission_summary.txt.gz" -o "/data/clinvar_summary" --minimal

# upon completion, this will have generated files in the data directory:
# - data/clinvar_summary.json - a JSON file containing the summarised data entries, one json object per line
# - data/clinvar_summary_unfiltered.vcf.bgz - a bgzipped VCF which can be used in VEP annotation
# - data/clinvar_summary.vcf.bgz - a bgzipped file containing the pathogenic SNV entries in VCF format
# - data/clinvar_summary.ht - a Hail Table containing the summarised data entries

# annotate this data with bcftools
# You'll need to provide an appropriate reference genome fasta (see note above, line 23)
bcftools csq -f /data/ref.fa -e 'CHROM=="chrY"' -g bcftools_data/GRCh38.gff3.gz /data/clinvar_summary.vcf.bgz -o /data/annotated_output.vcf
bcftools +split-vep /data/annotated_output.vcf -d -s :missense -f "%transcript\t%amino_acid_change\t%allele_id\t%gold_stars\n" > /data/annotated_snv.tsv

# use the consequence annotations to generate a PM5 table for use in annotation
pm5_table -i /data/annotated_snv.tsv -o /data/output

# upon completion, this will generate files in the data directory:
# - data/pm5.json - a JSON file containing the PM5 results, one JSON object per line
# - data/pm5.ht - a Hail table containing the same PM5 results
