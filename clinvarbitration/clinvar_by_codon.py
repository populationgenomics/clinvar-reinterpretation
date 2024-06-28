"""
Method file for re-sorting clinvar annotations by codon

Takes a VCF of annotated Pathogenic Clinvar Variants
re-indexes the data to be queryable on Transcript and Codon
writes the resulting Hail Table to the specified path

Data as input for this script should be a VCF, annotated by VEP 110
Compatibility with other versions of VEP is not guaranteed

This makes the assumption that the annotated data here
has been generated by summarise_clinvar_entries.py:

- SNV only
- Clinvar Pathogenic only
- ClinVar decision/alleles/gold stars are in INFO
"""

import json
import logging
from argparse import ArgumentParser
from collections import defaultdict

import hail as hl
from cyvcf2 import VCF, Variant


def pull_vep_from_header(vcf: VCF) -> list[str]:
    """
    yank the CSQ line out of the VCF header
    """
    for element in vcf.header_iter():
        if element['HeaderType'] == 'INFO' and element['ID'] == 'CSQ':
            return list(entry.lower() for entry in element['Description'].split('Format: ')[-1].rstrip('"').split('|'))
    raise IndexError('CSQ element not found in header')


def variant_consequences(variant: Variant, csq_header: list[str]) -> list[dict[str, str]]:
    """
    extracts the consequences for each transcript in this variant

    Args:
        variant (Variant):
        csq_header ():

    Returns:
        a list of all CSQ entries, cast as a dict
    """

    consequences: list[dict[str, str]] = []
    for csq in variant.INFO['CSQ'].split(','):
        csq_dict = dict(zip(csq_header, csq.split('|'), strict=True))
        if 'missense_variant' in csq_dict['consequence']:
            consequences.append(csq_dict)
    return consequences


def cli_main():
    """
    alternative access point with CLI arguments
    """
    logging.basicConfig(level=logging.INFO)
    parser = ArgumentParser()
    parser.add_argument('-i', help='Path to the annotated VCF')
    parser.add_argument('-o', help='Root to export PM5 table and JSON to')
    args = parser.parse_args()
    main(input_vcf=args.i, output_root=args.o)


def main(input_vcf: str, output_root: str):
    """

    Args:
        input_vcf (str): path to an input vcf
        output_root ():
    """

    # crack open a cold VCF, and have a sip
    vcf_reader = VCF(input_vcf)

    # find the header encoding all the VEP fields
    header_csq = pull_vep_from_header(vcf_reader)

    clinvar_dict = defaultdict(set)

    # iterate over the variants
    for variant in vcf_reader:
        # extract the clinvar details (added in previous script)
        clinvar_allele = variant.INFO['allele_id']
        clinvar_stars = variant.INFO['gold_stars']
        clinvar_key = f'{clinvar_allele}:{clinvar_stars}'

        # iterate over all missense consequences
        for csq_dict in variant_consequences(variant, header_csq):
            # add this clinvar entry in relation to the protein consequence
            protein_key = f"{csq_dict['ensp']}:{csq_dict['protein_position']}"
            clinvar_dict[protein_key].add(clinvar_key)

    # save the dictionary locally
    json_out_path = f'{output_root}.json'
    with open(json_out_path, 'w') as f:
        for key, value in clinvar_dict.items():
            new_dict = {'newkey': key, 'clinvar_alleles': '+'.join(value)}
            f.write(f'{json.dumps(new_dict)}\n')

    logging.info(f'JSON written to {json_out_path}')

    # now set a schema to read that into a table... if you want hail
    schema = hl.dtype('struct{newkey:str,clinvar_alleles:str}')

    # import the table, and transmute to top-level attributes
    ht = hl.import_table(json_out_path, no_header=True, types={'f0': schema})
    ht = ht.transmute(**ht.f0)
    ht = ht.key_by(ht.newkey)

    # write out
    ht.write(f'{output_root}.ht', overwrite=True)
    logging.info(f'Hail Table written to {output_root}.ht')


if __name__ == '__main__':
    cli_main()
