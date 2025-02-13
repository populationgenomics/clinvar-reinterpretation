from datetime import datetime
from functools import cache
from os.path import join

from cpg_utils import Path, to_path
from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import get_batch

from cpg_flow.stage import MultiCohortStage, StageInput, StageOutput, stage
from cpg_flow.targets import MultiCohort


@cache
def get_output_folder():
    """
    get the folder to use for this run
    sits in the cpg-common bucket, in a folder named YY-MM
    """

    return to_path(
        join(
            config_retrieve(['storage', 'common', 'analysis']),
            'clinvarbitration',
            datetime.now().strftime('%y-%m'),  # noqa: DTZ005
        ),
    )


@stage
class CopyLatestClinvarFiles(MultiCohortStage):
    def expected_outputs(self, mc: MultiCohort) -> dict[str, Path]:
        return {
            'submission_file': get_output_folder() / 'submission_summary.txt.gz',
            'variant_file': get_output_folder() / 'variant_summary.txt.gz',
        }

    def queue_jobs(self, mc: MultiCohort, inputs: StageInput) -> StageOutput:
        """
        run a wget copy of the relevant files into GCP
        """
        bash_job = get_batch().new_bash_job('wget latest ClinVar raw files')
        bash_job.image(config_retrieve(['workflow', 'driver_image']))

        directory = 'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/'
        sub_file = 'submission_summary.txt.gz'
        var_file = 'variant_summary.txt.gz'

        bash_job.command(f'wget -q {directory}{sub_file} -O {bash_job.subs}')
        bash_job.command(f'wget -q {directory}{var_file} -O {bash_job.vars}')

        outputs = self.expected_outputs(mc)

        get_batch().write_output(bash_job.subs, str(outputs['submission_file']))
        get_batch().write_output(bash_job.vars, str(outputs['variant_file']))

        return self.make_outputs(data=outputs, jobs=bash_job, target=mc)


@stage(required_stages=[CopyLatestClinvarFiles], analysis_type='clinvarbitration', analysis_keys=['clinvar_decisions'])
class GenerateNewClinvarSummary(MultiCohortStage):
    def expected_outputs(self, mc: MultiCohort) -> dict[str, Path]:
        """
        a couple of files and a HT as Paths
        """
        return {
            'clinvar_decisions': get_output_folder() / 'clinvar_decisions.ht.tar.gz',
            'snv_vcf': get_output_folder() / 'pathogenic_snvs.vcf.bgz',
        }

    def queue_jobs(self, mc: MultiCohort, inputs: StageInput) -> StageOutput:
        # relatively RAM intensive, short running task
        job = get_batch().new_job('Run ClinvArbitration Summary')
        job.image(config_retrieve(['workflow', 'driver_image'])).memory('highmem').cpu('2')

        # declare a resource group, leave the HT path implicit
        job.declare_resource_group(output={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'})

        # get the expected outputs
        outputs = self.expected_outputs(mc)

        if sites_to_blacklist := config_retrieve(['workflow', 'site_blacklist'], []):
            blacklist_sites = ' '.join(f'"{site}"' for site in sites_to_blacklist)
            blacklist_string = f' -b {blacklist_sites}'
        else:
            blacklist_string = ''

        var_file = get_batch().read_input(inputs.as_str(mc, CopyLatestClinvarFiles, 'variant_file'))
        sub_file = get_batch().read_input(inputs.as_str(mc, CopyLatestClinvarFiles, 'submission_file'))

        # resummary is an entrypoint alias for scripts/resummarise_clinvar.py
        job.command(f'resummary -v {var_file} -s {sub_file} -o {job.output} --minimal {blacklist_string}')

        # copy out compressed
        job.command(f'tar -czf {job.compressed} {job.output}.ht')

        # selectively copy back some outputs
        get_batch().write_output(job.output, str(outputs['snv_vcf']).removesuffix('.vcf.bgz'))
        get_batch().write_output(job.compressed, str(outputs['clinvar_decisions']))

        return self.make_outputs(target=mc, data=outputs, jobs=job)


@stage(required_stages=[GenerateNewClinvarSummary])
class AnnotateClinvarSnvsWithBcftools(MultiCohortStage):
    """
    take the vcf output from the clinvar stage, and apply consequence annotations
    this uses BCFtools for speed and re-deployability, instead of VEP
    """

    def expected_outputs(self, mc: MultiCohort) -> Path:
        return get_output_folder() / 'annotated_snv.tsv'

    def queue_jobs(self, mc: MultiCohort, inputs: StageInput) -> StageOutput:
        outputs = self.expected_outputs(mc)

        # need a genome reference - mandatory argument
        ref_fasta = get_batch().read_input(config_retrieve(['workflow', 'ref_fasta']))
        snv_vcf = get_batch().read_input(inputs.as_str(mc, GenerateNewClinvarSummary, 'snv_vcf'))

        genome_build = config_retrieve(['workflow', 'genome_build'], 'GRCh38')
        if genome_build not in ['GRCh37', 'GRCh38']:
            raise ValueError('Unsupported genome build specified')

        gff3 = f'/clinvarbitration/bcftools_data/{genome_build}.gff3.gz'

        job = get_batch().new_job('Annotate Clinvar SNVs with BCFtools', attributes={'tool': 'bcftools'})
        job.image(config_retrieve(['workflow', 'driver_image']))
        job.storage('10G')

        # -g is the GFF3 file, -f is the reference fasta
        # TODO check if this works correctly
        #  - it exists on a base mismatch between... VCF and reference? GFF3 and reference?
        # -e 'CHROM=="chrY"' is a partial solve, it just leaves the variants unannotated
        job.command(f"""bcftools csq -f {ref_fasta} -e 'CHROM=="chrY"' -g {gff3} {snv_vcf} -o annotated_output.vcf""")

        # split the bcftools CSQ fields, filter to missense, and write out a tab-delimited file
        # -d - duplicate, writes each transcript onto a new line
        # -s :missense - only keep Consequence==missense variants
        # -f - format string - tab delimited, Transcript, Amino Acid Change, ClinVar allele ID, ClinVar gold stars
        job.command(
            f'bcftools +split-vep annotated_output.vcf -d -s :missense -f "%transcript\t%amino_acid_change\t%allele_id\t%gold_stars\n" > {job.output}',  # noqa: E501
        )

        get_batch().write_output(job.output, str(outputs))

        return self.make_outputs(target=mc, jobs=job, data=outputs)


@stage(required_stages=[AnnotateClinvarSnvsWithBcftools])
class Pm5TableGeneration(MultiCohortStage):
    def expected_outputs(self, mc: MultiCohort) -> dict[str, Path]:
        """
        a single HT, compressed into a tarball
        """
        return {
            'pm5_ht': get_output_folder() / 'clinvar_pm5.ht.tar.gz',
            'pm5_json': get_output_folder() / 'clinvar_pm5.json',
        }

    def queue_jobs(self, mc: MultiCohort, inputs: StageInput) -> StageOutput:
        # declare a resource group, but don't copy the whole thing back
        job = get_batch().new_job('Run ClinvArbitration PM5')
        job.image(config_retrieve(['workflow', 'driver_image']))

        # get the expected outputs
        outputs = self.expected_outputs(mc)

        annotated_snvs = get_batch().read_input(inputs.as_str(mc, AnnotateClinvarSnvsWithBcftools))

        job.command(f'pm5_table -i {annotated_snvs} -o output')
        job.command(f'mv output.json {job.json}')
        get_batch().write_output(job.json, str(outputs['pm5_json']))

        # compress the HT and remove as a single file
        job.command(f'tar -czf {job.ht} output.ht')
        get_batch().write_output(job.ht, str(outputs['pm5_ht']))

        return self.make_outputs(target=mc, data=outputs, jobs=job)


def populate_job_meta(output_file: str):
    """
    populates some metadata for the job

    Args:
        output_file (str): path to the output file

    Returns:
        dict of metadata
    """

    from clinvarbitration import __version__ as clinvarbitration_version

    print(f'Generating output meta for {output_file}')
    return {'image': config_retrieve(['workflow', 'driver_image']), 'version': clinvarbitration_version}


@stage(
    required_stages=[
        GenerateNewClinvarSummary,
        AnnotateClinvarSnvsWithBcftools,
        Pm5TableGeneration,
    ],
    analysis_type='clinvarbitration',
    update_analysis_meta=populate_job_meta,
)
class PackageForRelease(MultiCohortStage):
    """
    takes the data created so far, and packages it up for release
    """

    def expected_outputs(self, multicohort: MultiCohort) -> Path:
        return get_output_folder() / 'clinvarbitration.tar.gz'

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput:
        """
        Localise all the previously generated data into a folder
        tarball it, and write out as a single file
        """
        tar_output = self.expected_outputs(multicohort)

        job = get_batch().new_job('Package ClinvArbitration data for release')
        job.storage('10G')
        job.image(config_retrieve(['workflow', 'driver_image']))

        # find paths to the previous outputs
        pm5 = inputs.as_dict(multicohort, Pm5TableGeneration)
        clinvar_decisions = inputs.as_str(multicohort, GenerateNewClinvarSummary, 'clinvar_decisions')

        annotated_tsv = get_batch().read_input(inputs.as_str(multicohort, AnnotateClinvarSnvsWithBcftools))
        pm5_ht = get_batch().read_input(str(pm5['pm5_ht']))
        pm5_json = get_batch().read_input(str(pm5['pm5_json']))
        decisions_ht = get_batch().read_input(clinvar_decisions)

        # create a new folder, move the files into it
        # unpack all the already compressed HTs into it
        # the compress once, containing all files for distribution
        job.command(
            f"""
            mkdir clinvarbitration_data
            mv {annotated_tsv} clinvarbitration_data/annotated_snvs.tsv
            mv {pm5_json} clinvarbitration_data/pm5.json
            tar -xzf {pm5_ht} -C clinvarbitration_data
            tar -xzf {decisions_ht} -C clinvarbitration_data
            tar -czf {job.output} clinvarbitration_data
        """,
        )
        get_batch().write_output(job.output, str(tar_output))
        return self.make_outputs(multicohort, data=tar_output, jobs=job)
