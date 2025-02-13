"""
point of entry script to run this pipeline
"""

from argparse import ArgumentParser

from clinvarbitration.stages import (
    CopyLatestClinvarFiles,
    GenerateNewClinvarSummary,
    AnnotateClinvarSnvsWithBcftools,
    Pm5TableGeneration,
    PackageForRelease,
)
from cpg_flow.workflow import run_workflow


def cli_main():
    parser = ArgumentParser()
    parser.add_argument('--dry-run', action='store_true', help='Print the commands that would be run')
    args = parser.parse_args()
    main(dry_run=args.dry_run)


def main(dry_run: bool = False):
    run_workflow(
        stages=[
            CopyLatestClinvarFiles,
            GenerateNewClinvarSummary,
            AnnotateClinvarSnvsWithBcftools,
            Pm5TableGeneration,
            PackageForRelease,
        ],
        dry_run=dry_run,
    )


if __name__ == '__main__':
    cli_main()
