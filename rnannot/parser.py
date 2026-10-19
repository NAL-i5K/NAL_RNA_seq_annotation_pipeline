import argparse
import datetime
import re
import sys


def threads_type(value):
    """Accept 'auto' (use all CPUs allocated to the job) or a positive integer."""
    value = str(value).strip().lower()
    if value in ('auto', '0', ''):
        return 'auto'
    try:
        n = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError('must be a positive integer or "auto", got "{}"'.format(value))
    if n < 1:
        raise argparse.ArgumentTypeError('must be a positive integer or "auto", got "{}"'.format(value))
    return str(n)


def gggsss_type(value):
    """Validate a user-supplied organism code: letters/digits only, stored lowercase."""
    value = value.strip().lower()
    if not re.fullmatch(r'[a-z0-9]+', value):
        raise argparse.ArgumentTypeError(
            'must contain only letters and digits, got "{}"'.format(value))
    return value


def parse_args(argv):
    parser = argparse.ArgumentParser(description='Easy to use pipeline built for large-scale RNA-seq mapping with a genome assembly')
    # Input: either an SRA metadata tsv (-i) or a directory of local FASTQ files (-d)
    parser.add_argument('-i', '--input', dest='input', type=str, default=None,
                        help="A tsv file with a list of SRA runs' information (from an NCBI SRA metadata pull).")
    parser.add_argument('-d', '--fastq_dir', dest='fastq_dir', default=None,
                        help='Directory of local FASTQ files (.fastq/.fq, optionally .gz). Use instead of -i.')
    parser.add_argument('-s', '--species', dest='species', default=None,
                        help='Scientific name, e.g. "Apis mellifera" or "Apis mellifera ligustica". Required with -d.')
    parser.add_argument('--platform', dest='platform', default='ILLUMINA',
                        help='Sequencing platform for local files. Default: %(default)s')
    parser.add_argument('--model', dest='model', default='',
                        help='Instrument model for local files, e.g. "Illumina HiSeq 2500" or "Element AVITI". '
                             'Selects the adapter set; if empty, adapters are detected with BBMerge (paired) '
                             'or taken from BBMap (single).')
    parser.add_argument('--gggsss', dest='gggsss', type=gggsss_type, default=None,
                        help='Organism code for output file names. Default: derived from the scientific name '
                             '(gggsss for genus+species, ggssii when a subspecies is present).')
    parser.add_argument('-g', '--genome', dest='genome', help='A fasta file to align with.')
    parser.add_argument('-n', '--name', nargs='?',
                        default=datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S"),
                        help='name of the output folder, if not specified, use the time of start')
    parser.add_argument('-o', '--outdir', dest='outdir', nargs='?', default='.',
                        help='directory of output folder at, if not specified, use current folder')
    parser.add_argument('-a', '--assembly', dest='assembly', default='AssemblyName', help='The assembly name is used for naming output file.')
    parser.add_argument('-t', '--tempFile', dest='tempFile', default=False, action='store_true', help='if specified, intermediate output bam files will be kept')
    parser.add_argument('-m', '--MaximumSRA', dest='MaximumSRA', type=int, default=10,
                        help='The maximum amount of the sra files downloaded from NCBI (SRA mode only). The default is 10')
    parser.add_argument('-T', '--threads_num', dest='threads_num', type=threads_type, default='auto',
                        help="Number of threads to use, or 'auto' to use all CPUs allocated to the job "
                             "(e.g. by Slurm --cpus-per-task). Values above the allocation are capped. "
                             "Default: %(default)s")
    parser.add_argument('-R', '--Run_prefix', dest='Run_prefix', default=False, action='store_true', help='if specified, use Run number as output prefix')
    args = parser.parse_args(argv)
    if bool(args.input) == bool(args.fastq_dir):
        parser.error('specify exactly one of -i/--input or -d/--fastq_dir')
    if args.fastq_dir and (not args.species or len(args.species.split()) < 2):
        parser.error('-d/--fastq_dir requires -s/--species "Genus species"')
    if not args.genome:
        parser.error('-g/--genome is required')
    return args
