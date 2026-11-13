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


def mem_type(value):
    """Accept 'auto' (use the job's memory limit) or a positive number of GB."""
    value = str(value).strip().lower()
    if value in ('auto', '0', ''):
        return 'auto'
    try:
        n = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError('must be a number of GB or "auto", got "{}"'.format(value))
    if n <= 0:
        raise argparse.ArgumentTypeError('must be a positive number of GB or "auto", got "{}"'.format(value))
    return str(n)


def prefix_type(value):
    """Validate a user-supplied filename prefix: letters, digits, underscores; case kept."""
    value = value.strip()
    if not re.fullmatch(r'[A-Za-z0-9_]+', value):
        raise argparse.ArgumentTypeError(
            'must contain only letters, digits and underscores, got "{}"'.format(value))
    return value


def sample_type(value):
    """Validate a sample name given with --sample (as detected from -d, without read suffix)."""
    value = value.strip()
    if not value or '/' in value:
        raise argparse.ArgumentTypeError(
            'must be a sample name as detected in -d (no path), got "{}"'.format(value))
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

    # Local mode: one output set per sample instead of one merged set
    per_sample = parser.add_argument_group(
        'per-sample outputs (local mode only)',
        'By default all samples in -d are merged into one BAM. With these options each sample '
        '(e.g. each tissue) gets its own BAM, CSI, bigWig and BED in <outdir>/<name>/<sample>/.')
    per_sample.add_argument('--per_sample', dest='per_sample', default=False, action='store_true',
                            help='align each sample in -d separately and write a full output set for each')
    per_sample.add_argument('--sample', dest='samples', action='append', default=None,
                            type=sample_type, metavar='NAME',
                            help='process only this sample from -d (repeatable, e.g. '
                                 '--sample Zacosmia_maculata_Head). Implies --per_sample')

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
    parser.add_argument('-M', '--mem_gb', dest='mem_gb', type=mem_type, default='auto',
                        help="Memory available to the pipeline in GB, or 'auto' to use the job's "
                             "Slurm/cgroup memory limit. Java heaps (BBTools) are capped at 80%% of this. "
                             "Values above the limit are capped. Default: %(default)s")

    # Output naming: at most one of these may be given
    naming = parser.add_argument_group(
        'output file naming',
        'Use at most one of these options. Default: Genus_species[_subspecies] from the '
        'scientific name. Except with -R, the prefix is followed by '
        '_<assembly>_RNA-Seq-alignments_<date>. With --per_sample, the sample label '
        '(e.g. Head) is added after the prefix.')
    name_opts = naming.add_mutually_exclusive_group()
    name_opts.add_argument('-R', '--Run_prefix', dest='Run_prefix', default=False, action='store_true',
                           help='use the run/sample name as the whole output name '
                                '(the first run, or each sample with --per_sample)')
    name_opts.add_argument('--prefix', dest='prefix', type=prefix_type, default=None,
                           help='file name prefix for outputs, e.g. Apis_mellifera_ligustica')
    name_opts.add_argument('--input_prefix', dest='input_prefix', default=False, action='store_true',
                           help='derive the prefix from the input file names: the common leading '
                                'part of the sample names with -d, or the TSV file name with -i')

    args = parser.parse_args(argv)

    if bool(args.input) == bool(args.fastq_dir):
        parser.error('specify exactly one of -i/--input or -d/--fastq_dir')
    if args.fastq_dir and (not args.species or len(args.species.split()) < 2):
        parser.error('-d/--fastq_dir requires -s/--species "Genus species"')
    if not args.genome:
        parser.error('-g/--genome is required')

    # --sample implies --per_sample; both need local FASTQ input
    if args.samples:
        args.samples = list(dict.fromkeys(args.samples))   # drop duplicates, keep order
        args.per_sample = True
    if args.per_sample and not args.fastq_dir:
        parser.error('--per_sample/--sample require -d/--fastq_dir')

    # -R / --prefix / --input_prefix conflicts are rejected by the mutually exclusive group

    return args
