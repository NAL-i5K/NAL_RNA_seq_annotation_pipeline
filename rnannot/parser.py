import argparse
import datetime
import sys


def threads_type(value):
    """Accept 'auto' (use all CPUs allocated to the job) or a positive integer."""
    value = str(value).strip().lower()
    if value in ('auto', '0', ''):
        return 'auto'
    n = int(value)  # raises ValueError for non-numbers -> argparse prints a clean error
    if n < 1:
        raise ValueError('threads must be a positive integer or "auto"')
    return str(n)


def parse_args(argv):
    parser = argparse.ArgumentParser(description='Easy to use pipeline built for large-scale RNA-seq mapping with a genome assembly')
    parser.add_argument('-i', '--input', dest='input', type=str, help="A tsv file with a list of SRA runs' information.")
    parser.add_argument('-g', '--genome', dest='genome', help='A fasta file to align with.')
    parser.add_argument('-n', '--name', nargs='?',
                        default=datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S"),
                        help='name of the output folder, if not specified, use the time of start')
    parser.add_argument('-o', '--outdir', dest='outdir', nargs='?', default='.',
                        help='directory of output folder at, if not specified, use current folder')
    parser.add_argument('-a', '--assembly', dest='assembly', default='AssemblyName', help='The assembly name is used for naming output file.')
    parser.add_argument('-t', '--tempFile', dest='tempFile', default=False, action='store_true', help='if specified, intermediate output bam files will be kept')
    parser.add_argument('-m', '--MaximumSRA', dest='MaximumSRA', type=int, default=10, help='The maximum amount of the sra files downloaded from NCBI. The default is 10')
    parser.add_argument('-T', '--threads_num', dest='threads_num', type=threads_type, default='auto',
                        help="Number of threads to use, or 'auto' to use all CPUs allocated to the job "
                             "(e.g. by Slurm --cpus-per-task). Values above the allocation are capped. "
                             "Default: %(default)s")
    parser.add_argument('-R', '--Run_prefix', dest='Run_prefix', default=False, action='store_true', help='if specified, use Run number as output prefix')
    args = parser.parse_args(argv)
    return args
