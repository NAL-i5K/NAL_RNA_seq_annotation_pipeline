import argparse
import datetime
import sys


def parse_args(argv):
    parser = argparse.ArgumentParser(description='Easy to use pipeline built for large-scale RNA-seq mapping with a genome assembly')
    parser.add_argument('-i', '--input', dest='input', type=str, help="A tsv file with a list of SRA runs' information.")
    parser.add_argument('-g', '--genome', dest='genome', help='A fasta file to align with.')
    parser.add_argument('-n', '--name', nargs='?',
                        default=datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S"),
                        help='name of the output folder, if not specified, use the time of start')
    parser.add_argument('-o', '--outdir', dest='outdir', nargs='?', default='.',
                        help='directory of output folder at, if not specified, use current folder')
    parser.add_argument('-d', '--downsample', dest='downsample', default=False,action='store_true', help='if specified, a downsampled bam file will be downsampled')
    if len(argv) == 0:
        parser.print_help(sys.stderr)
    args = parser.parse_args(argv)
    return args
