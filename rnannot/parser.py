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
    parser.add_argument('-a', '--assembly', dest='assembly', default='AssemblyName', help='The assembly name is used for naming output file.')
    parser.add_argument('-t', '--tempFile', dest='tempFile', default=False, action='store_true', help='if specified, intermediate output bam files will be kept')
    parser.add_argument('-m', '--MaximumSRA', dest='MaximumSRA', type=int, default=10, help='The maximum amount of the sra files downloaded from NCBI. The default is 10')
    parser.add_argument('-T', '--threads_num', dest='threads_num', type=str, default='24', help='The number of threads used in this program. The default is 24')
    parser.add_argument('-R', '--Run_prefix', dest='Run_prefix', default=False, action='store_true', help='if specified, use Run number as output prefix')
    args = parser.parse_args(argv)
    return args
