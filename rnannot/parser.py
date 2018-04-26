import argparse
import datetime
import sys


def parse_args(argv):
    parser = argparse.ArgumentParser(description='Easy to use pipeline built for large-scale RNA-seq mapping with a genome assembly')
    parser.add_argument('-f', '--file', dest='file', type=str,
                        help='A sra file that you want mapped to an genome')
    parser.add_argument('-g', '--genome', dest='genome')
    parser.add_argument('-p', '--platform', dest='platform', help='Platform of sequencing')
    parser.add_argument('-a', '--adapter', dest='adapter', nargs='?', help='input adaptor file')
    parser.add_argument('-m', '--model', dest='model', help='Model of the seqencer')
    parser.add_argument('-l', '--layout', help='single end or paried-end')
    parser.add_argument('-n', '--name', nargs='?',
                        default=datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S"),
                        help='name of the output folder, if not specified, use the time of start')
    parser.add_argument('-o', '--outdir', dest='outdir', nargs='?', default='.',
                        help='directory of output folder at, if not specified, use current folder')
    args = parser.parse_args(argv)
    if args.platform is None or args.name is None:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return args
