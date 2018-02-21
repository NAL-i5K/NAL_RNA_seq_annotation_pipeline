import argparse
import datetime


def parse_args(argv):
    parser = argparse.ArgumentParser(description='Easy to use pipeline built for large-scale RNA-seq mapping with a genome assembly')
    parser.add_argument('files', metavar='files', type=str, nargs='+',
                        help='input files')
    parser.add_argument('--platform', nargs='?', default='ILLUMINA',
                        help='Platform of sequencing')
    parser.add_argument('--adaptor', nargs='?', default=None,
                        help='input adaptor file')
    parser.add_argument('--model', nargs='?', default=None,
                        help='Model of the seqencer')
    parser.add_argument('--layout', nargs='?', default=None)
    parser.add_argument('--name', nargs='?',
                        default=datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S"),
                        help='name of the output folder, if not specified, use the time of start')
    parser.add_argument('--outdir', nargs='?',
                        default='',
                        help='directory of output folder at, if not specified, use current folder')
    # option to determine it's a sra file or not
    # TODO
    parser.add_argument('--sra', help='If this argument is specified, the input should be a SRA run',
                        action='store_true')
    args = parser.parse_args(argv)
    if args.platform is None:
        raise Exception('No platform provided')

    if args.name is None:
        raise Exception('No name provided')

    # if args.sra is None:
    #     raise Exception('No sra run accession provide')

    # handle invalid input arguments
    if len(args.files) not in [1, 2]:
        raise Exception('Too many input files')
    elif args.sra and not (len(args.files) == 1):
        raise Exception('Expected only one file for sra')

    return args
