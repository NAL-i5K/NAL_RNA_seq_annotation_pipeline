import os
from os import path
from parser import parse_args
from sys import argv

args = parse_args(argv[1:])  # exclude the script name


# if it's a sra run, split first
# TODO

cwd = os.getcwd()
# create the output folder
os.mkdir(path.join(cwd, args.name))

# Run Multi QC first
# TODO

# determine trimming
if args.adaptor is None:
    if args.platform == 'ILLUMINA' and (args.model.startswith('Illumina HiSeq') or args.model.startswith('Illumina MiSeq')):
        # Use Trimmomatic TruSeq3 as adaptor
        # TODO
        pass
    elif args.platform == 'ILLUMINA' and args.modlel.startswith('Illumina Genome Analyzer'):
        # Use Trimmomatic TruSeq2 as adaptor
        # TODO
        pass
    else:
        if args.layout == 'paired':
            # Use BBTool (BBMerge) to determine the adaptor first, then run the Trimmomatic
            # TODO
            pass
        else:  # single end
            # TODO
            pass
else:  # use adaptor file if user has provided
    # TODO
    pass
