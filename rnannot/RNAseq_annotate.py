import os
from os import path
from parser import parse_args
from sys import argv
from utils import get_trimmomatic_jar_path, get_fastqc_path, get_adapter_path, get_hisat2_command_path
import subprocess
from zipfile import ZipFile
import gzip
import shutil

# parse the arguments, exclude the script name
args = parse_args(argv[1:])

# convert many arguments to absolute path
if not path.isabs(args.outdir):
    args.outdir = path.abspath(args.outdir)
if not path.isabs(args.genome):
    args.genome = path.abspath(args.genome)
if not path.isabs(args.file):
    args.file = path.abspath(args.file)

# create the output folder
output_prefix = path.join(args.outdir, args.name)
os.mkdir(output_prefix)

if args.genome.endswith('.gz'):
    new_genome_file_name = path.join(output_prefix, path.basename(args.genome).rstrip('.gz'))
    with gzip.open(args.genome, 'rb') as f_in:
        with open(new_genome_file_name, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    args.genome =  new_genome_file_name


sra_file_name = path.basename(args.file)
genome_file_name = path.basename(args.genome)
# convert SRA file to fastq file(s)
subprocess.run(['fastq-dump', '--split-files', '-O', output_prefix, args.file])


# run FastQC first, then use Trimmomatic to do trimming
fastqc_path = get_fastqc_path()
trimmomatic_jar_path = get_trimmomatic_jar_path()
if args.layout == 'single':
    subprocess.run([fastqc_path, '--outdir', output_prefix, path.join(output_prefix, sra_file_name + '_1.fastq')])
    with ZipFile(path.join(output_prefix, sra_file_name + '_1_fastqc.zip'), 'r') as zip_ref:
        zip_ref.extractall(output_prefix)
    if args.adaptor is None:
        if args.platform == 'ILLUMINA' and (args.model.startswith('Illumina HiSeq') or args.model.startswith('Illumina MiSeq')):
            subprocess.run(['java', '-jar', trimmomatic_jar_path, 'SE', path.join(output_prefix, sra_file_name + '_1.fastq'), path.join(output_prefix, 'output.fastq.gz'), 'ILLUMINACLIP:' + get_adapter_path('TruSeq3-SE.fa') + ':2:30:10', 'LEADING:3', 'TRAILING:3', 'SLIDINGWINDOW:4:15', 'MINLEN:36'])
        elif args.platform == 'ILLUMINA' and args.model.startswith('Illumina Genome Analyzer'):
            subprocess.run(['java', '-jar', trimmomatic_jar_path, 'SE', path.join(output_prefix, sra_file_name + '_1.fastq'), path.join(output_prefix, 'output.fastq'), 'ILLUMINACLIP:' + get_adapter_path('TruSeq2-SE.fa') + ':2:30:10', 'LEADING:3', 'TRAILING:3', 'SLIDINGWINDOW:4:15', 'MINLEN:36'])
        else:  # TODO: handle the cases that use other platforms and models, maybe we can guess adaptor from FastQC output
            pass
    else:  # TODO: use adaptor file if user has provided
        pass
    subprocess.run([get_hisat2_command_path('hisat2-build'), args.genome, path.join(output_prefix, genome_file_name)])
    subprocess.run([get_hisat2_command_path('hisat2'), '-x', path.join(output_prefix, genome_file_name), '-U', path.join(output_prefix, 'output.fastq'), '-S', path.join(output_prefix, 'output.sam')])
elif args.layout == 'paired':  # TODO: handle paired-end case, run FastQC first
    if args.adaptor is None:
        if args.platform == 'ILLUMINA' and (args.model.startswith('Illumina HiSeq') or args.model.startswith('Illumina MiSeq')):
            # TODO: Use Trimmomatic TruSeq3 as adaptor
            pass
        elif args.platform == 'ILLUMINA' and args.modlel.startswith('Illumina Genome Analyzer'):
            # TODO: Use Trimmomatic TruSeq2 as adaptor
            pass
        else:
            # TODO: Use BBTool (BBMerge) to determine the adaptor first, then run the Trimmomatic
            pass
    else:  # TODO: use adaptor file if user has provided
        pass
