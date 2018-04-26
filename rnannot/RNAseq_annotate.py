import os
from os import path
from parser import parse_args
from sys import argv
from utils import get_trimmomatic_jar_path, get_fastqc_path, get_adapter_path, get_hisat2_command_path
import subprocess
from zipfile import ZipFile
import gzip
import shutil



def run_pipeline(file, genome, outdir, name, layout, platform, model, adapter):
    # create the output folder
    output_prefix = path.join(outdir, name)
    os.mkdir(output_prefix)
    
    # decompress the gz file, becasue some of tools don't accerpt .gz compressed files
    if genome.endswith('.gz'):
        new_genome_file_name = path.join(output_prefix, path.basename(genome).rstrip('.gz'))
        with gzip.open(genome, 'rb') as f_in:
            with open(new_genome_file_name, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        genome =  new_genome_file_name
    sra_file_name = path.basename(file)
    genome_file_name = path.basename(genome)

    # convert SRA file to fastq file(s)
    subprocess.run(['fastq-dump', '--split-files', '-O', output_prefix, file])

    # run FastQC first, then use Trimmomatic to do trimming
    fastqc_path = get_fastqc_path()
    trimmomatic_jar_path = get_trimmomatic_jar_path()
    if layout == 'single':
        subprocess.run([fastqc_path, '--outdir', output_prefix, path.join(output_prefix, sra_file_name + '_1.fastq')])
        with ZipFile(path.join(output_prefix, sra_file_name + '_1_fastqc.zip'), 'r') as zip_ref:
            zip_ref.extractall(output_prefix)
        if adapter is None:
            if platform == 'ILLUMINA' and (model.startswith('Illumina HiSeq') or model.startswith('Illumina MiSeq')):
                subprocess.run(['java', '-jar', trimmomatic_jar_path, 'SE', path.join(output_prefix, sra_file_name + '_1.fastq'), path.join(output_prefix, 'output.fastq.gz'), 'ILLUMINACLIP:' + get_adapter_path('TruSeq3-SE.fa') + ':2:30:10', 'LEADING:3', 'TRAILING:3', 'SLIDINGWINDOW:4:15', 'MINLEN:36'])
            elif platform == 'ILLUMINA' and model.startswith('Illumina Genome Analyzer'):
                subprocess.run(['java', '-jar', trimmomatic_jar_path, 'SE', path.join(output_prefix, sra_file_name + '_1.fastq'), path.join(output_prefix, 'output.fastq'), 'ILLUMINACLIP:' + get_adapter_path('TruSeq2-SE.fa') + ':2:30:10', 'LEADING:3', 'TRAILING:3', 'SLIDINGWINDOW:4:15', 'MINLEN:36'])
            else:  # TODO: handle the cases that use other platforms and models, maybe we can guess adapter from FastQC output
                pass
        else:  # TODO: use adapter file if user has provided
            pass
        subprocess.run([get_hisat2_command_path('hisat2-build'), genome, path.join(output_prefix, genome_file_name)])
        subprocess.run([get_hisat2_command_path('hisat2'), '-x', path.join(output_prefix, genome_file_name), '-U', path.join(output_prefix, 'output.fastq'), '-S', path.join(output_prefix, 'output.sam')])
    elif layout == 'paired':  # TODO: handle paired-end case, run FastQC first
        if adapter is None:
            if platform == 'ILLUMINA' and (model.startswith('Illumina HiSeq') or model.startswith('Illumina MiSeq')):
                # TODO: Use Trimmomatic TruSeq3 as adapter
                pass
            elif platform == 'ILLUMINA' and modlel.startswith('Illumina Genome Analyzer II'):
                # TODO: Use Trimmomatic TruSeq2 as adapter
                pass
            else:
                # TODO: Use BBTool (BBMerge) to determine the adapter first, then run the Trimmomatic
                pass
        else:  # TODO: use adapter file if user has provided
            pass


if __name__ == '__main__':
    # parse the arguments, exclude the script name
    args = parse_args(argv[1:])
    # convert many arguments to absolute path
    if not path.isabs(args.outdir):
        args.outdir = path.abspath(args.outdir)
    if not path.isabs(args.genome):
        args.genome = path.abspath(args.genome)
    if not path.isabs(args.file):
        args.file = path.abspath(args.file)
    if args.adapter is not None and not path.isabs(args.adapter):
        args.adapter = path.abspath(args.adapter)

    run_pipeline(file=args.file, genome=args.genome, outdir=args.outdir, name=args.name, layout=args.layout, platform=args.platform, model=args.model, adapter=args.adapter)
