#!/usr/bin/env python3
import os
from os import path
from rnannot.parser import parse_args
from sys import argv, exit
from rnannot.utils import get_trimmomatic_jar_path, get_fastqc_path, get_trimmomatic_adapter_path, get_hisat2_command_path, get_bbmap_command_path, get_bbmap_adapter_path, get_gatk_jar_path, get_picard_jar_path
import subprocess
from zipfile import ZipFile
import gzip
import shutil
from itertools import islice
from six.moves import urllib
import datetime
import random

def run_pipeline(file, genome, outdir, name, layout, platform, model, download_link):
    # create the output folder
    output_prefix = path.join(outdir, name)
    os.mkdir(output_prefix)

    if platform == 'ABI_SOLID':
        return (
            False,
            'Currently, the colorspace data from ABI_SOLID is not supported', '')

    # check if SRA file exist or download it first
    if not path.exists(file):
        print('downloading sra files...')
        urllib.request.urlretrieve(download_link, file)

    # convert SRA file to fastq file(s)
    sra_file_name = path.basename(file)
    print('Unpacking the SRA file: {} ...'.format(file))
    f_stdout = open(
        path.join(output_prefix, sra_file_name + '.fastq-dump.log'), 'w')
    f_stderr = open(
        path.join(output_prefix, sra_file_name + '.fastq-dump.errlog'), 'w')
    subprocess.run(
        [
            'fastq-dump', '--dumpbase', '--split-files', '-O', output_prefix,
            file
        ],
        stdout=f_stdout,
        stderr=f_stderr)
    f_stdout.close()
    f_stderr.close()
    # Check if the SRA file is correct or not first
    if layout == 'PAIRED' and (not path.exists(
            path.join(
                output_prefix, sra_file_name + '_1.fastq')) or not path.exists(
                    path.join(output_prefix, sra_file_name + '_2.fastq'))):
        return (
            False,
            "run {} doesn't have paired data. It's not processed.".format(run), '')

    # Run FastQC first
    # Then, use Trimmomatic to do trimming
    # In the last step, perfom the alignment using HISAT2
    fastqc_path = get_fastqc_path()
    trimmomatic_jar_path = get_trimmomatic_jar_path()
    if layout == 'SINGLE':
        print('QC ...')
        f_stdout = open(
            path.join(output_prefix, sra_file_name + '_1.fastqc.log'), 'w')
        f_stderr = open(
            path.join(output_prefix, sra_file_name + '_1.fastqc.errlog'), 'w')
        subprocess.run(
            [
                fastqc_path, '--outdir', output_prefix,
                path.join(output_prefix, sra_file_name + '_1.fastq')
            ],
            stdout=f_stdout,
            stderr=f_stderr)
        with ZipFile(
                path.join(output_prefix, sra_file_name + '_1_fastqc.zip'),
                'r') as zip_ref:
            zip_ref.extractall(output_prefix)
        f_stdout.close()
        f_stderr.close()
        # Trimming with trimmomatic
        print('Trimming ...')
        f_stdout = open(
            path.join(output_prefix, sra_file_name + '.trimmomatic.log'), 'w')
        f_stderr = open(
            path.join(output_prefix, sra_file_name + '.trimmomatic.errlog'),
            'w')
        if platform == 'ILLUMINA' and (model.startswith('Illumina HiSeq')
                                       or model.startswith('Illumina MiSeq')):
            subprocess.run(
                [
                    'java', '-jar', trimmomatic_jar_path, 'SE',
                    path.join(output_prefix, sra_file_name + '_1.fastq'),
                    path.join(output_prefix, 'output.fastq'), 'ILLUMINACLIP:' +
                    get_trimmomatic_adapter_path('TruSeq3-SE.fa') + ':2:30:10',
                    'LEADING:30', 'TRAILING:30', 'SLIDINGWINDOW:4:15',
                    'MINLEN:36', 'TOPHRED33'
                ],
                stdout=f_stdout,
                stderr=f_stderr)
        elif platform == 'ILLUMINA' and model.startswith(
                'Illumina Genome Analyzer II'):
            subprocess.run(
                [
                    'java', '-jar', trimmomatic_jar_path, 'SE',
                    path.join(output_prefix, sra_file_name + '_1.fastq'),
                    path.join(output_prefix, 'output.fastq'), 'ILLUMINACLIP:' +
                    get_trimmomatic_adapter_path('TruSeq2-SE.fa') + ':2:30:10',
                    'LEADING:30', 'TRAILING:30', 'SLIDINGWINDOW:4:15',
                    'MINLEN:36', 'TOPHRED33'
                ],
                stdout=f_stdout,
                stderr=f_stderr)
        else:
            # Use adapter file from BBMap for other platforms and models.
            subprocess.run(
                [
                    'java', '-jar', trimmomatic_jar_path, 'SE',
                    path.join(output_prefix, sra_file_name + '_1.fastq'),
                    path.join(output_prefix, 'output.fastq'),
                    'ILLUMINACLIP:' + get_bbmap_adapter_path() + ':2:30:10',
                    'LEADING:30', 'TRAILING:30', 'SLIDINGWINDOW:4:15',
                    'MINLEN:36', 'TOPHRED33'
                ],
                stdout=f_stdout,
                stderr=f_stderr)
        f_stdout.close()
        f_stderr.close()
        # Normalizing
        print('Normalizing ...')   
        subprocess.run([get_bbmap_command_path('bbnorm.sh'), 'in1=' + path.join(output_prefix, 'output.fastq'), 
                        'out=' + path.join(output_prefix, 'normalized.fastq'), 'target=' + '100', 'threads=' + '24'])
        return (True, '', 'single')

    elif layout == 'PAIRED':
        print('QC ...')
        f_stdout = open(
            path.join(output_prefix, sra_file_name + '_1.fastqc.log'), 'w')
        f_stderr = open(
            path.join(output_prefix, sra_file_name + '_1.fastqc.errlog'), 'w')
        subprocess.run(
            [
                fastqc_path, '--outdir', output_prefix,
                path.join(output_prefix, sra_file_name + '_1.fastq')
            ],
            stdout=f_stdout,
            stderr=f_stderr)
        f_stdout.close()
        f_stderr.close()
        f_stdout = open(
            path.join(output_prefix, sra_file_name + '_2.fastqc.log'), 'w')
        f_stderr = open(
            path.join(output_prefix, sra_file_name + '_2.fastqc.errlog'), 'w')
        subprocess.run(
            [
                fastqc_path, '--outdir', output_prefix,
                path.join(output_prefix, sra_file_name + '_2.fastq')
            ],
            stdout=f_stdout,
            stderr=f_stderr)
        with ZipFile(
                path.join(output_prefix, sra_file_name + '_1_fastqc.zip'),
                'r') as zip_ref:
            zip_ref.extractall(output_prefix)
        with ZipFile(
                path.join(output_prefix, sra_file_name + '_2_fastqc.zip'),
                'r') as zip_ref:
            zip_ref.extractall(output_prefix)
        # Trimming with trimmomatic
        print('Trimming ...')
        f_stdout = open(
            path.join(output_prefix, sra_file_name + '.trimmomatic.log'), 'w')
        f_stderr = open(
            path.join(output_prefix, sra_file_name + '.trimmomatic.errlog'),
            'w')
        if platform == 'ILLUMINA' and (model.startswith('Illumina HiSeq')
                                       or model.startswith('Illumina MiSeq')):
            subprocess.run(
                [
                    'java', '-jar', trimmomatic_jar_path, 'PE',
                    path.join(output_prefix, sra_file_name + '_1.fastq'),
                    path.join(output_prefix, sra_file_name + '_2.fastq'),
                    path.join(output_prefix, 'output_1.fastq'),
                    path.join(output_prefix, 'output_1_un.fastq'),
                    path.join(output_prefix, 'output_2.fastq'),
                    path.join(output_prefix,
                              'output_2_un.fastq'), 'ILLUMINACLIP:' +
                    get_trimmomatic_adapter_path('TruSeq3-PE.fa') + ':2:30:10',
                    'LEADING:30', 'TRAILING:30', 'SLIDINGWINDOW:4:15',
                    'MINLEN:36', 'TOPHRED33'
                ],
                stdout=f_stdout,
                stderr=f_stderr)
        elif platform == 'ILLUMINA' and model.startswith(
                'Illumina Genome Analyzer II'):
            subprocess.run(
                [
                    'java', '-jar', trimmomatic_jar_path, 'PE',
                    path.join(output_prefix, sra_file_name + '_1.fastq'),
                    path.join(output_prefix, sra_file_name + '_2.fastq'),
                    path.join(output_prefix, 'output_1.fastq'),
                    path.join(output_prefix, 'output_1_un.fastq'),
                    path.join(output_prefix, 'output_2.fastq'),
                    path.join(output_prefix,
                              'output_2_un.fastq'), 'ILLUMINACLIP:' +
                    get_trimmomatic_adapter_path('TruSeq2-PE.fa') + ':2:30:10',
                    'LEADING:30', 'TRAILING:30', 'SLIDINGWINDOW:4:15',
                    'MINLEN:36', 'TOPHRED33'
                ],
                stdout=f_stdout,
                stderr=f_stderr)
        else:
            # use BBTool (BBMerge) to determine the adapter first, then run the Trimmomatic
            f_bbmap_stdout = open(
                path.join(output_prefix, sra_file_name + '.bbmap.log'), 'w')
            f_bbmap_stderr = open(
                path.join(output_prefix, sra_file_name + '.bbmap.errlog'), 'w')
            subprocess.run(
                [
                    get_bbmap_command_path('bbmerge.sh'), 'in1=' + path.join(
                        output_prefix, sra_file_name + '_1.fastq'), 'in2=' +
                    path.join(output_prefix, sra_file_name + '_2.fastq'),
                    'outa=' + path.join(output_prefix, 'adapters.fa')
                ],
                stdout=f_bbmap_stdout,
                stderr=f_bbmap_stderr)
            f_bbmap_stdout.close()
            f_bbmap_stderr.close()
            subprocess.run(
                [
                    'java', '-jar', trimmomatic_jar_path, 'PE',
                    path.join(output_prefix, sra_file_name + '_1.fastq'),
                    path.join(output_prefix, sra_file_name + '_2.fastq'),
                    path.join(output_prefix, 'output_1.fastq'),
                    path.join(output_prefix, 'output_1_un.fastq'),
                    path.join(output_prefix, 'output_2.fastq'),
                    path.join(output_prefix, 'output_2_un.fastq'),
                    'ILLUMINACLIP:' + path.join(output_prefix, 'adapters.fa') +
                    ':2:30:10', 'LEADING:30', 'TRAILING:30',
                    'SLIDINGWINDOW:4:15', 'MINLEN:36', 'TOPHRED33'
                ],
                stdout=f_stdout,
                stderr=f_stderr)
        f_stdout.close()
        f_stderr.close()
        # Normalizing
        print('Normalizing ...')
        subprocess.run([get_bbmap_command_path('bbnorm.sh'), 'in1=' + path.join(output_prefix,'output_1.fastq'), 'in2=' + path.join(output_prefix, 'output_2.fastq'),
                       'out1=' + path.join(output_prefix, 'normalized_1.fastq'), 'out2=' + path.join(output_prefix, 'normalized_2.fastq'), 'target=' + '100', 'threads=' + '24'])
        return (True, '', 'paired')


def merge_files(files, outdir):  # merge sam files
    print('Combing the sam/bam files ...')
    f_stdout = open(path.join(outdir, 'out.log'), 'a')
    f_stderr = open(path.join(outdir, 'out.errlog'), 'a')
    args = [
        'java', '-jar',
        get_picard_jar_path(), 'MergeSamFiles',
        'O=' + path.join(outdir, 'output.bam')
    ]
    args += ['I=' + f for f in files]
    subprocess.run(args, stdout=f_stdout, stderr=f_stderr)
    print('Finished combining the sam/bam files')


def check_ref_files(ref_path):
    if path.exists(ref_path + '.fai') and path.exists('.dict'):
        return True
    else:
        return False


def read_sam_errors(file_path):
    warns = set()
    errors = set()
    with open(file_path) as f:
        for line in islice(f, 4, None):
            temp = line.split('\t')[0]
            if 'ERROR' in temp:
                errors.add(temp.lstrip('ERROR:'))
            elif 'WARNING' in temp:
                warns.add(temp.lstrip('WARNING:'))
            elif temp == '\n':
                break
    return (errors, warns)

if __name__ == '__main__':
    
    # parse the arguments, exclude the script name
    args = parse_args(argv[1:])
    # convert many arguments to absolute path
    if not path.isabs(args.outdir):
        args.outdir = path.abspath(args.outdir)
    if not path.isabs(args.input):
        args.input = path.abspath(args.input)
    if not path.isabs(args.genome):
        args.genome = path.abspath(args.genome)
   
    os.mkdir(path.join(args.outdir, args.name))
    if args.genome.endswith('.gz'):
        new_genome_file_name = path.join(
            args.outdir, args.name,
            path.basename(args.genome).rstrip('.gz'))
        with gzip.open(args.genome, 'rb') as f_in:
            with open(new_genome_file_name, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        args.genome = new_genome_file_name
    
    with open(args.input) as f:
        col_names = f.readline().rstrip('\n').split('\t')
        run_ind = col_names.index('Run')
        platform_ind = col_names.index('Platform')
        model_ind = col_names.index('Model')
        layout_ind = col_names.index('LibraryLayout')
        download_ind = col_names.index('download_path')
        scientific_name_ind = col_names.index('ScientificName')
        print('Checking the input tsv file: {}'.format(args.input))
        for ind, name in zip([run_ind, platform_ind, model_ind, layout_ind, download_ind, scientific_name_ind],
                             ['Run', 'Platform', 'Model', 'LibraryLayout', 'download_path', 'ScientificName']):
            if ind == -1:
                print('{} column is missing in input tsv file.'.format(name))
                exit(1)
        runs = []
        platforms = []
        models = []
        layouts = []
        download_links = []
        scientific_names = []
        for line in f:
            temp = line.rstrip('\n').split('\t')
            runs.append(temp[run_ind])
            platforms.append(temp[platform_ind])
            models.append(temp[model_ind])
            layouts.append(temp[layout_ind])
            download_links.append(temp[download_ind])
            scientific_names.append(temp[scientific_name_ind])
    
        # check the amount of sra files in tsv
        if len(runs) > args.MaximumSRA:
            print('The amount of sra files is more than {}'.format(args.MaximumSRA))
            print('Randomly pick {} sra files for downloading'.format(args.MaximumSRA)) 
            runs_temp = []
            platforms_temp = []
            models_temp = []
            layouts_temp = []
            download_links_temp = []
            scientific_names_temp = []
            random_sra = random.sample(range(0,len(runs)-1), args.MaximumSRA)
            # randomly pick the maximum amount of sra files to download
            for i in random_sra:
                runs_temp.append(runs[i])
                platforms_temp.append(platforms[i])
                models_temp.append(models[i])
                layouts_temp.append(layouts[i])
                download_links_temp.append(download_links[i])
                scientific_names_temp.append(scientific_names[i])
            runs = runs_temp
            platforms = platforms_temp
            models = models_temp
            layouts = layouts_temp
            download_links = download_links_temp
            scientific_names = scientific_names_temp
    # output runs/scientific_name/assembly_name to Source.txt
    with open(path.join(args.outdir, args.name, 'Source.txt'), 'w') as outfile:
        outfile.write(scientific_names[0] + '\n')
        outfile.write(args.assembly + '\n')
        for run in runs:
            outfile.write(run + '\n')
 
    single_files_for_merge = []
    paired1_files_for_merge = []
    paired2_files_for_merge = []
    for run, platform, model, layout, download_link in zip(runs, platforms, models, layouts, download_links):
        print('Processing the file: {}'.format(run))
        if not path.isabs(run):
            run = path.abspath(run)
        run_file_name = path.basename(run)
        return_status, err_message, return_layout = run_pipeline(
            file=run,
            genome=args.genome,
            outdir=path.join(args.outdir, args.name),
            name=run_file_name,
            layout=layout,
            platform=platform,
            model=model,
            download_link=download_link
            )
        if return_status:
            if return_layout == 'single':
                single_files_for_merge.append(
                    path.join(args.outdir, args.name, run_file_name, 'normalized.fastq'))
            if return_layout == 'paired':
                paired1_files_for_merge.append(
                    path.join(args.outdir, args.name, run_file_name, 'normalized_1.fastq'))
                paired2_files_for_merge.append(
                    path.join(args.outdir, args.name, run_file_name, 'normalized_2.fastq'))
        else:
            print(err_message)
    
    # merge normalized fastq files together
    with open(path.join(args.outdir, args.name, 'merged_normalized.fastq'), 'w') as outfile:
        for fname in single_files_for_merge:
            with open(fname, 'r') as infile:
                 for line in infile:
                     outfile.write(line)
    with open(path.join(args.outdir, args.name, 'merged_normalized_1.fastq'), 'w') as outfile:
        for fname in paired1_files_for_merge:
            with open(fname, 'r') as infile:
                 for line in infile:
                     outfile.write(line)
    with open(path.join(args.outdir, args.name, 'merged_normalized_2.fastq'), 'w') as outfile:
        for fname in paired2_files_for_merge:
            with open(fname, 'r') as infile:
                 for line in infile:
                     outfile.write(line)
    # Decompress the gz file, becasue some of tools don't accept .gz compressed files
    genome = args.genome
    if genome.endswith('.gz'):
        new_genome_file_name = path.join(args.outdir, args.name,
                                         path.basename(genome).rstrip('.gz'))
        with gzip.open(genome, 'rb') as f_in:
            with open(new_genome_file_name, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        genome = new_genome_file_name
    genome_file_name = path.basename(genome)    
    # Aligning SINGLE fastq with HISAT2 
    print('Aligning single file...')
    f_stdout = open(
        path.join(args.outdir, args.name, genome_file_name + '_single.hisat2.log'), 'w')
    f_stderr = open(
        path.join(args.outdir, args.name, genome_file_name + '_single.hisat2.errlog'), 'w')
    subprocess.run(
        [
            get_hisat2_command_path('hisat2-build'), genome,
            path.join(args.outdir, args.name, genome_file_name)
        ],
        stdout=f_stdout,
        stderr=f_stderr)
    subprocess.run(
        [
            get_hisat2_command_path('hisat2'), '--no-mixed', '--no-discordant', '-p', '12', '-x',
            path.join(args.outdir, args.name, genome_file_name), '-U',
            path.join(args.outdir, args.name, 'merged_normalized.fastq'), '-S',
            path.join(args.outdir, args.name, 'single_output.sam')
        ],
        stdout=f_stdout,
        stderr=f_stderr)
    f_stdout.close()
    f_stderr.close()
    # Aligning PAIRED fastq with HISAT2
    print('Aligning paired file...')
    f_stdout = open(
        path.join(args.outdir, args.name, genome_file_name + '_paired.hisat2-build.log'), 'w')
    f_stderr = open(
        path.join(args.outdir, args.name, genome_file_name + '_paired.hisat2-build.errlog'), 'w')
    subprocess.run(
        [
            get_hisat2_command_path('hisat2-build'), genome,
            path.join(args.outdir, args.name, genome_file_name)
        ],
        stdout=f_stdout,
        stderr=f_stderr)
    f_stdout.close()
    f_stderr.close()
    f_stdout = open(
        path.join(args.outdir, args.name, genome_file_name + '_paired.hisat2.log'), 'w')
    f_stderr = open(
        path.join(args.outdir, args.name, genome_file_name + '_paired.hisat2.errlog'), 'w')
    subprocess.run(
        [
            get_hisat2_command_path('hisat2'), '--no-mixed', '--no-discordant', '-p', '12', '-x',
            path.join(args.outdir, args.name, genome_file_name), '-1',
            path.join(args.outdir, args.name, 'merged_normalized_1.fastq'), '-2',
            path.join(args.outdir, args.name, 'merged_normalized_2.fastq'), '-S',
            path.join(args.outdir, args.name, 'paired_output.sam')
        ],
        stdout=f_stdout,
        stderr=f_stderr)
    f_stdout.close()
    f_stderr.close()
    # sort and convert to the bam file-single
    f_stdout = open(
        path.join(args.outdir, args.name, genome_file_name + '_single.samtools.log'), 'w')
    f_stderr = open(
        path.join(args.outdir, args.name, genome_file_name + '_single.samtools.errlog'), 'w')
    subprocess.run(
        [
            'samtools', 'sort', '-@', '12', '-o',
            path.join(args.outdir, args.name, 'single_output.bam'), '-O', 'bam', '-T',
            path.join(args.outdir, args.name, 'single_output'),
            path.join(args.outdir, args.name, 'single_output.sam') 
        ],
        stdout=f_stdout,
        stderr=f_stderr)
    f_stdout.close()
    f_stderr.close()
    # sort and convert to the bam file-paired
    f_stdout = open(
        path.join(args.outdir, args.name, genome_file_name + '_paired.samtools.log'), 'w')
    f_stderr = open(
        path.join(args.outdir, args.name, genome_file_name + '_paired.samtools.errlog'), 'w')
    subprocess.run(
        [
            'samtools', 'sort', '-@', '12', '-o',
            path.join(args.outdir, args.name, 'paired_output.bam'), '-O', 'bam', '-T',
            path.join(args.outdir, args.name, 'paired_output'),
            path.join(args.outdir, args.name, 'paired_output.sam')
        ],
        stdout=f_stdout,
        stderr=f_stderr)
    f_stdout.close()
    f_stderr.close()
    # combine the sam/bam files together and convert to BAM file
    files_for_merge = [path.join(args.outdir, args.name, 'single_output.bam'), path.join(args.outdir, args.name, 'paired_output.bam')]
    merge_files(files_for_merge, path.join(args.outdir, args.name))
    # sorting bam file
    print('Start sorting bam file...')
    bam_dir = path.join(args.outdir, args.name, 'output.bam')
    output_dir = path.join(args.outdir, args.name, 'output.sorted.bam')
    subprocess.run(['samtools', 'sort', '-@', '12', bam_dir, '-o', output_dir])
    # converting dowsampled bam to bigwig (includes indexing bam file)
    print('Generating bigwig file from bam file...')
    bam_dir = path.join(args.outdir, args.name, 'output.sorted.bam')
    bigwig_dir = path.join(args.outdir, args.name, 'output.bigwig')
    subprocess.run(['python3', 'bam_to_bigwig.py', bam_dir, '-o', bigwig_dir])
    # rename bam and bigwig file to [gggsss]_[assembly_name]_RNA-Seq-alignments_[datetime]
    temp = scientific_names[0].split(" ")
    gene_name = temp[0]
    species_name = temp[1]
    new_name = gene_name[0:3] + species_name[0:3]  + '_' + args.assembly + '_RNA-Seq-alignments_' + datetime.datetime.now().strftime("%Y-%m-%d")
    os.rename(path.join(args.outdir, args.name, 'output.sorted.bam'), path.join(args.outdir, args.name, new_name + '.bam'))
    os.rename(path.join(args.outdir, args.name, 'output.sorted.bam.bai'), path.join(args.outdir, args.name, new_name + '.bam.bai'))
    os.rename(path.join(args.outdir, args.name, 'output.bigwig'), path.join(args.outdir, args.name, new_name + '.bigwig'))
    #  remove intermediate files
    if not args.tempFile:
        os.remove(path.join(args.outdir, args.name, 'output.bam'))
        os.remove(path.join(args.outdir, args.name, 'single_output.bam'))
        os.remove(path.join(args.outdir, args.name, 'paired_output.bam'))
        os.remove(path.join(args.outdir, args.name, 'single_output.sam'))
        os.remove(path.join(args.outdir, args.name, 'paired_output.sam'))
    print('Finished processing')
    
