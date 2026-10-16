#!/usr/bin/env python3
import os
from os import path
from rnannot.parser import parse_args
from sys import argv, exit
from rnannot.utils import get_trimmomatic_jar_path, get_fastqc_path, get_trimmomatic_adapter_path, get_hisat2_command_path, get_bbmap_command_path, get_bbmap_adapter_path, get_picard_jar_path, get_bam_to_bigwig_path, get_regtools_path
import subprocess
from zipfile import ZipFile
import gzip
import shutil
from itertools import islice
import urllib.request
import urllib.parse
import urllib.error
import datetime
import random
import glob


def resolve_threads(requested):
    """Return the thread count as a string, based on the CPUs this job may use.

    os.sched_getaffinity reflects the Slurm allocation (inherited inside
    Apptainer); os.cpu_count would report every CPU on the node.
    """
    try:
        available = len(os.sched_getaffinity(0))
    except AttributeError:  # non-Linux fallback
        available = os.cpu_count() or 1
    if requested is None or str(requested).strip().lower() in ('', '0', 'auto'):
        return str(available)
    n = int(requested)
    if n > available:
        print('Requested {} threads but only {} CPUs are available; using {}'.format(
            n, available, available))
        n = available
    return str(max(n, 1))


def run_pipeline(run, genome, outdir, layout, platform, model, threads):
    # create the output folder
    output_prefix = path.join(outdir, run)
    os.mkdir(output_prefix)
    if platform == 'ABI_SOLID':
        return (
            False,
            'Currently, the colorspace data from ABI_SOLID is not supported', '')
    # download SRA files and convert SRA file to fastq file(s)
    print('downloading sra files...')
    sra_file_name = run
    f_stdout = open(
        path.join(output_prefix, sra_file_name + '.fastq-dump.log'), 'w')
    f_stderr = open(
        path.join(output_prefix, sra_file_name + '.fastq-dump.errlog'), 'w')
    if layout == 'SINGLE':
        pwd = os.getcwd()
        subprocess.run(
            [
                'fasterq-dump', '-O', output_prefix, '-t', pwd,
                run
            ],
            stdout=f_stdout,
            stderr=f_stderr)
        f_stdout.close()
        f_stderr.close()
    elif layout == 'PAIRED':
        pwd = os.getcwd()
        subprocess.run(
            [
                'fasterq-dump', '-O', output_prefix, '-t', pwd,
                run
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
    # In the last step, perfom normalizing using bbnorm
    fastqc_path = get_fastqc_path()
    trimmomatic_jar_path = get_trimmomatic_jar_path()
    if layout == 'SINGLE':
        print('QC ...')
        f_stdout = open(
            path.join(output_prefix, sra_file_name + '.fastqc.log'), 'w')
        f_stderr = open(
            path.join(output_prefix, sra_file_name + '.fastqc.errlog'), 'w')
        subprocess.run(
            [
                fastqc_path, '--outdir', output_prefix,
                path.join(output_prefix, sra_file_name + '.fastq')
            ],
            stdout=f_stdout,
            stderr=f_stderr)
        with ZipFile(
                path.join(output_prefix, sra_file_name + '_fastqc.zip'),
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
                                       or model.startswith('Illumina MiSeq')) \
               or model.startswith('Element AVITI'):
            subprocess.run(
                [
                    'java', '-jar', trimmomatic_jar_path, 'SE', '-phred33',
                    path.join(output_prefix, sra_file_name + '.fastq'),
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
                    'java', '-jar', trimmomatic_jar_path, 'SE', '-phred33',
                    path.join(output_prefix, sra_file_name + '.fastq'),
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
                    'java', '-jar', trimmomatic_jar_path, 'SE', '-phred33',
                    path.join(output_prefix, sra_file_name + '.fastq'),
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
        subprocess.run([get_bbmap_command_path('bbnorm.sh'), 'qin=33', 'in=' + path.join(output_prefix, 'output.fastq'),
                        'out=' + path.join(output_prefix, 'normalized.fastq'), 'target=' + '100', 'prefilter=t',
                        'threads=' + threads], check=True)
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
                                       or model.startswith('Illumina MiSeq')) \
               or model.startswith('Element AVITI'):
            subprocess.run(
                [
                    'java', '-jar', trimmomatic_jar_path, 'PE', '-phred33',
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
                    'java', '-jar', trimmomatic_jar_path, 'PE',  '-phred33',
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
                    'java', '-jar', trimmomatic_jar_path, 'PE', '-phred33',
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
        subprocess.run([get_bbmap_command_path('repair.sh'), 'in1=' + path.join(output_prefix,'output_1.fastq'), 'in2=' + path.join(output_prefix, 'output_2.fastq'),
                       'out1=' + path.join(output_prefix, 'fixed_output_1.fastq'), 'out2=' + path.join(output_prefix, 'fixed_output_2.fastq'), 'outsingle=' + path.join(output_prefix, 'orphans.fastq')], check=True)
        subprocess.run([get_bbmap_command_path('bbnorm.sh'), 'qin=33', 'in1=' + path.join(output_prefix,'fixed_output_1.fastq'), 'in2=' + path.join(output_prefix, 'fixed_output_2.fastq'),
                       'out1=' + path.join(output_prefix, 'normalized_1.fastq'), 'out2=' + path.join(output_prefix, 'normalized_2.fastq'), 'target=' + '100', 'prefilter=t',
                       'threads=' + threads], check=True)
        return (True, '', 'paired')


def merge_files(files, outdir, java_mem='4g'):  # merge sam/bam files
    print('Combining the sam/bam files ...')
    args = [
        'java', '-Xmx' + java_mem, '-jar',
        get_picard_jar_path(), 'MergeSamFiles',
        'O=' + path.join(outdir, 'output.bam'),
    ]
    args += ['I=' + f for f in files]
    with open(path.join(outdir, 'out.log'), 'a') as f_stdout, \
         open(path.join(outdir, 'out.errlog'), 'a') as f_stderr:
        result = subprocess.run(args, stdout=f_stdout, stderr=f_stderr)
    if result.returncode != 0:
        raise RuntimeError(
            'Picard MergeSamFiles failed (exit code {}); see {}'.format(
                result.returncode, path.join(outdir, 'out.errlog')))
    print('Finished combining the sam/bam files')


def build_hisat2_index(genome, index_prefix, threads):
    """Build the HISAT2 index once per run; skip if it already exists."""
    # .ht2 for genomes < ~4 Gb, .ht2l for larger genomes
    if path.exists(index_prefix + '.1.ht2') or path.exists(index_prefix + '.1.ht2l'):
        print('HISAT2 index found, skipping build')
        return
    print('Building HISAT2 index...')
    # Log names use '_' (not '.') so the cleanup glob genome_file_name + '.*' won't delete them
    with open(index_prefix + '_hisat2-build.log', 'w') as f_out, \
         open(index_prefix + '_hisat2-build.errlog', 'w') as f_err:
        subprocess.run([get_hisat2_command_path('hisat2-build'), '-p', str(threads),
                        genome, index_prefix],
                       stdout=f_out, stderr=f_err, check=True)


def is_coordinate_sorted(bam_path):
    """True if the BAM header's @HD line declares SO:coordinate."""
    header = subprocess.run(['samtools', 'view', '-H', bam_path],
                            capture_output=True, text=True, check=True).stdout
    for line in header.splitlines():
        if line.startswith('@HD'):
            return 'SO:coordinate' in line.split('\t')
    return False


def ensure_coordinate_sorted(bam_path, threads):
    """Sort bam_path in place only if its header says it isn't coordinate-sorted."""
    if is_coordinate_sorted(bam_path):
        print('BAM is already coordinate-sorted, skipping sort')
        return
    print('Sorting bam file...')
    tmp_path = bam_path + '.sorting.bam'
    subprocess.run(['samtools', 'sort', '-@', str(threads), '-o', tmp_path, bam_path],
                   check=True)
    os.replace(tmp_path, bam_path)


def make_bigwig(bam_path, bigwig_path, threads):
    """Run bam_to_bigwig.py (unmodified) and verify that it actually worked."""
    # 1. CSI index, so bam_to_bigwig.py skips its own pysam.index() (.bai) step
    subprocess.run(['samtools', 'index', '-c', '-@', str(threads), bam_path],
                   check=True)
    stale_bai = bam_path + '.bai'
    if path.exists(stale_bai):
        os.remove(stale_bai)   # a .bai would make the script re-index
    # 2. Remove old output so the script never reaches its interactive prompt
    if path.exists(bigwig_path):
        os.remove(bigwig_path)
    # 3. Run the script; stdin=DEVNULL guarantees it can never hang on input()
    subprocess.run(['python3', get_bam_to_bigwig_path(), bam_path,
                    '-o', bigwig_path],
                   stdin=subprocess.DEVNULL)
    # 4. The script exits 0 even on failure, so check the result directly
    if not path.exists(bigwig_path) or path.getsize(bigwig_path) == 0:
        prefix = path.splitext(bam_path)[0]
        for leftover in (prefix + '-filtered.bam', prefix + '.wig',
                         prefix + '.sizes'):
            if path.exists(leftover):
                os.remove(leftover)
        raise RuntimeError('bam_to_bigwig.py did not produce ' + bigwig_path)


def check_ref_files(ref_path):
    dict_path = path.splitext(ref_path)[0] + '.dict'
    return path.exists(ref_path + '.fai') and path.exists(dict_path)


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
    # use the CPUs actually allocated to this job (Slurm/Apptainer aware)
    args.threads_num = resolve_threads(getattr(args, 'threads_num', None))
    print('Using {} threads'.format(args.threads_num))
    # convert many arguments to absolute path
    if not path.isabs(args.outdir):
        args.outdir = path.abspath(args.outdir)
    if not path.isabs(args.input):
        args.input = path.abspath(args.input)
    if not path.isabs(args.genome):
        args.genome = path.abspath(args.genome)

    os.mkdir(path.join(args.outdir, args.name))
    # Decompress the gz file, becasue some of tools don't accept .gz compressed files
    genome = args.genome
    if genome.endswith('.gz'):
        new_genome_file_name = path.join(args.outdir, args.name,
                                         path.basename(genome)[:-3])  # drop exactly '.gz'
        with gzip.open(genome, 'rb') as f_in:
            with open(new_genome_file_name, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        genome = new_genome_file_name
    genome_file_name = path.basename(genome)
    with open(args.input) as f:
        col_names = f.readline().rstrip('\n').split('\t')
        run_ind = col_names.index('Run')
        platform_ind = col_names.index('Platform')
        model_ind = col_names.index('Model')
        layout_ind = col_names.index('LibraryLayout')
        scientific_name_ind = col_names.index('ScientificName')
        print('Checking the input tsv file: {}'.format(args.input))
        for ind, name in zip([run_ind, platform_ind, model_ind, layout_ind, scientific_name_ind],
                             ['Run', 'Platform', 'Model', 'LibraryLayout', 'ScientificName']):
            if ind == -1:
                print('{} column is missing in input tsv file.'.format(name))
                exit(1)
        runs = []
        platforms = []
        models = []
        layouts = []
        scientific_names = []
        for line in f:
            temp = line.rstrip('\n').split('\t')
            runs.append(temp[run_ind])
            platforms.append(temp[platform_ind])
            models.append(temp[model_ind])
            layouts.append(temp[layout_ind])
            scientific_names.append(temp[scientific_name_ind])

        # check the amount of sra files in tsv
        if len(runs) > args.MaximumSRA:
            print('The amount of sra files is more than {}'.format(args.MaximumSRA))
            print('Randomly pick {} sra files for downloading'.format(args.MaximumSRA))
            runs_temp = []
            platforms_temp = []
            models_temp = []
            layouts_temp = []
            scientific_names_temp = []
            random_sra = random.sample(range(0,len(runs)), args.MaximumSRA)
            # randomly pick the maximum amount of sra files to download
            for i in random_sra:
                runs_temp.append(runs[i])
                platforms_temp.append(platforms[i])
                models_temp.append(models[i])
                layouts_temp.append(layouts[i])
                scientific_names_temp.append(scientific_names[i])
            runs = runs_temp
            platforms = platforms_temp
            models = models_temp
            layouts = layouts_temp
            scientific_names = scientific_names_temp
    # output runs/scientific_name/assembly_name to Source.txt
    date = datetime.datetime.now().strftime("%Y-%m-%d")
    with open(path.join(args.outdir, args.name, 'Source.txt'), 'w') as outfile:
        outfile.write(scientific_names[0] + '\n')
        outfile.write(args.assembly + '\n')
        outfile.write(date + '\n')
        for run in runs:
            outfile.write(run + '\n')

    single_files_for_merge = []
    paired1_files_for_merge = []
    paired2_files_for_merge = []
    for run, platform, model, layout in zip(runs, platforms, models, layouts):
        print('Processing the file: {}'.format(run))
        return_status, err_message, return_layout = run_pipeline(
            run=run,
            genome=genome,
            outdir=path.join(args.outdir, args.name),
            layout=layout,
            platform=platform,
            model=model,
            threads=args.threads_num,
            )
        if return_status:
            if return_layout == 'single':
                single_files_for_merge.append(
                    path.join(args.outdir, args.name, run, 'normalized.fastq'))
            if return_layout == 'paired':
                paired1_files_for_merge.append(
                    path.join(args.outdir, args.name, run, 'normalized_1.fastq'))
                paired2_files_for_merge.append(
                    path.join(args.outdir, args.name, run, 'normalized_2.fastq'))
        else:
            print(err_message)

    # merge normalized fastq files together
    if len(single_files_for_merge) > 0:
        with open(path.join(args.outdir, args.name, 'merged_normalized.fastq'), 'w') as outfile:
            for fname in single_files_for_merge:
                with open(fname, 'r') as infile:
                    for line in infile:
                        outfile.write(line)
    if len(paired1_files_for_merge) > 0 and len(paired2_files_for_merge) > 0:
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

    # HISAT2 index prefix (built once, shared by single and paired alignment)
    index_prefix = path.join(args.outdir, args.name, genome_file_name)

    # Aligning SINGLE fastq with HISAT2
    file_for_aligned = path.join(args.outdir, args.name, 'merged_normalized.fastq')
    if path.exists(file_for_aligned):
        build_hisat2_index(genome, index_prefix, args.threads_num)
        print('Aligning single file...')
        f_stdout = open(
            path.join(args.outdir, args.name, genome_file_name + '_single.hisat2.log'), 'w')
        f_stderr = open(
            path.join(args.outdir, args.name, genome_file_name + '_single.hisat2.errlog'), 'w')
        subprocess.run(
            [
                get_hisat2_command_path('hisat2'), '--dta-cufflinks', '--no-mixed', '--no-discordant',
                '-p', args.threads_num, '-x', index_prefix,
                '-U', file_for_aligned,
                '-S', path.join(args.outdir, args.name, 'single_output.sam')
            ],
            stdout=f_stdout,
            stderr=f_stderr)
        f_stdout.close()
        f_stderr.close()
    # Aligning PAIRED fastq with HISAT2
    file_for_aligned_1 = path.join(args.outdir, args.name, 'merged_normalized_1.fastq')
    file_for_aligned_2 = path.join(args.outdir, args.name, 'merged_normalized_2.fastq')
    if path.exists(file_for_aligned_1) and path.exists(file_for_aligned_2):
        build_hisat2_index(genome, index_prefix, args.threads_num)   # skipped if already built
        print('Aligning paired file...')
        f_stdout = open(
            path.join(args.outdir, args.name, genome_file_name + '_paired.hisat2.log'), 'w')
        f_stderr = open(
            path.join(args.outdir, args.name, genome_file_name + '_paired.hisat2.errlog'), 'w')
        subprocess.run(
            [
                get_hisat2_command_path('hisat2'), '--dta-cufflinks', '--no-mixed', '--no-discordant',
                '-p', args.threads_num, '-x', index_prefix,
                '-1', file_for_aligned_1,
                '-2', file_for_aligned_2,
                '-S', path.join(args.outdir, args.name, 'paired_output.sam')
            ],
            stdout=f_stdout,
            stderr=f_stderr)
        f_stdout.close()
        f_stderr.close()
    # sort and convert to the bam file-single
    if path.exists(path.join(args.outdir, args.name, 'single_output.sam')):
        f_stdout = open(
            path.join(args.outdir, args.name, genome_file_name + '_single.samtools.log'), 'w')
        f_stderr = open(
            path.join(args.outdir, args.name, genome_file_name + '_single.samtools.errlog'), 'w')
        subprocess.run(
            [
                'samtools', 'sort', '-@', args.threads_num, '-o',
                path.join(args.outdir, args.name, 'single_output.bam'), '-O', 'bam', '-T',
                path.join(args.outdir, args.name, 'single_output'),
                path.join(args.outdir, args.name, 'single_output.sam')
            ],
            stdout=f_stdout,
            stderr=f_stderr)
        f_stdout.close()
        f_stderr.close()
    # sort and convert to the bam file-paired
    if path.exists(path.join(args.outdir, args.name, 'paired_output.sam')):
        f_stdout = open(
            path.join(args.outdir, args.name, genome_file_name + '_paired.samtools.log'), 'w')
        f_stderr = open(
            path.join(args.outdir, args.name, genome_file_name + '_paired.samtools.errlog'), 'w')
        subprocess.run(
            [
                'samtools', 'sort', '-@', args.threads_num, '-o',
                path.join(args.outdir, args.name, 'paired_output.bam'), '-O', 'bam', '-T',
                path.join(args.outdir, args.name, 'paired_output'),
                path.join(args.outdir, args.name, 'paired_output.sam')
            ],
            stdout=f_stdout,
            stderr=f_stderr)
        f_stdout.close()
        f_stderr.close()
    # combine the sam/bam files together and convert to BAM file
    if path.exists(path.join(args.outdir, args.name, 'single_output.bam')) and path.exists(path.join(args.outdir, args.name, 'paired_output.bam')):
        print('Start merging single and paired output files...')
        files_for_merge = [path.join(args.outdir, args.name, 'single_output.bam'), path.join(args.outdir, args.name, 'paired_output.bam')]
        merge_files(files_for_merge, path.join(args.outdir, args.name))
    elif path.exists(path.join(args.outdir, args.name, 'single_output.bam')):
        print('The layout of all files is SINGLE, skip merging process...')
        os.rename(path.join(args.outdir, args.name, 'single_output.bam'), path.join(args.outdir, args.name, 'output.bam'))
    elif path.exists(path.join(args.outdir, args.name, 'paired_output.bam')):
        print('The layout of all files is PAIRED, skip merging process...')
        os.rename(path.join(args.outdir, args.name, 'paired_output.bam'), path.join(args.outdir, args.name, 'output.bam'))
    # output.bam is already coordinate-sorted; only sort if the header says otherwise
    final_bam = path.join(args.outdir, args.name, 'output.bam')
    ensure_coordinate_sorted(final_bam, args.threads_num)
    # index + bigwig (creates output.bam.csi and output.bigwig)
    print('Generating bigwig file...')
    make_bigwig(final_bam,
                path.join(args.outdir, args.name, 'output.bigwig'),
                args.threads_num)
    # rename bam and bigwig file to [gggsss]_[assembly_name]_RNA-Seq-alignments_[datetime]
    temp = scientific_names[0].split(" ")
    gene_name = temp[0]
    species_name = temp[1]
    if not args.Run_prefix:
        new_name = gene_name[0:3] + species_name[0:3]  + '_' + args.assembly + '_RNA-Seq-alignments_' + date
    if args.Run_prefix:
        new_name = runs[0]
    os.rename(final_bam, path.join(args.outdir, args.name, new_name + '.bam'))
    os.rename(final_bam + '.csi', path.join(args.outdir, args.name, new_name + '.bam.csi'))
    os.rename(path.join(args.outdir, args.name, 'output.bigwig'), path.join(args.outdir, args.name, new_name + '.bigwig'))
    # generate bed file
    print('Extracting junctions...')
    subprocess.run([get_regtools_path(), 'junctions', 'extract', '-m', '20', '-s', 'XS', '-o', path.join(args.outdir, args.name, new_name + '.bed'), path.join(args.outdir, args.name, new_name + '.bam')],
                   check=True)
    #  remove intermediate files
    # 1. Clean up individual run directories
    for run in runs:
        run_dir = path.join(args.outdir, args.name, run)
        # Remove ALL .fastq files (fixes the orphans/fixed_output_2 issue)
        for fastq_file in glob.glob(path.join(run_dir, "*.fastq")):
            os.remove(fastq_file)
        # Clean up FastQC zips and folders
        for fqc_item in glob.glob(path.join(run_dir, f"{run}_*_fastqc*")):
            if path.isdir(fqc_item):
                shutil.rmtree(fqc_item)
            elif path.exists(fqc_item):
                os.remove(fqc_item)
    # 2. Clean up top-level files if tempFile is False
    if not args.tempFile:
        base_dir = path.join(args.outdir, args.name)
        # List of specific files to remove if they exist
        top_level_files = [
            'single_output.sam', 'single_output.bam',
            'paired_output.sam', 'paired_output.bam',
            'merged_normalized_1.fastq', 'merged_normalized_2.fastq',
            'merged_normalized.fastq'
        ]
        for f in top_level_files:
            f_path = path.join(base_dir, f)
            if path.exists(f_path):
                os.remove(f_path)
        # Remove genome index files (e.g., HiSat2 .ht2 files)
        hs2_filelist = glob.glob(path.join(base_dir, genome_file_name + '.*'))
        for hs2_file in hs2_filelist:
            os.remove(hs2_file)
    print('Finished processing')
