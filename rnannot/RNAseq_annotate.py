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
import re


FASTQ_EXTS = ('.fastq.gz', '.fq.gz', '.fastq', '.fq')
# sample name, then _R1/_R2, _1/_2, .R1/.R2 or .1/.2, optionally followed by Illumina's _001
PAIR_RE = re.compile(r'^(?P<sample>.+?)[._]R?(?P<read>[12])(?:_001)?$', re.IGNORECASE)

# Rank words skipped when parsing scientific names (same rules as i5k_naming.py)
RANK_MARKERS = {'subsp.', 'subsp', 'ssp.', 'ssp', 'var.', 'var', 'subvar.',
                'f.', 'forma', 'cf.', 'aff.', 'x', '×'}


def parse_scientific_name(scientific_name):
    """Return [genus, species] or [genus, species, subspecies], cleaned for file names."""
    parts = [p for p in scientific_name.split() if p.lower() not in RANK_MARKERS]
    parts = [re.sub(r'[^A-Za-z0-9]', '', p) for p in parts]
    parts = [p for p in parts if p]
    if len(parts) < 2:
        raise ValueError('Cannot parse genus and species from "{}"'.format(scientific_name))
    if len(parts) > 3:
        parts = parts[:2] + [parts[-1]]   # keep the lowest infraspecific epithet
    return parts


def make_name_prefix(scientific_name):
    """Genus_species[_subspecies] from a scientific name, safe for filenames."""
    return '_'.join(parse_scientific_name(scientific_name))


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


def _cgroup_limit_bytes():
    """Smallest memory limit on this process's cgroup path (v2 or v1), or None.

    Slurm usually sets the limit on the job cgroup, not the step/task cgroup the
    process lives in, so walk from the process's cgroup up to the root.
    """
    try:
        with open('/proc/self/cgroup') as f:
            lines = f.read().splitlines()
    except OSError:
        return None
    limits = []
    for line in lines:
        parts = line.split(':', 2)
        if len(parts) != 3:
            continue
        hier, controllers, cg_path = parts
        if hier == '0' and controllers == '':                 # cgroup v2
            root, fname = '/sys/fs/cgroup', 'memory.max'
        elif 'memory' in controllers.split(','):              # cgroup v1
            root, fname = '/sys/fs/cgroup/memory', 'memory.limit_in_bytes'
        else:
            continue
        rel = cg_path.strip('/')
        while True:
            candidate = path.join(root, rel, fname) if rel else path.join(root, fname)
            try:
                with open(candidate) as f:
                    value = f.read().strip()
                if value != 'max' and int(value) < 2**60:     # huge value = no limit
                    limits.append(int(value))
            except (OSError, ValueError):
                pass
            if not rel:
                break
            rel = path.dirname(rel)
    return min(limits) if limits else None


def resolve_mem_gb(requested):
    """Usable memory in GB: the job's cgroup limit, the requested value, or physical RAM."""
    phys = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES') / 1024**3
    limit = _cgroup_limit_bytes()
    if limit:
        available = limit / 1024**3
    else:
        available = phys
    if requested is None or str(requested).strip().lower() in ('', '0', 'auto'):
        if not limit:
            print('WARNING: no cgroup memory limit found; assuming all {:.0f} GB of physical RAM. '
                  'Use -M to set the memory available to this run.'.format(phys))
        return available
    n = float(requested)
    if n > available:
        print('Requested {} GB but only {:.0f} GB are available; using {:.0f}'.format(
            n, available, available))
        n = available
    return n


def java_xmx(mem_gb, fraction=0.8):
    """-Xmx flag for BBTools: a fraction of the job's memory, leaving JVM/OS headroom."""
    return '-Xmx{}g'.format(max(1, int(mem_gb * fraction)))


def stage_local_fastq(src, dest):
    """Place a local FASTQ at dest: symlink plain files, decompress gzipped ones."""
    if src.endswith('.gz'):
        with gzip.open(src, 'rb') as f_in, open(dest, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    else:
        os.symlink(src, dest)   # cleanup removes the link, never the original


def discover_local_fastq(fastq_dir):
    """Scan fastq_dir and return [(sample, 'SINGLE'|'PAIRED', fq1, fq2_or_None), ...]."""
    fastq_dir = path.abspath(fastq_dir)
    if not path.isdir(fastq_dir):
        print('FASTQ directory not found: {}'.format(fastq_dir))
        exit(1)
    reads = {}      # sample -> {'1': path, '2': path, 'SE': path}
    errors = []
    for name in sorted(os.listdir(fastq_dir)):
        full = path.join(fastq_dir, name)
        ext = next((e for e in FASTQ_EXTS if name.lower().endswith(e)), None)
        if ext is None or not path.isfile(full):
            continue
        stem = name[:-len(ext)]
        m = PAIR_RE.match(stem)
        sample, read = (m.group('sample'), m.group('read')) if m else (stem, 'SE')
        slot = reads.setdefault(sample, {})
        if read in slot:
            errors.append('{}: more than one file for read {}: {} and {}'.format(
                sample, read, path.basename(slot[read]), name))
        slot[read] = full

    samples = []
    for sample, slot in sorted(reads.items()):
        if 'SE' in slot and ('1' in slot or '2' in slot):
            errors.append('{}: has both single-end and paired files'.format(sample))
        elif '1' in slot and '2' in slot:
            samples.append((sample, 'PAIRED', slot['1'], slot['2']))
        else:
            if '2' in slot:
                print('WARNING: {} has only a read-2 file; processing as single-end'.format(sample))
            only = slot.get('SE') or slot.get('1') or slot.get('2')
            samples.append((sample, 'SINGLE', only, None))
    if not samples:
        errors.append('No FASTQ files found in {}'.format(fastq_dir))
    if errors:
        print('Problems found in the FASTQ directory:')
        for e in errors:
            print('  ' + e)
        exit(1)

    print('Detected {} local sample(s):'.format(len(samples)))
    for sample, layout, fq1, fq2 in samples:
        print('  {:30s} {:7s} {}{}'.format(sample, layout, path.basename(fq1),
                                           ', ' + path.basename(fq2) if fq2 else ''))
    return samples


def run_pipeline(run, genome, outdir, layout, platform, model, threads, mem_gb, local_fastq=None):
    # create the output folder
    output_prefix = path.join(outdir, run)
    os.mkdir(output_prefix)
    if platform == 'ABI_SOLID':
        return (
            False,
            'Currently, the colorspace data from ABI_SOLID is not supported', '')
    sra_file_name = run
    if local_fastq is not None:
        print('staging local fastq files...')
        fq1, fq2 = local_fastq
        if layout == 'SINGLE':
            stage_local_fastq(fq1, path.join(output_prefix, sra_file_name + '.fastq'))
        else:
            stage_local_fastq(fq1, path.join(output_prefix, sra_file_name + '_1.fastq'))
            stage_local_fastq(fq2, path.join(output_prefix, sra_file_name + '_2.fastq'))
    else:
        # download SRA files and convert SRA file to fastq file(s)
        print('downloading sra files...')
        with open(path.join(output_prefix, sra_file_name + '.fastq-dump.log'), 'w') as f_stdout, \
             open(path.join(output_prefix, sra_file_name + '.fastq-dump.errlog'), 'w') as f_stderr:
            subprocess.run(['fasterq-dump', '-O', output_prefix, '-t', os.getcwd(), run],
                           stdout=f_stdout, stderr=f_stderr)
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
        # Element runs are listed under their own platform in SRA, so match them by model
        if (platform == 'ILLUMINA' and (model.startswith('Illumina HiSeq')
                                        or model.startswith('Illumina MiSeq'))) \
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
        # Normalizing (explicit -Xmx: BBTools would otherwise size the heap from the node's RAM)
        print('Normalizing ...')
        subprocess.run([get_bbmap_command_path('bbnorm.sh'), java_xmx(mem_gb), 'qin=33',
                        'in=' + path.join(output_prefix, 'output.fastq'),
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
        f_stdout.close()
        f_stderr.close()
        # Trimming with trimmomatic
        print('Trimming ...')
        f_stdout = open(
            path.join(output_prefix, sra_file_name + '.trimmomatic.log'), 'w')
        f_stderr = open(
            path.join(output_prefix, sra_file_name + '.trimmomatic.errlog'),
            'w')
        # Element runs are listed under their own platform in SRA, so match them by model
        if (platform == 'ILLUMINA' and (model.startswith('Illumina HiSeq')
                                        or model.startswith('Illumina MiSeq'))) \
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
                    get_bbmap_command_path('bbmerge.sh'), java_xmx(mem_gb), 'in1=' + path.join(
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
        # Normalizing (explicit -Xmx: BBTools would otherwise size the heap from the node's RAM)
        print('Normalizing ...')
        subprocess.run([get_bbmap_command_path('repair.sh'), java_xmx(mem_gb),
                        'in1=' + path.join(output_prefix, 'output_1.fastq'), 'in2=' + path.join(output_prefix, 'output_2.fastq'),
                        'out1=' + path.join(output_prefix, 'fixed_output_1.fastq'), 'out2=' + path.join(output_prefix, 'fixed_output_2.fastq'),
                        'outsingle=' + path.join(output_prefix, 'orphans.fastq')], check=True)
        subprocess.run([get_bbmap_command_path('bbnorm.sh'), java_xmx(mem_gb), 'qin=33',
                        'in1=' + path.join(output_prefix, 'fixed_output_1.fastq'), 'in2=' + path.join(output_prefix, 'fixed_output_2.fastq'),
                        'out1=' + path.join(output_prefix, 'normalized_1.fastq'), 'out2=' + path.join(output_prefix, 'normalized_2.fastq'),
                        'target=' + '100', 'prefilter=t', 'threads=' + threads], check=True)
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
    # use the CPUs and memory actually allocated to this job (Slurm/Apptainer aware)
    args.threads_num = resolve_threads(getattr(args, 'threads_num', None))
    print('Using {} threads'.format(args.threads_num))
    args.mem_gb = resolve_mem_gb(getattr(args, 'mem_gb', None))
    print('Using {:.0f} GB memory limit (BBTools {})'.format(args.mem_gb, java_xmx(args.mem_gb)))
    # convert many arguments to absolute path
    if not path.isabs(args.outdir):
        args.outdir = path.abspath(args.outdir)
    if args.input and not path.isabs(args.input):
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

    if args.fastq_dir:
        # Local mode: metadata comes from the command line, files from the directory
        samples = discover_local_fastq(args.fastq_dir)
        runs = [s[0] for s in samples]
        layouts = [s[1] for s in samples]
        local_files = [(s[2], s[3]) for s in samples]
        platforms = [args.platform] * len(runs)
        models = [args.model] * len(runs)
        scientific_names = [args.species] * len(runs)
    else:
        # SRA mode: metadata TSV from the NCBI pull
        with open(args.input) as f:
            col_names = f.readline().rstrip('\n').split('\t')
            required = ['Run', 'Platform', 'Model', 'LibraryLayout', 'ScientificName']
            missing = [c for c in required if c not in col_names]
            if missing:
                print('Missing column(s) in input tsv file: {}'.format(', '.join(missing)))
                exit(1)
            run_ind = col_names.index('Run')
            platform_ind = col_names.index('Platform')
            model_ind = col_names.index('Model')
            layout_ind = col_names.index('LibraryLayout')
            scientific_name_ind = col_names.index('ScientificName')
            print('Checking the input tsv file: {}'.format(args.input))
            runs = []
            platforms = []
            models = []
            layouts = []
            scientific_names = []
            for line in f:
                if not line.strip():
                    continue   # skip blank lines
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
            random_sra = random.sample(range(0, len(runs)), args.MaximumSRA)
            # randomly pick the maximum amount of sra files to download
            runs = [runs[i] for i in random_sra]
            platforms = [platforms[i] for i in random_sra]
            models = [models[i] for i in random_sra]
            layouts = [layouts[i] for i in random_sra]
            scientific_names = [scientific_names[i] for i in random_sra]
        local_files = [None] * len(runs)

    # output runs/scientific_name/assembly_name to Source.txt
    date = datetime.datetime.now().strftime("%Y-%m-%d")
    with open(path.join(args.outdir, args.name, 'Source.txt'), 'w') as outfile:
        outfile.write(scientific_names[0] + '\n')
        outfile.write(args.assembly + '\n')
        outfile.write(date + '\n')
        for run in runs:
            outfile.write(run + '\n')

    # Decide the output name now, so a bad name fails before any processing
    # [Genus_species(_subspecies)]_[assembly_name]_RNA-Seq-alignments_[date]
    distinct_names = sorted(set(scientific_names))
    if len(distinct_names) > 1:
        print('WARNING: runs have different scientific names: {}; using "{}" for naming'.format(
            '; '.join(distinct_names), scientific_names[0]))
    if args.Run_prefix:
        new_name = runs[0]
    else:
        if args.prefix:
            name_prefix = args.prefix
            print('Using supplied file name prefix: {}'.format(name_prefix))
        else:
            name_prefix = make_name_prefix(scientific_names[0])
            print('Derived file name prefix: {}'.format(name_prefix))
        new_name = name_prefix + '_' + args.assembly + '_RNA-Seq-alignments_' + date
    print('Output files will be named: {}.*'.format(new_name))

    single_files_for_merge = []
    paired1_files_for_merge = []
    paired2_files_for_merge = []
    for run, platform, model, layout, local_fastq in zip(runs, platforms, models, layouts, local_files):
        print('Processing the file: {}'.format(run))
        return_status, err_message, return_layout = run_pipeline(
            run=run,
            genome=genome,
            outdir=path.join(args.outdir, args.name),
            layout=layout,
            platform=platform,
            model=model,
            threads=args.threads_num,
            mem_gb=args.mem_gb,
            local_fastq=local_fastq,
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
    # rename to the name chosen above
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
        # Remove ALL .fastq files (fixes the orphans/fixed_output_2 issue);
        # for local plain FASTQs these are symlinks, so the originals are untouched
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
