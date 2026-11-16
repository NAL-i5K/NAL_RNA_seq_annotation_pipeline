# NAL RNA-Seq Annotation Pipeline
**Current release: 2026.11.13**

An RNA-Seq alignment pipeline that produces BAM, CSI index, bigWig coverage, and
splice-junction BED tracks for genome assemblies, and loads them into JBrowse/Apollo
for the i5k Workspace. It is built on
[SRA Toolkit](https://github.com/ncbi/sra-tools),
[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/),
[Trimmomatic](https://github.com/usadellab/Trimmomatic),
[BBMap](https://sourceforge.net/projects/bbmap/) (BBNorm, BBMerge, repair.sh),
[HISAT2](https://github.com/DaehwanKimLab/hisat2),
[samtools](https://github.com/samtools/samtools),
[Picard](https://broadinstitute.github.io/picard/),
[regtools](https://github.com/griffithlab/regtools),
[RSEM](https://github.com/deweylab/RSEM) (`rsem-bam2wig`),
UCSC `wigToBigWig`, and
[bam_to_bigwig](https://github.com/NAL-i5K/bam_to_bigwig).

The pipeline is distributed as a Docker image, which is converted to an Apptainer
image for HPC use. Releases use calendar versions (`YYYY.MM.DD`).

## Contents

- [Overview](#overview)
- [Output modes: merged or per-sample](#output-modes-merged-or-per-sample)
- [Naming](#naming)
- [Container contents](#container-contents)
- [Installation](#installation)
- [Usage](#usage)
- [Running on SCINet with Apptainer and SLURM](#running-on-scinet-with-apptainer-and-slurm)
- [Loading into JBrowse/Apollo](#loading-into-jbrowseapollo)
- [Examples](#examples)
- [Troubleshooting](#troubleshooting)
- [Notes](#notes)
- [Tests](#tests)

## Overview

The pipeline has two parts.

**1. rnannot (alignment), run on HPC (e.g. SCINet Ceres/Atlas)**

1. Get reads: download SRA runs (`fasterq-dump`) or use a directory of local FASTQ files
2. FastQC
3. Adapter and quality trimming (Trimmomatic)
4. Read normalization (BBNorm, `target=100`)
5. Build the HISAT2 index (once per run, multithreaded)
6. Align with HISAT2 (`--dta-cufflinks`, which writes XS strand tags)
7. Sort, and merge single- and paired-end alignments (merged mode)
8. CSI index, bigWig coverage (`bam_to_bigwig.py`), and junction BED
   (`regtools junctions extract -m 20 -s XS`)

**2. JBrowse loading, run on the Apollo server**

- `add_trackList.py` copies the outputs into the species directory and adds the BAM,
  bigWig, and junction tracks to `trackList.json`.
- `move_data.py` (legacy) syncs loaded data to apollo-node1.
- Both use the shared naming rules in `i5k_naming.py`, which must be importable by the
  two scripts (see [Troubleshooting](#troubleshooting)).

## Output modes: merged or per-sample

### Merged (default)

All runs or samples are combined into **one** set of outputs in `<outdir>/<name>/`:

| File | Description |
|---|---|
| `<prefix>_<assembly>_RNA-Seq-alignments_<date>.bam` | Coordinate-sorted alignments |
| `<prefix>_<assembly>_RNA-Seq-alignments_<date>.bam.csi` | CSI index |
| `<prefix>_<assembly>_RNA-Seq-alignments_<date>.bigwig` | Coverage track |
| `<prefix>_<assembly>_RNA-Seq-alignments_<date>.bed` | Splice junctions |
| `Source.txt` | Scientific name, assembly name, date, and list of runs/samples |
| `<run>/` | Per-run logs and FastQC reports |

### Per-sample (`--per_sample`, local mode only)

Each sample in `-d` (for example, each tissue) gets its **own** BAM, CSI, bigWig, and
BED in its own folder:

```
<outdir>/<name>/
├── Source.txt                                   # all samples in this run
├── Zacosmia_maculata_Abdomen/
│   ├── Zacosmia_maculata_Abdomen_<assembly>_RNA-Seq-alignments_<date>.bam
│   ├── ...bam.csi  ...bigwig  ...bed
│   ├── Source.txt                               # this sample only (use for loading)
│   └── *.log, *.errlog, FastQC reports
├── Zacosmia_maculata_Head/
└── Zacosmia_maculata_Thorax/
```

How it works:

- The HISAT2 index is built once and shared by all samples in the job.
- Samples are processed one at a time. Each sample's FASTQ copies are deleted before
  the next sample starts, so peak disk use is about one sample's worth of
  intermediates.
- If one sample fails, the others still run. The job then exits with an error and
  lists the failed samples.
- `--sample NAME` processes only the named sample(s). It can be repeated, and it implies
  `--per_sample`. This is useful for running samples in parallel as a Slurm array.

`<date>` is the day the pipeline runs. With `-t/--tempFile`, intermediate SAM/BAM
files, merged normalized FASTQs, the HISAT2 index, and the decompressed genome are kept.

## Naming

### Output file names

Choose **at most one** of these options. Combining them is an error.

| Option | Merged mode | Per-sample mode |
|---|---|---|
| *(none)* | `Genus_species[_subspecies]_<assembly>_RNA-Seq-alignments_<date>` | `Genus_species[_subspecies]_<label>_<assembly>_RNA-Seq-alignments_<date>` |
| `--prefix NAME` | `NAME_<assembly>_RNA-Seq-alignments_<date>` | `NAME_<label>_<assembly>_RNA-Seq-alignments_<date>` |
| `--input_prefix` | `<prefix from input names>_<assembly>_RNA-Seq-alignments_<date>` | `<prefix from input names>_<label>_<assembly>_RNA-Seq-alignments_<date>` |
| `-R, --Run_prefix` | First run ID or sample name, used as the whole name | Each sample name, used as the whole name |

In per-sample mode, `<label>` is the sample name with the leading prefix or species
name removed. For example, `Zacosmia_maculata_Head` gives the label `Head`. If nothing
is left (one sample whose name is the prefix), no label is added. If two samples would
get the same name, the run stops before processing.

All output names are printed at the start of the run, before any processing. Check
them in the job log.

#### From the scientific name (default)

| Scientific name | File prefix |
|---|---|
| `Apis mellifera` | `Apis_mellifera` |
| `Apis mellifera ligustica` | `Apis_mellifera_ligustica` |
| `Brassica rapa subsp. pekinensis` | `Brassica_rapa_pekinensis` (rank words such as `subsp.` and `var.` are dropped) |
| `Osmia lignaria_lignaria` | `Osmia_lignaria_lignaria` (underscores are read as spaces) |

- Case is not changed. SRA metadata is already correctly capitalized. With `-d`, type
  the name properly, e.g. `-s "Apis mellifera"`.
- Underscores and extra whitespace in the name are converted to single spaces before
  parsing, and the cleaned name is what goes into `Source.txt`.
- Names with a strain or cultivar label keep it as the last part
  (`Saccharomyces cerevisiae S288C` → `Saccharomyces_cerevisiae_S288C`). Use `--prefix`
  to choose something else.

#### `--prefix`

Used as given. Only letters, digits, and underscores are allowed. Case is kept.

#### `--input_prefix`

- **Local mode (`-d`):** the prefix is the leading part shared by all sample names.
  Sample names are the file names with the extension and the read suffix
  (`_R1`, `_2`, `_R1_001`, ...) removed. They are compared as whole words split on `_`,
  `.`, and `-`, so the prefix never ends in the middle of a word. With only one sample
  (including a single `--sample`), the prefix is the whole sample name.
- **SRA mode (`-i`):** the prefix is the TSV file name without `.tsv`, `.txt`, `.csv`,
  or `.gz`.
- Characters other than letters, digits, and underscores become `_`.
- If the sample names share no leading word, the run stops and asks for `--prefix`.

| Input | Prefix |
|---|---|
| `Zacosmia_maculata_{Abdomen,Head,Thorax}_R{1,2}.fastq.gz` | `Zacosmia_maculata` |
| `sampleA_S1_L001_R1_001.fastq.gz` … `sampleA_S1_L002_R2_001.fastq.gz` | `sampleA_S1` |
| A single sample, `sampleB.fastq.gz` | `sampleB` |
| `Ctrl_rep1_R1.fq.gz` … and `Heat_rep1_R1.fq.gz` … | Stops: no common prefix |
| `-i Apis_runs.tsv` | `Apis_runs` |
| `-i 1049336.tsv` | `1049336` (usually not what you want; use the default or `--prefix`) |

When the files already start with `Genus_species`, as in the Zacosmia example,
`--input_prefix` gives the same result as the default.

### JBrowse species directories (loaders)

The i5k directory convention is changing. Legacy species keep their `gggsss`
directories, and new species use full names. The loaders detect which applies:

| Situation | Species directory |
|---|---|
| A legacy directory already exists (`apimel`, or `apmeli` for a trinomial) | The legacy directory |
| New species | `Genus_species[_subspecies]` |
| `--species_dir` given | Exactly that |
| Both legacy and full-name directories exist | Stops; choose with `--species_dir` |
| More than one legacy code matches | Stops; choose with `--species_dir` |

A trinomial is checked against both legacy codes. For example,
`Osmia lignaria lignaria` checks `osmlig` and `oslili`. If a directory for the
binomial already exists (`osmlig`), the subspecies is loaded there unless you pass
`--species_dir`.

`--species_dir` (aliases `-g`, `--gggsss`) accepts letters, digits, and underscores,
so both `apimel` and `Saccharomyces_eubayanus` are valid. If an underscore is rejected
with "must contain only letters and digits", you are running a pre-2026.11.13 copy of
the loader scripts. Replace `add_trackList.py`, `move_data.py`, and `i5k_naming.py`
with the current versions.

The analyses folder is `<species_dir>-RNA-Seq_v1.0`, shared by all RNA-Seq datasets for
an assembly. The dated file names, plus the sample label in per-sample mode, keep
datasets distinct.

**Always run `add_trackList.py --dry_run` the first time you load a species.** The
script prints which directory it chose and why.

## Container contents

Tool versions are set as `ARG`s at the top of the `Dockerfile`.

| Component | Version |
|---|---|
| Base image | Ubuntu 24.04 |
| Python | 3.12 (virtual environment at `/opt/venv`) |
| Java | OpenJDK 21 (`default-jre-headless`) |
| samtools / pysam | 1.24 / 0.24.1 |
| SRA Toolkit | 3.4.1 |
| FastQC | 0.12.1 |
| Trimmomatic | 0.39 |
| HISAT2 | 2.2.1 |
| BBMap | 39.33 |
| Picard | 3.5.0 |
| regtools | 1.0.0 |
| RSEM, NCBI EDirect | Ubuntu 24.04 packages |
| wigToBigWig | UCSC binary, installed with rsync from `hgdownload.soe.ucsc.edu` (SHA-256 recorded in the release notes) |
| bam_to_bigwig | pinned commit (`BAM_TO_BIGWIG_COMMIT`) |

The exact versions in each release are listed in `release-<version>.txt`.

## Installation

### Build the Docker image

```shell
git clone https://github.com/NAL-i5K/NAL_RNA_seq_annotation_pipeline.git
cd NAL_RNA_seq_annotation_pipeline
docker build --pull --no-cache --progress=plain -t nal-rnaseq:2026.11.13 . 2>&1 | tee build-2026.11.13.log
```

The build needs outbound TCP port 873 (rsync) for the wigToBigWig download.

### Convert to Apptainer

On the machine where the image was built:

```shell
apptainer build nal-rnaseq_2026.11.13.sif docker-daemon:nal-rnaseq:2026.11.13
```

Building on one machine and running on another:

```shell
docker save nal-rnaseq:2026.11.13 -o nal-rnaseq_2026.11.13.tar
# copy the .tar to the cluster, then:
apptainer build nal-rnaseq_2026.11.13.sif docker-archive:nal-rnaseq_2026.11.13.tar
```

After **any** code change, rebuild the Docker image **and** reconvert the SIF. Then
confirm that the new code is inside the SIF (see [Tests](#tests)).

### Python package only (advanced)

```shell
pip install .        # or: pip install -e .  so edits take effect without reinstalling
```

All external tools must then be installed separately. Set `RNANNOT_LIB` to the
directory containing them (the same layout as `rnannot/lib` in the container), and put
`samtools`, `fasterq-dump`, `rsem-bam2wig`, and `wigToBigWig` on `PATH`.

## Usage

### Get SRA metadata

```shell
download_sra_metadata.py [-t TAXID] [-o OUTPUT]
  -t, --taxid TAXID        find all RNA-Seq SRA runs for a taxid
  -o, --output OUTPUT      output tsv (default: current folder)

download_sra_metadata_by_accessions.py [-a ACCESSION ...] [-o OUTPUT ...]
  -a, --accession          one or more SRA accessions (space-separated)
  -o, --output             output file for each accession (default: ACCESSION.tsv)
```

The output tsv is the `-i` input for `RNAseq_annotate.py`.

### Run the pipeline

```shell
RNAseq_annotate.py (-i INPUT | -d FASTQ_DIR -s "Genus species") -g GENOME
                   [--per_sample | --sample NAME ...]
                   [-R | --prefix PREFIX | --input_prefix] [options]
```

**Input**

| Option | Description |
|---|---|
| `-i, --input` | SRA metadata tsv (from `download_sra_metadata*.py`) |
| `-d, --fastq_dir` | Directory of local FASTQ files. Use instead of `-i`. Requires `-s` |
| `-s, --species` | Scientific name, e.g. `"Apis mellifera"` or `"Apis mellifera ligustica"`. Required with `-d`, even with `--input_prefix` (it is written to `Source.txt`) |
| `--platform` | Platform for local files. Default `ILLUMINA` |
| `--model` | Instrument model for local files, e.g. `"Illumina HiSeq 2500"` or `"Element AVITI"`. Selects the adapter set |
| `-g, --genome` | Genome FASTA (`.gz` accepted). Required |

**Per-sample outputs (local mode only)**

| Option | Description |
|---|---|
| `--per_sample` | Align each sample in `-d` separately and write a BAM, CSI, bigWig, and BED for each |
| `--sample NAME` | Process only this sample (repeatable). Implies `--per_sample`. Names are as detected from `-d`, without the read suffix, e.g. `Zacosmia_maculata_Head` |

**Output location and resources**

| Option | Description |
|---|---|
| `-o, --outdir` | Parent directory for the output folder (must already exist). Default: current folder |
| `-n, --name` | Output folder name (must **not** already exist). Default: start time |
| `-a, --assembly` | Assembly name, used in output file names. Default `AssemblyName` |
| `-t, --tempFile` | Keep intermediate files |
| `-m, --MaximumSRA` | SRA mode only: randomly select at most this many runs. Default 10 |
| `-T, --threads_num` | Number of threads, or `auto` (default) to use all CPUs allocated to the job. Values above the allocation are capped |
| `-M, --mem_gb` | Memory in GB, or `auto` (default) to use the job's Slurm/cgroup memory limit. Values above the limit are capped |

**Output naming** (use at most one; see [Naming](#naming))

| Option | Description |
|---|---|
| `--prefix` | Output file prefix (letters, digits, underscores) |
| `--input_prefix` | Take the prefix from the input file names: the common start of the sample names with `-d`, or the TSV file name with `-i` |
| `-R, --Run_prefix` | Use the run ID or sample name as the whole output name |

Exactly one of `-i` or `-d` must be given.

### What happens at startup

Before any processing, the pipeline:

1. Prints the path of the script being run (`RNAseq_annotate.py: /opt/...`), the
   thread count, the memory limit and its source.
2. Stops if `<outdir>/<name>/` already exists.
3. Reads the TSV or scans the FASTQ directory, and stops on input errors.
4. Decides and prints the output name(s), and stops on naming errors.
5. Only then creates `<outdir>/<name>/` and `Source.txt`.

An input or naming error therefore leaves nothing behind, so you can fix it and rerun
with the same `-n`.

### Input: SRA metadata tsv

Tab-separated, with a header row containing at least:

| Column | Use |
|---|---|
| `Run` | SRA run accession |
| `Platform` | e.g. `ILLUMINA`, `ELEMENT`. `ABI_SOLID` (colorspace) is skipped |
| `Model` | Instrument model. Selects the adapter set |
| `LibraryLayout` | `SINGLE` or `PAIRED` |
| `ScientificName` | Used for output names and `Source.txt` |

If the runs have different scientific names (e.g. several subspecies from a taxid
pull), the pipeline prints a warning and uses the first name. Use `--prefix` to set
the name yourself. SRA mode always produces merged outputs.

### Input: local FASTQ directory

All `.fastq`, `.fq`, `.fastq.gz`, and `.fq.gz` files directly in the directory are
used. Subdirectories are not searched. Pairs are detected from file names:

| File names | Detected as |
|---|---|
| `Zacosmia_maculata_Head_R1.fastq.gz` + `..._R2.fastq.gz` | PAIRED, sample `Zacosmia_maculata_Head` |
| `sampleA_S1_L001_R1_001.fastq.gz` + `..._R2_001.fastq.gz` | PAIRED, sample `sampleA_S1_L001` |
| `sampleA_R1.fq.gz` + `sampleA_R2.fq.gz` | PAIRED, sample `sampleA` |
| `sampleA_1.fastq` + `sampleA_2.fastq` | PAIRED, sample `sampleA` |
| `sampleB.fastq.gz` | SINGLE, sample `sampleB` |

- Detected samples are printed at startup. The run stops before processing if there
  are duplicate read files or a sample has both single and paired files.
- A sample with only a read-2 file is processed as single-end, with a warning.
- Original files are only read, never modified or deleted, so read-only directories
  work. Plain FASTQs are symlinked. Gzipped FASTQs are decompressed into the sample
  folder and removed at cleanup.
- In merged mode, multi-lane samples are processed separately and merged before
  alignment. In per-sample mode, each lane is its own sample.
- `-m/--MaximumSRA` does not apply. Every file (or every `--sample`) is used.

Preview the samples and output names without running the pipeline:

```shell
apptainer exec --cleanenv --bind "$FASTQ" nal-rnaseq_2026.11.13.sif python3 -c "
from rnannot.RNAseq_annotate import discover_local_fastq, make_input_name_prefix, make_name_prefix, sample_label
samples = [s[0] for s in discover_local_fastq('$FASTQ')]
sp = make_name_prefix('$SPECIES')
print('--input_prefix gives:', make_input_name_prefix(samples))
print('per-sample labels:', {s: sample_label(s, [sp]) for s in samples})"
```

### Adapter selection

| Platform / Model | Adapters |
|---|---|
| `ILLUMINA` + `Illumina HiSeq*` or `Illumina MiSeq*`; any `Element AVITI*` model | TruSeq3 |
| `ILLUMINA` + `Illumina Genome Analyzer II*` | TruSeq2 |
| Anything else (including empty `--model`) | Paired: detected with BBMerge. Single: BBMap `adapters.fa` |

## Running on SCINet with Apptainer and SLURM

Paths used inside the container must be bound with `--bind`. Use `readlink -f` so
symlinked project paths resolve correctly.

### One job (merged, or all samples one after another)

```shell
#!/bin/bash
#SBATCH --job-name=rnannot
#SBATCH --cpus-per-task=16        # at least 12: bam_to_bigwig.py uses 12 threads internally
#SBATCH --mem=64G
#SBATCH --time=48:00:00           # per-sample mode: allow time for every sample
#SBATCH --output=rnannot_%j.out   # stdout and stderr go to this file if --error is not set

set -euo pipefail                         # stop the script if a command fails
export APPTAINERENV_PYTHONUNBUFFERED=1    # keep log output in order (works with --cleanenv)

module load apptainer
SIF=$(readlink -f /project/nal_genomics/tools/nal-rnaseq/nal-rnaseq_2026.11.13.sif)
WORK=$(readlink -f /project/nal_genomics/$USER/rnaseq/zacosmia)
GENOME=$(readlink -f /path/to/genome.fna.gz)

# --- SRA mode (merged) ---
# apptainer exec --cleanenv --bind "$WORK" --bind "$(dirname "$GENOME")" --pwd "$WORK" "$SIF" \
#     RNAseq_annotate.py -i "$WORK/7460.tsv" -g "$GENOME" -a Amel_HAv3.1 -o "$WORK"

# --- Local mode, one output set per tissue ---
FASTQ=$(readlink -f /path/to/raw_reads/Zacosmia_maculata)
apptainer exec --cleanenv --bind "$WORK" --bind "$FASTQ" --bind "$(dirname "$GENOME")" \
    --pwd "$WORK" "$SIF" \
    RNAseq_annotate.py -d "$FASTQ" -s "Zacosmia maculata" -g "$GENOME" -a "<assembly>" \
        --per_sample -o "$WORK" -n Zacosmia_maculata_rnaseq
```

Remove `--per_sample` to merge all samples into one output set instead.

### Slurm array (samples in parallel)

Each task processes one sample from the same directory. Each task needs its own `-n`
and builds its own HISAT2 index.

```shell
#SBATCH --array=0-2               # one task per sample (3 here)
# ... same header and variables as above ...

mapfile -t SAMPLES < <(ls "$FASTQ" | sed -nE 's/_R1\.fastq\.gz$//p' | sort)
S=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

apptainer exec --cleanenv --bind "$WORK" --bind "$FASTQ" --bind "$(dirname "$GENOME")" \
    --pwd "$WORK" "$SIF" \
    RNAseq_annotate.py -d "$FASTQ" -s "Zacosmia maculata" -g "$GENOME" -a "<assembly>" \
        --sample "$S" -o "$WORK" -n "$S"
```

The `sed` pattern assumes `_R1.fastq.gz`. Adjust it for `_R1_001.fastq.gz`, `.fq.gz`,
and so on. Use `--array=0-2%1` to run one task at a time if disk space is limited.

### Threads (`-T`) and memory (`-M`)

- By default the pipeline uses every CPU and all the memory Slurm allocated to the job.
  CPUs come from the process's CPU affinity. Memory is the smaller of the Slurm
  allocation and the cgroup memory limit. When `--cleanenv` removes the Slurm
  variables, they are read from the parent batch shell.
- The job log shows `Using N threads`, `Memory limit source: ...`, and
  `Using N GB memory limit (BBTools -Xmx...)`. These should match `--cpus-per-task`
  and `--mem`.
- BBTools (BBNorm, BBMerge, repair.sh) get an explicit `-Xmx` of 80% of the limit.
  Trimmomatic gets 80% of at most 10 GB. `samtools sort` memory per thread is set so
  the total stays within half the limit (100 MB–4 GB per thread).
- Inside a Slurm job, if the memory limit can't be determined, the run stops and asks
  for `-M`. Outside Slurm, all physical RAM is assumed, with a warning.
- Always set `--cpus-per-task`. Without it, Slurm allocates one CPU and the pipeline
  runs single-threaded.
- If you wrap the command in `srun`, repeat the CPU request:
  `srun --cpus-per-task="$SLURM_CPUS_PER_TASK" apptainer exec ...`

### Monitoring

- The first log line shows which `RNAseq_annotate.py` is running. Check it after
  rebuilding the SIF.
- Failures in alignment, sorting, indexing, bigWig, or regtools stop the run (or, in
  per-sample mode, that sample) with an error in the job log.
- Per-sample logs (`hisat2.errlog`, `samtools.errlog`, `*.trimmomatic.errlog`, ...)
  are in each sample folder.
- Check resource use after a run:
  `sacct -j <jobid> --format=JobID,State,ExitCode,Elapsed,MaxRSS,ReqMem` or `seff <jobid>`.
- A plain `hisat2-build` needs roughly 8 GB of RAM for a human-sized genome.

## Loading into JBrowse/Apollo

### add_trackList.py

```shell
add_trackList.py -op OUTPUT_PATH -j JBROWSE_SIF -bam BAM -bigwig BIGWIG \
                 -csi CSI -bed BED -track TRACKLIST -s SOURCE \
                 [--species_dir DIR] [--prefix PREFIX] [--overwrite] [--dry_run]
```

| Option | Description |
|---|---|
| `-op, --output_path` | Species data root, e.g. `/app/data/other_species` |
| `-j, --jbrowse_path` | JBrowse container (run with `singularity exec`) |
| `-bam`, `-bigwig`, `-csi`, `-bed` | Pipeline output files |
| `-track, --input_track` | `trackList.json` path (a `.bak` backup is written first) |
| `-s, --Source` | `Source.txt` from the pipeline run. In per-sample mode, use the one in the sample folder |
| `--species_dir` (aliases `-g`, `--gggsss`) | Species directory (letters, digits, underscores). Default: existing legacy directory, otherwise full name |
| `--prefix` | File prefix used by the pipeline. Give the value passed to the pipeline's `--prefix`, or, after an `--input_prefix` run, the prefix it printed. Not needed when outputs use the default species prefix |
| `--overwrite` | Replace files that already exist in the analyses folder |
| `--dry_run` | Print the resolved paths and exit without copying or editing anything |

Files are copied to
`<output_path>/<species_dir>/<assembly>/scaffold/analyses/<species_dir>-RNA-Seq_v1.0/`.
Track labels are the BAM name without `.bam`, plus `_coverage` and `_junctions` for the
bigWig and BED tracks. The script warns if the BAM name doesn't start with the expected
prefix.

**Per-sample outputs** are loaded once per sample:

```shell
RUN=$WORK/Zacosmia_maculata_rnaseq
for D in "$RUN"/*/; do
  B=$(ls "$D"*.bam 2>/dev/null) || continue
  add_trackList.py -op /app/data/other_species -j /path/to/jbrowse.sif \
    -bam "$B" -csi "$B.csi" -bigwig "${B%.bam}.bigwig" -bed "${B%.bam}.bed" \
    -track <trackList.json> -s "$D/Source.txt" --dry_run
done
```

Each tissue gets its own BAM, coverage, and junction tracks in the shared analyses
folder. Remove `--dry_run` once the paths look right.

### move_data.py (legacy)

```shell
move_data.py -Node1a user@apollo-node1.nal.usda.gov -s Source.txt \
             [--species_dir DIR] [--prefix PREFIX] [--label LABEL] \
             [--base_path /app/data/other_species] [--dry_run]
```

Syncs the analyses folder and junction track to apollo-node1, creates the
`jbrowse/data/analyses` symlink if it is missing, and copies `trackList.json` last.
Use the same `--species_dir` and `--prefix` values given to `add_trackList.py`. Use
`--label` only if the pipeline ran with `-R/--Run_prefix`.

### Legacy and new species

- Legacy species loaded before this release keep working. Their directories, dated
  analyses folders, and `trackList.json` entries are untouched.
- New species directories must be registered in Apollo under their new paths. That
  happens outside these scripts.
- To migrate a legacy species, rename its directory and update the `urlTemplate` /
  `csiUrlTemplate` paths in `trackList.json`. Once only the full-name directory exists,
  the loaders pick it automatically.

## Examples

### SRA mode

```shell
# 1. Metadata
download_sra_metadata.py -t 1049336 -o 1049336.tsv

# 2. Genome
wget "https://i5k.nal.usda.gov/data/Arthropoda/ephdan-(Ephemera_danica)/BCM-After-Atlas/1.Genome%20Assembly/BCM-After-Atlas/Scaffolds/Edan07162013.scaffolds.fa.gz"

# 3. Pipeline (inside the container, see the SLURM example above)
RNAseq_annotate.py -i 1049336.tsv -g Edan07162013.scaffolds.fa.gz -a Edan_2.0
#   -> Ephemera_danica_Edan_2.0_RNA-Seq-alignments_<date>.bam / .bam.csi / .bigwig / .bed

# 4. Check where the loader will put it (legacy directory ephdan exists, so it is used)
add_trackList.py -op /app/data/other_species -j /path/to/jbrowse.sif \
  -bam Ephemera_danica_Edan_2.0_RNA-Seq-alignments_<date>.bam \
  -csi Ephemera_danica_Edan_2.0_RNA-Seq-alignments_<date>.bam.csi \
  -bigwig Ephemera_danica_Edan_2.0_RNA-Seq-alignments_<date>.bigwig \
  -bed Ephemera_danica_Edan_2.0_RNA-Seq-alignments_<date>.bed \
  -track /app/data/other_species/ephdan/Edan_2.0/jbrowse/data/trackList.json \
  -s Source.txt --dry_run

# 5. Load (same command without --dry_run)
```

### Local FASTQ, one output set per tissue

```shell
# $FASTQ contains:
#   Zacosmia_maculata_Abdomen_R1.fastq.gz  Zacosmia_maculata_Abdomen_R2.fastq.gz
#   Zacosmia_maculata_Head_R1.fastq.gz     Zacosmia_maculata_Head_R2.fastq.gz
#   Zacosmia_maculata_Thorax_R1.fastq.gz   Zacosmia_maculata_Thorax_R2.fastq.gz
RNAseq_annotate.py -d "$FASTQ" -s "Zacosmia maculata" -g genome.fa.gz \
    -a <assembly> --per_sample -o "$WORK" -n Zacosmia_maculata_rnaseq
```

The log should show, within seconds:

```
Detected 3 local sample(s):
  Zacosmia_maculata_Abdomen      PAIRED  ...
Derived file name prefix: Zacosmia_maculata
Per-sample mode: 3 sample(s), each with its own BAM/CSI/bigWig/BED
Per-sample output names:
  Zacosmia_maculata_Abdomen -> Zacosmia_maculata_Abdomen/Zacosmia_maculata_Abdomen_<assembly>_RNA-Seq-alignments_<date>.*
  Zacosmia_maculata_Head    -> Zacosmia_maculata_Head/Zacosmia_maculata_Head_<assembly>_RNA-Seq-alignments_<date>.*
  Zacosmia_maculata_Thorax  -> Zacosmia_maculata_Thorax/Zacosmia_maculata_Thorax_<assembly>_RNA-Seq-alignments_<date>.*
Building HISAT2 index...
```

### Local FASTQ, all tissues merged

```shell
RNAseq_annotate.py -d "$FASTQ" -s "Zacosmia maculata" -g genome.fa.gz \
    -a <assembly> -o "$WORK" -n Zacosmia_maculata_merged
#   -> Zacosmia_maculata_<assembly>_RNA-Seq-alignments_<date>.*
```

## Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| `AttributeError: 'Namespace' object has no attribute 'input_prefix'` (or `per_sample`) | An older `parser.py` is being imported | Replace `rnannot/parser.py`, rebuild the image and SIF, and check `RNAseq_annotate.py --help` |
| Run ends seconds after the naming lines, with no `Per-sample mode:` line | An older `RNAseq_annotate.py` is in the SIF, or `--per_sample` isn't on the command line | Check the first log line (script path) and the job script's command line |
| `ERROR: output folder already exists` | `-n` folder from an earlier run | Remove `<outdir>/<name>/` or choose another `-n` |
| `argument ...: not allowed with argument ...` | More than one naming option given | Use only one of `-R`, `--prefix`, `--input_prefix` |
| `--per_sample/--sample require -d/--fastq_dir` | Per-sample options with `-i` | Per-sample mode is local-only |
| `Sample(s) not found` | `--sample` name doesn't match a detected sample | Use a name from the `Detected samples:` list (no read suffix or extension) |
| `input sample names have no common prefix` | `--input_prefix` with unrelated sample names | Use `--prefix` |
| Loader warns `bam file ... does not start with "Osmia_lignarialignaria_"` | Line 1 of an older `Source.txt` has an underscore (`Osmia lignaria_lignaria`) | Fix line 1 to `Osmia lignaria lignaria`, or update `i5k_naming.py` to read underscores as spaces |
| Log lines out of order, or a traceback above earlier output | Python output buffering in a combined Slurm log | Set `APPTAINERENV_PYTHONUNBUFFERED=1` |
| Job reports success but nothing was produced | The job script kept going after Python failed | Add `set -euo pipefail`; check `sacct` ExitCode |
| Loader imports the wrong `i5k_naming.py` | The loader adds its parent's parent directory to `sys.path` | Print `i5k_naming.__file__` to check which copy is loaded; remove stale copies |

## Notes

- **SRA Toolkit configuration.** Current SRA Toolkit releases don't need interactive
  configuration. If downloads fail, run `vdb-config --restore-defaults`.
- **Disk space.** `fasterq-dump` needs roughly 17× the accession size during
  conversion. Gzipped local FASTQs are decompressed, and trimmed, repaired, and
  normalized copies are written alongside them. The pipeline also writes uncompressed
  SAM files before sorting. Plan for several times the compressed input size; per-sample
  mode needs this for one sample at a time. Run in a project or scratch directory,
  not `$HOME`.
- **Normalization.** Reads are normalized with BBNorm before alignment. The bigWig
  shows normalized depth for annotation evidence, not expression levels, so coverage
  is not comparable between tissues as a measure of expression.
- **Local samples in JBrowse metadata.** The "Data source" field links each entry in
  `Source.txt` to an SRA search. For local-mode runs these are sample names, so the
  links won't resolve.
- **Older `Source.txt` files.** Files created by hand before subspecies support may
  have the scientific name in file-name form (e.g. `Osmia lignaria_lignaria`). Fix
  line 1 before loading.
- **Legacy Docker Hub image.** `k2025242322/i5k_rna_seq_annotation_pipeline` predates
  this release. Use the image built from this repository.

## Tests

```shell
python -m unittest -f tests/test_parser.py
```

Quick checks after rebuilding the SIF:

```shell
# Options are present
apptainer exec --cleanenv nal-rnaseq_2026.11.13.sif RNAseq_annotate.py --help \
  | grep -E 'input_prefix|per_sample|--sample'

# The right code is loaded
apptainer exec --cleanenv nal-rnaseq_2026.11.13.sif python3 -c "
import rnannot.parser as p, rnannot.RNAseq_annotate as r
print(p.__file__); print(r.__file__)
print(all(hasattr(r, f) for f in ['make_input_name_prefix', 'sample_label', 'align_sample']))"

# Naming rules
apptainer exec --cleanenv nal-rnaseq_2026.11.13.sif python3 -c "
from rnannot.RNAseq_annotate import make_name_prefix, make_input_name_prefix, sample_label, resolve_mem_gb, java_xmx
for n in ['Saccharomyces eubayanus', 'Apis mellifera ligustica', 'Brassica rapa subsp. pekinensis', 'Osmia lignaria_lignaria']:
    print(n, '->', make_name_prefix(n))
print(make_input_name_prefix(['Zacosmia_maculata_Abdomen', 'Zacosmia_maculata_Head', 'Zacosmia_maculata_Thorax']))
print(sample_label('Zacosmia_maculata_Head', ['Zacosmia_maculata']))
m = resolve_mem_gb('auto'); print('memory limit GB:', round(m), java_xmx(m))"
```

Expected output: `Saccharomyces_eubayanus`, `Apis_mellifera_ligustica`,
`Brassica_rapa_pekinensis`, `Osmia_lignaria_lignaria`, `Zacosmia_maculata`, `Head`,
then the job's memory limit (run inside a Slurm job to see the cgroup limit).

Option checks (these should fail with an error):

```shell
# conflicting naming options
apptainer exec --cleanenv nal-rnaseq_2026.11.13.sif \
  RNAseq_annotate.py -d . -s "Apis mellifera" -g x.fa --input_prefix -R
# per-sample mode with SRA input
apptainer exec --cleanenv nal-rnaseq_2026.11.13.sif \
  RNAseq_annotate.py -i runs.tsv -g x.fa --per_sample
```
