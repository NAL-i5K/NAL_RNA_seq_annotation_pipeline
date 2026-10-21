# NAL RNA-Seq Annotation Pipeline

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
image for HPC use.

## Overview

The pipeline has two parts.

**1. rnannot (alignment), run on HPC (e.g. SCINet Ceres/Atlas)**

1. Get reads, either by downloading SRA runs (`fasterq-dump`) or from a directory of local FASTQ files
2. FastQC
3. Adapter and quality trimming (Trimmomatic)
4. Read normalization (BBNorm, `target=100`)
5. Build the HISAT2 index (once per run, multithreaded)
6. Align with HISAT2 (`--dta-cufflinks`, which writes XS strand tags)
7. Sort and merge single- and paired-end alignments
8. CSI index, bigWig coverage (`bam_to_bigwig.py`), and junction BED (`regtools junctions extract -s XS`)

Outputs, in `<outdir>/<name>/`:

| File | Description |
|---|---|
| `<gggsss>_<assembly>_RNA-Seq-alignments_<date>.bam` | Coordinate-sorted alignments |
| `<gggsss>_<assembly>_RNA-Seq-alignments_<date>.bam.csi` | CSI index |
| `<gggsss>_<assembly>_RNA-Seq-alignments_<date>.bigwig` | Coverage track |
| `<gggsss>_<assembly>_RNA-Seq-alignments_<date>.bed` | Splice junctions |
| `Source.txt` | Scientific name, assembly name, date, and the list of runs/samples |

With `-t/--tempFile`, intermediate SAM/BAM files, merged normalized FASTQs, and the
HISAT2 index are kept.

**2. JBrowse loading, run on the Apollo server**

- `add_trackList.py` copies the outputs into the species directory and adds the BAM,
  bigWig, and junction tracks to `trackList.json`.
- `move_data.py` (legacy) syncs the loaded data to apollo-node1.

## Output naming (gggsss)

The organism code is derived from the scientific name and is always lowercase:

| Scientific name | Code | Rule |
|---|---|---|
| `Apis mellifera` | `apimel` | first 3 letters of genus + species |
| `Apis mellifera ligustica` | `apmeli` | first 2 letters each of genus, species, subspecies |
| `Brassica rapa subsp. pekinensis` | `brrape` | rank words (`subsp.`, `var.`, ...) are skipped |

To set the code yourself, use `--gggsss` in `RNAseq_annotate.py` and **the same value**
with `-g` in `add_trackList.py` and `move_data.py`. Do this for names that end in a
strain or isolate label (e.g. `Saccharomyces cerevisiae S288C`), which would otherwise
be treated as a subspecies.

`-R/--Run_prefix` names the outputs after the first run ID instead.

The downstream scripts depend on this naming format. If you change it,
`add_trackList.py` and `move_data.py` must be updated to match.

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
| wigToBigWig | UCSC binary |
| bam_to_bigwig | pinned commit (`BAM_TO_BIGWIG_COMMIT`) |

GATK is no longer included because the pipeline never used it.

## Installation

### Build the Docker image

```shell
git clone https://github.com/NAL-i5K/NAL_RNA_seq_annotation_pipeline.git
cd NAL_RNA_seq_annotation_pipeline
docker build -t nal-rnaseq:updated .
```

### Convert to Apptainer

On the machine where the image was built:

```shell
apptainer build nal-rnaseq.sif docker-daemon:nal-rnaseq:updated
```

The tag must be included. Building on one machine and running on another:

```shell
docker save nal-rnaseq:updated -o nal-rnaseq.tar
# copy nal-rnaseq.tar to the cluster, then:
apptainer build nal-rnaseq.sif docker-archive:nal-rnaseq.tar
```

After any code change, rebuild the Docker image **and** reconvert the SIF.

### Python package only (advanced)

```shell
pip install .
```

All external tools must then be installed separately. Set `RNANNOT_LIB` to the
directory containing them (the same layout as `rnannot/lib` in the container), and
make sure `samtools`, `fasterq-dump`, `rsem-bam2wig`, and `wigToBigWig` are on `PATH`.

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
RNAseq_annotate.py (-i INPUT | -d FASTQ_DIR -s "Genus species") -g GENOME [options]
```

| Option | Description |
|---|---|
| `-i, --input` | SRA metadata tsv (from `download_sra_metadata*.py`) |
| `-d, --fastq_dir` | Directory of local FASTQ files. Use instead of `-i`. Requires `-s` |
| `-s, --species` | Scientific name, e.g. `"Apis mellifera"` or `"Apis mellifera ligustica"`. Required with `-d` |
| `--platform` | Platform for local files. Default `ILLUMINA` |
| `--model` | Instrument model for local files, e.g. `"Illumina HiSeq 2500"` or `"Element AVITI"`. Selects the adapter set (see below) |
| `--gggsss` | Organism code for output names. Default: derived from the scientific name |
| `-g, --genome` | Genome FASTA (`.gz` accepted). Required |
| `-n, --name` | Output folder name (must not already exist). Default: start time |
| `-o, --outdir` | Parent directory for the output folder (must already exist). Default: current folder |
| `-a, --assembly` | Assembly name, used in output file names |
| `-t, --tempFile` | Keep intermediate files |
| `-m, --MaximumSRA` | SRA mode only: randomly select at most this many runs. Default 10 |
| `-T, --threads_num` | Number of threads, or `auto` (default) to use all CPUs allocated to the job. Values above the allocation are capped |
| `-R, --Run_prefix` | Name outputs after the first run ID instead of the gggsss prefix |

Exactly one of `-i` or `-d` must be given.

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
pull), the pipeline prints a warning and uses the first name.

### Input: local FASTQ directory

All `.fastq`, `.fq`, `.fastq.gz`, and `.fq.gz` files in the directory are processed.
Pairs are detected from the file names:

| File names | Detected as |
|---|---|
| `sampleA_S1_L001_R1_001.fastq.gz` + `..._R2_001.fastq.gz` | PAIRED, sample `sampleA_S1_L001` |
| `sampleA_R1.fq.gz` + `sampleA_R2.fq.gz` | PAIRED, sample `sampleA` |
| `sampleA_1.fastq` + `sampleA_2.fastq` | PAIRED, sample `sampleA` |
| `sampleB.fastq.gz` | SINGLE, sample `sampleB` |

- The detected samples are printed at startup. The run stops before any processing if
  there are duplicate read files or a sample has both single and paired files.
- Plain FASTQs are symlinked (originals are never modified or deleted). Gzipped FASTQs
  are decompressed into the run folder and removed at cleanup.
- Multi-lane samples are processed as separate samples and merged before alignment.
- `-m/--MaximumSRA` does not apply. Every file is used.

Preview the pairing without running the pipeline:

```shell
apptainer exec --cleanenv --bind "$FASTQ" nal-rnaseq.sif python3 -c \
  "from rnannot.RNAseq_annotate import discover_local_fastq; discover_local_fastq('$FASTQ')"
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

```shell
#!/bin/bash
#SBATCH --job-name=rnannot
#SBATCH --cpus-per-task=16        # at least 12: bam_to_bigwig.py uses 12 threads internally
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=rnannot_%j.out
#SBATCH --error=rnannot_%j.err

module load apptainer

SIF=$(readlink -f /project/nal_genomics/$USER/tools/nal-rnaseq.sif)
WORK=$(readlink -f /project/nal_genomics/$USER/rnaseq/apis)

# SRA mode
apptainer exec --cleanenv --bind "$WORK" --pwd "$WORK" "$SIF" \
    RNAseq_annotate.py -i "$WORK/7460.tsv" -g "$WORK/genome.fna.gz" \
        -a Amel_HAv3.1 -o "$WORK"

# Local mode (bind the FASTQ directory too)
# FASTQ=$(readlink -f /project/nal_genomics/$USER/raw_reads/apis)
# apptainer exec --cleanenv --bind "$WORK" --bind "$FASTQ" --pwd "$WORK" "$SIF" \
#     RNAseq_annotate.py -d "$FASTQ" -s "Apis mellifera ligustica" \
#         --model "Illumina HiSeq 2500" -g "$WORK/genome.fna.gz" -a Amel_HAv3.1 -o "$WORK"
```

Notes:
- `--cleanenv` removes Slurm environment variables inside the container. With
  `-T auto`, the thread count is still detected correctly from the job's CPU
  allocation.
- Failures in the final steps (sort check, indexing, bigWig, regtools) stop the job
  with a traceback in the `.err` file. Check that file as well as `.out`.
- Check memory use after a run with
  `sacct -j <jobid> --format=JobID,Elapsed,MaxRSS,ReqMem,State`.
- A plain `hisat2-build` needs roughly 8 GB of RAM for a human-sized genome.

## Loading into JBrowse/Apollo

### add_trackList.py

```shell
add_trackList.py -op OUTPUT_PATH -j JBROWSE_SIF -bam BAM -bigwig BIGWIG \
                 -csi CSI -bed BED -track TRACKLIST -s SOURCE [-g GGGSSS] [--overwrite]
```

| Option | Description |
|---|---|
| `-op, --output_path` | Species data root, e.g. `/app/data/other_species` |
| `-j, --jbrowse_path` | JBrowse container (run with `singularity exec`) |
| `-bam`, `-bigwig`, `-csi`, `-bed` | Pipeline output files |
| `-track, --input_track` | `trackList.json` path (a `.bak` backup is written first) |
| `-s, --Source` | `Source.txt` from the pipeline run |
| `-g, --gggsss` | Organism code. Use the same value given to the pipeline's `--gggsss` |
| `--overwrite` | Replace files that already exist in the analyses folder |

Files are copied to
`<output_path>/<gggsss>/<assembly>/scaffold/analyses/<gggsss>-RNA-Seq_v1.0/`.
This folder is shared by all RNA-Seq datasets for an assembly. The dated file names keep
datasets distinct. Track labels are the BAM name without `.bam`, plus `_coverage` and
`_junctions` for the bigWig and BED tracks.

The script warns if the BAM name does not start with the expected gggsss code.

### move_data.py (legacy)

```shell
move_data.py -Node1a user@apollo-node1.nal.usda.gov -s Source.txt \
             [-g GGGSSS] [--label LABEL] [--base_path /app/data/other_species]
```

Syncs the analyses folder and junction track to apollo-node1, creates the
`jbrowse/data/analyses` symlink if it is missing, and copies `trackList.json` last.
Use `--label` only if the pipeline ran with `-R/--Run_prefix`.

## Example

```shell
# 1. Metadata
download_sra_metadata.py -t 1049336 -o 1049336.tsv

# 2. Genome
wget "https://i5k.nal.usda.gov/data/Arthropoda/ephdan-(Ephemera_danica)/BCM-After-Atlas/1.Genome%20Assembly/BCM-After-Atlas/Scaffolds/Edan07162013.scaffolds.fa.gz"

# 3. Pipeline (inside the container, see the SLURM example above)
RNAseq_annotate.py -i 1049336.tsv -g Edan07162013.scaffolds.fa.gz -a Edan_2.0

# 4. Load into JBrowse
add_trackList.py -op /app/data/other_species -j /path/to/jbrowse.sif \
  -bam ephdan_Edan_2.0_RNA-Seq-alignments_2026-10-20.bam \
  -csi ephdan_Edan_2.0_RNA-Seq-alignments_2026-10-20.bam.csi \
  -bigwig ephdan_Edan_2.0_RNA-Seq-alignments_2026-10-20.bigwig \
  -bed ephdan_Edan_2.0_RNA-Seq-alignments_2026-10-20.bed \
  -track /app/data/other_species/ephdan/Edan_2.0/jbrowse/data/trackList.json \
  -s Source.txt
```

## Notes

- **SRA Toolkit configuration.** Current SRA Toolkit releases don't require interactive
  configuration. If downloads fail, run `vdb-config --restore-defaults` (outside
  `--cleanenv` runs this writes `~/.ncbi`).
- **Disk space.** `fasterq-dump` needs roughly 17× the accession size during
  conversion. The pipeline also writes uncompressed SAM files before sorting. Run in a
  project or scratch directory, not `$HOME`.
- **Normalization.** Reads are normalized with BBNorm before alignment. The bigWig shows
  normalized depth for annotation evidence, not expression levels.
- **Local samples in JBrowse metadata.** The "Data source" field links each entry in
  `Source.txt` to an SRA search. For local-mode runs these are sample names, so the
  links won't resolve.
- **Legacy Docker Hub image.** `k2025242322/i5k_rna_seq_annotation_pipeline` predates
  this update. Build the image from this repository instead.

## Tests

```shell
python -m unittest -f tests/test_parser.py
```

Quick checks after rebuilding the SIF:

```shell
apptainer exec --cleanenv nal-rnaseq.sif RNAseq_annotate.py --help
apptainer exec --cleanenv nal-rnaseq.sif python3 -c "
from rnannot.RNAseq_annotate import parse_scientific_name, make_gggsss
for n in ['Saccharomyces eubayanus', 'Apis mellifera ligustica']:
    print(n, '->', make_gggsss(parse_scientific_name(n)))"
```
