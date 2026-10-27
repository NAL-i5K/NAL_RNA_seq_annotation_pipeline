# NAL RNA-Seq Annotation Pipeline

**Current release: 2026.10.26**

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

## Overview

The pipeline has two parts.

**1. rnannot (alignment), run on HPC (e.g. SCINet Ceres/Atlas)**

1. Get reads: download SRA runs (`fasterq-dump`) or use a directory of local FASTQ files
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
| `<Genus_species>_<assembly>_RNA-Seq-alignments_<date>.bam` | Coordinate-sorted alignments |
| `<Genus_species>_<assembly>_RNA-Seq-alignments_<date>.bam.csi` | CSI index |
| `<Genus_species>_<assembly>_RNA-Seq-alignments_<date>.bigwig` | Coverage track |
| `<Genus_species>_<assembly>_RNA-Seq-alignments_<date>.bed` | Splice junctions |
| `Source.txt` | Scientific name, assembly name, date, and list of runs/samples |

With `-t/--tempFile`, intermediate SAM/BAM files, merged normalized FASTQs, and the
HISAT2 index are kept.

**2. JBrowse loading, run on the Apollo server**

- `add_trackList.py` copies the outputs into the species directory and adds the BAM,
  bigWig, and junction tracks to `trackList.json`.
- `move_data.py` (legacy) syncs loaded data to apollo-node1.
- Both use the shared naming rules in `i5k_naming.py`, which must be in the same
  directory as the two scripts.

## Naming

### Output files (pipeline)

| Scientific name | File prefix |
|---|---|
| `Apis mellifera` | `Apis_mellifera` |
| `Apis mellifera ligustica` | `Apis_mellifera_ligustica` |
| `Brassica rapa subsp. pekinensis` | `Brassica_rapa_pekinensis` (rank words such as `subsp.` and `var.` are dropped) |

- The name is used exactly as written (case is not changed). SRA metadata is already
  correctly capitalized. With `-d`, type the name properly, e.g. `-s "Apis mellifera"`.
- Names with a strain or cultivar label keep it as the last part
  (`Saccharomyces cerevisiae S288C` → `Saccharomyces_cerevisiae_S288C`). Use `--prefix`
  to choose something else.
- `-R/--Run_prefix` names the outputs after the first run ID instead.

### JBrowse species directories (loaders)

The i5k directory convention is changing. Legacy species keep their `gggsss`
directories, and new species use full names. The loaders detect which applies:

| Situation | Species directory |
|---|---|
| A legacy directory already exists (`apimel`, or `apmeli` for a trinomial) | The legacy directory |
| New species | `Genus_species[_subspecies]` |
| `--species_dir` given | Exactly that |
| Both legacy and full-name directories exist | Stops; choose with `--species_dir` |

`--species_dir` (aliases `-g`, `--gggsss`) accepts letters, digits, and underscores,
so both `apimel` and `Saccharomyces_eubayanus` are valid. If an underscore is rejected
with "must contain only letters and digits", you are running a pre-2026.10.26 copy of
the loader scripts. Replace `add_trackList.py`, `move_data.py`, and `i5k_naming.py`
with the current versions.

The analyses folder is `<species_dir>-RNA-Seq_v1.0`, shared by all RNA-Seq datasets for
an assembly. The dated file names keep datasets distinct.

**Always run `add_trackList.py --dry_run` the first time you load a species.** Short
legacy codes can collide (two species with the same `gggsss`). The script prints which
directory it chose and why.

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
docker build --pull --no-cache --progress=plain -t nal-rnaseq:2026.10.26 . 2>&1 | tee build-2026.10.26.log
```

The build needs outbound TCP port 873 (rsync) for the wigToBigWig download.

### Convert to Apptainer

On the machine where the image was built:

```shell
apptainer build nal-rnaseq_2026.10.26.sif docker-daemon:nal-rnaseq:2026.10.26
```

Building on one machine and running on another:

```shell
docker save nal-rnaseq:2026.10.26 -o nal-rnaseq_2026.10.26.tar
# copy the .tar to the cluster, then:
apptainer build nal-rnaseq_2026.10.26.sif docker-archive:nal-rnaseq_2026.10.26.tar
```

After any code change, rebuild the Docker image **and** reconvert the SIF.

### Python package only (advanced)

```shell
pip install .
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
RNAseq_annotate.py (-i INPUT | -d FASTQ_DIR -s "Genus species") -g GENOME [options]
```

| Option | Description |
|---|---|
| `-i, --input` | SRA metadata tsv (from `download_sra_metadata*.py`) |
| `-d, --fastq_dir` | Directory of local FASTQ files. Use instead of `-i`. Requires `-s` |
| `-s, --species` | Scientific name, e.g. `"Apis mellifera"` or `"Apis mellifera ligustica"`. Required with `-d` |
| `--platform` | Platform for local files. Default `ILLUMINA` |
| `--model` | Instrument model for local files, e.g. `"Illumina HiSeq 2500"` or `"Element AVITI"`. Selects the adapter set |
| `--prefix` | Output file prefix. Default: `Genus_species[_subspecies]` from the scientific name |
| `-g, --genome` | Genome FASTA (`.gz` accepted). Required |
| `-n, --name` | Output folder name (must not already exist). Default: start time |
| `-o, --outdir` | Parent directory for the output folder (must already exist). Default: current folder |
| `-a, --assembly` | Assembly name, used in output file names |
| `-t, --tempFile` | Keep intermediate files |
| `-m, --MaximumSRA` | SRA mode only: randomly select at most this many runs. Default 10 |
| `-T, --threads_num` | Number of threads, or `auto` (default) to use all CPUs allocated to the job. Values above the allocation are capped |
| `-M, --mem_gb` | Memory in GB, or `auto` (default) to use the job's Slurm/cgroup memory limit. BBTools Java heaps are capped at 80% of this. Values above the limit are capped |
| `-R, --Run_prefix` | Name outputs after the first run ID. Cannot be combined with `--prefix` |

Exactly one of `-i` or `-d` must be given. The output name is decided and printed
before any processing starts.

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
the name yourself.

### Input: local FASTQ directory

All `.fastq`, `.fq`, `.fastq.gz`, and `.fq.gz` files in the directory are processed.
Pairs are detected from file names:

| File names | Detected as |
|---|---|
| `sampleA_S1_L001_R1_001.fastq.gz` + `..._R2_001.fastq.gz` | PAIRED, sample `sampleA_S1_L001` |
| `sampleA_R1.fq.gz` + `sampleA_R2.fq.gz` | PAIRED, sample `sampleA` |
| `sampleA_1.fastq` + `sampleA_2.fastq` | PAIRED, sample `sampleA` |
| `sampleB.fastq.gz` | SINGLE, sample `sampleB` |

- Detected samples are printed at startup. The run stops before processing if there
  are duplicate read files or a sample has both single and paired files.
- Plain FASTQs are symlinked (originals are never modified or deleted). Gzipped FASTQs
  are decompressed into the run folder and removed at cleanup.
- Multi-lane samples are processed separately and merged before alignment.
- `-m/--MaximumSRA` does not apply. Every file is used.

Preview the pairing without running the pipeline:

```shell
apptainer exec --cleanenv --bind "$FASTQ" nal-rnaseq_2026.10.26.sif python3 -c \
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

SIF=$(readlink -f /project/nal_genomics/tools/nal-rnaseq/nal-rnaseq_2026.10.26.sif)
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

### Threads (`-T`) and memory (`-M`)

- By default the pipeline uses every CPU and all the memory Slurm allocated to the job.
  Both are detected from the process itself (CPU affinity and the cgroup memory
  limit), so detection works with `--cleanenv`.
- The job log shows `Using N threads` and `Using N GB memory limit (BBTools -Xmx...)`.
  These should match `--cpus-per-task` and `--mem`.
- BBTools (BBNorm, BBMerge, repair.sh) get an explicit `-Xmx` of 80% of the limit.
  Without it, BBTools size their heap from the whole node's RAM, which can exceed the
  job's allocation and cause an OOM kill.
- If the log warns that no cgroup memory limit was found (e.g. outside Slurm), pass
  `-M` with the memory you want the run to use.
- Always set `--cpus-per-task`. Without it, Slurm allocates one CPU and the pipeline
  runs single-threaded.
- If you wrap the command in `srun`, repeat the CPU request:
  `srun --cpus-per-task="$SLURM_CPUS_PER_TASK" apptainer exec ...`
- `samtools sort` memory is per thread, so raise `--mem` along with `--cpus-per-task`.

### Monitoring

- Failures in the final steps (sort check, indexing, bigWig, regtools) stop the job
  with a traceback in the `.err` file. Check that file as well as `.out`.
- Check resource use after a run:
  `sacct -j <jobid> --format=JobID,Elapsed,MaxRSS,ReqMem,State` or `seff <jobid>`.
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
| `-s, --Source` | `Source.txt` from the pipeline run |
| `--species_dir` (aliases `-g`, `--gggsss`) | Species directory (letters, digits, underscores). Default: existing legacy directory, otherwise full name |
| `--prefix` | File prefix used by the pipeline. Use the same value given to the pipeline's `--prefix` |
| `--overwrite` | Replace files that already exist in the analyses folder |
| `--dry_run` | Print the resolved paths and exit without copying or editing anything |

Files are copied to
`<output_path>/<species_dir>/<assembly>/scaffold/analyses/<species_dir>-RNA-Seq_v1.0/`.
Track labels are the BAM name without `.bam`, plus `_coverage` and `_junctions` for the
bigWig and BED tracks. The script warns if the BAM name doesn't start with the expected
prefix.

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

## Example

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

## Notes

- **SRA Toolkit configuration.** Current SRA Toolkit releases don't need interactive
  configuration. If downloads fail, run `vdb-config --restore-defaults`.
- **Disk space.** `fasterq-dump` needs roughly 17× the accession size during
  conversion, and the pipeline writes uncompressed SAM files before sorting. Run in a
  project or scratch directory, not `$HOME`.
- **Normalization.** Reads are normalized with BBNorm before alignment. The bigWig
  shows normalized depth for annotation evidence, not expression levels.
- **Local samples in JBrowse metadata.** The "Data source" field links each entry in
  `Source.txt` to an SRA search. For local-mode runs these are sample names, so the
  links won't resolve.
- **Legacy Docker Hub image.** `k2025242322/i5k_rna_seq_annotation_pipeline` predates
  this release. Use the image built from this repository.

## Tests

```shell
python -m unittest -f tests/test_parser.py
```

Quick checks after rebuilding the SIF:

```shell
apptainer exec --cleanenv nal-rnaseq_2026.10.26.sif RNAseq_annotate.py --help
apptainer exec --cleanenv nal-rnaseq_2026.10.26.sif python3 -c "
from rnannot.RNAseq_annotate import make_name_prefix, resolve_mem_gb, java_xmx
for n in ['Saccharomyces eubayanus', 'Apis mellifera ligustica', 'Brassica rapa subsp. pekinensis']:
    print(n, '->', make_name_prefix(n))
m = resolve_mem_gb('auto'); print('memory limit GB:', round(m), java_xmx(m))"
```

Expected output: `Saccharomyces_eubayanus`, `Apis_mellifera_ligustica`,
`Brassica_rapa_pekinensis`, then the job's memory limit (run the second check inside a
Slurm job to see the cgroup limit).
