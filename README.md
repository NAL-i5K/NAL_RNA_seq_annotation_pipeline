# NAL RNA-Seq Annotation Pipeline (Under Development)

[![Build Status](https://travis-ci.org/NAL-i5K/NAL_RNA_seq_annotation_pipeline.svg?branch=master)](https://travis-ci.org/NAL-i5K/NAL_RNA_seq_annotation_pipeline)

A RNA-Seq annotation pipeline based on [SRA Toolkit](https://github.com/ncbi/sra-tools), [fastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic), [HISAT2](https://github.com/infphilo/hisat2), [BBMap](https://sourceforge.net/projects/bbmap/), [picard](https://broadinstitute.github.io/picard/), [GATK3](https://github.com/broadgsa/gatk-protected), and [samtools](https://github.com/samtools/samtools). It's distributed as a python package.

## Prerequisite

- At least Python 3.5
- Java
- [SRA Toolkit](https://github.com/ncbi/sra-tools)
- [samtools](https://github.com/samtools/samtools)

## Installation

- `python setup.py install`. It will install a copy of FastQC, Trimmomatic, HISAT2, GATK3, and picard in this python package. You may need to add `--user` in arguments.

## Uninstallation

- `pip uninstall rnannot`

## Usage

``` shell
RNAseq_annotate.py [-h] [-i INPUT] [-g GENOME] [-n [NAME]]
                          [-o [OUTDIR]] [-d]

Easy to use pipeline built for large-scale RNA-seq mapping with a genome
assembly

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        A tsv file with a list of SRA runs' information.
  -g GENOME, --genome GENOME
                        A fasta file to align with.
  -n [NAME], --name [NAME]
                        name of the output folder, if not specified, use the
                        time of start
  -o [OUTDIR], --outdir [OUTDIR]
                        directory of output folder at, if not specified, use
                        current folder
  -d, --downsample      if specified, a downsampled bam file will be
                        downsampled
```

## Example

- `wget "https://i5k.nal.usda.gov/data/Arthropoda/ephdan-(Ephemera_danica)/Current%20Genome%20Assembly/1.Genome%20Assembly/BCM-After-Atlas/Scaffolds/Edan07162013.scaffolds.fa.gz"`
- `RNAseq_annotate.py -i ./example/1049336.tsv -g ./Edan07162013.scaffolds.fa.gz -d`

## Notes

- The input tsv should have at least five columns, including `Run`, `Platform`, `Model`, `LibraryLayout` (header must be presented), and `download_path`.
  - `Run` column represents the paths to the SRA files. You can use either relative path to your current directory or absolute path. To make less confusion, we recommned to use absolute path.
  - `Platform` column represents the sequencer's brand. We will recognize `ILLUMINA` and `ABI_SOLID` (although we will not process `ABI_SOLID`) in this field, because it determines the adapters used in the pipeline with `Model`.
  - `Model` column represents the sequencer's model. We will recongize `Illumina HiSeq ...`, `Illumina MiSeq ...`, and `Illumina Genome Analyzer II ...`, because it determines the adapters used in the pipeline with `Platform`. Don't forget that if you have any space in the name of the model. You need to escape them using quoting such as `"Illumina Hiseq 2000"` or back slash. 
  - `LibraryLayout` column represents what's the strategy of RNA-Seq experiment. It can be only `SINGLE` or `PAIRED`.
    - Currently, the `ABI_SOLID` sequencer is not supported.
    - For paired-end layout, `Trimmomatic` will produces four fastq files: forward\_paired, forward\_unpaired, reverse\_paired, reverse\_unpaired, but we will only use the paired data in alignment (by HISAT2)
  - `download_path` column represents where we can download the SRA files.
- When using on server, make sure you use the `JAVA_TOOL_OPTIONS` environment to set the maximum memory usage like `export JAVA_TOOL_OPTIONS="-Xmx2g"` when running the toolkit. You can also check an example [here](example/example_script.sh).

## Tests

### Test environment

- `python -m unittest -f tests/test_setup.py`

### Test parser

- `python -m unittest -f tests/test_parser.py`
