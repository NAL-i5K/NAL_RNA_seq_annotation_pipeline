# NAL RNA-Seq Annotation Pipeline

A RNA-Seq annotation pipeline based on SRA Toolkit, fastQC, Trimmomatic, HISAT2, BBMap, and samtools.

## Prerequisite

- Python 3.5
- Java
- SRA Toolkit
- samtools

## Installation

- `python setup.py`, this will install a copy of FastQC, Trimmomatic, and HISAT2 in the `rnannot/lib`.

## Example Usage

- `python3 ./rnannot/RNAseq_annotate.py -i ./example/197043.tsv -g ../GCA_000696855.1_Hvit_1.0_genomic.fna.gz`
- `python3 ./rnannot/RNAseq_annotate.py -i ./example/104688.tsv -g ../104688_ref_gapfilled_joined_lt9474.gt500.covgt10_chrMT_and_UN_refseq_IDs.fa.gz`

## Usage

## Notes

- The input tsv should have at least four columns, including `Run`, `Platform`, `Model`, and `LibraryLayout` (header must be presented).
  - `Run` column represents the paths to the SRA files. You can use either relative path to your current directory or absolute path. To make less confusion, we recommned to use absolute path.
  - `Platform` column represents the sequencer's brand. We will recognize `ILLUMINA` and `ABI_SOLID` (although we will not process `ABI_SOLID`) in this field, because it determines the adapters used in the pipeline with `Model`.
  - `Model` column represents the sequencer's model. We will recongize `Illumina HiSeq ...`, `Illumina MiSeq ...`, and `Illumina Genome Analyzer II ...`, because it determines the adapters used in the pipeline with `Platform`. Don't forget that if you have any space in the name of the model. You need to escape them using quoting such as `"Illumina Hiseq 2000"` or back slash. 
  - `LibraryLayout` column represents what's the strategy of RNA-Seq experiment. It can be only `SINGLE` or `PAIRED`.
  - Currently, the `ABI_SOLID` sequencer is not supported.
  - For paired-end layout, `Trimmomatic` will produces four fastq files: forward\_paired, forward\_unpaired, reverse\_paired, reverse\_unpaired, but we will only use the paired data in alignment (by HISAT2)


## Test Environment
- `python -m unittest -f test/test_setup.py`