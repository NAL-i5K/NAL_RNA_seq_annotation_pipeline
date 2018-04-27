# NAL RNA-Seq Annotation Pipeline

A RNA-Seq annotation pipeline based on SRA Toolkit, fastQC, Trimmomatic, HISAT2, and samtools.

## Prerequisite

- Python 3.5
- Java
- SRA Toolkit
- samtools

## Installation

- `python setup.py`, this will install a copy of FastQC, Trimmomatic, and HISAT2 in the `rnannot`.

## Example Usage

- `python3 ./rnannot/RNAseq_annotate.py -i ./example/197043.tsv -g ../GCA_000696855.1_Hvit_1.0_genomic.fna.gz`

## Test Environment
- `python -m unittest -f test/test_setup.py`
