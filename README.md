# NAL RNA-Seq Annotation Pipeline

A RNA-Seq annotation pipeline based on FastQC, Trimmomatic, and HISAT2.

## Prerequisite

- Python 3.5
- SRA Toolkit
- Java

## Installation

- `python setup.py`, this will install a copy of FastQC, Trimmomatic, and HISAT2 in the `rnannot`.

## Example Usage

- `python3 ./rnannot/RNAseq_annotate.py -m "Illumina Genome Analyzer II" -p ILLUMINA -l single -f /project/nal_genomics/leo/issue_270/output/sra_runs_data/SRR2095935 -g ../GCA_000696855.1_Hvit_1.0_genomic.fna.gz`

## Test Environment
- `python -m unittest -f test/test_setup.py`
