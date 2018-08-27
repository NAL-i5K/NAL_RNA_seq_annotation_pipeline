#!/bin/bash

#SBATCH --job-name="example RNAseq pipeline"

#SBATCH –p mem

#SBATCH -N 1

#SBATCH -n 4

 

#SBATCH -t 48:00:00

#SBATCH --mail-user=hsiaoyi0504@gmail.com

#SBATCH --mail-type=BEGIN,END,FAIL

#SBATCH -o "stdout"

#SBATCH -e "stderr"

#SBATCH --mem=200GB

export JAVA_TOOL_OPTIONS="-Xmx2g"

module load java/1.8.0_121
module load python_3/gcc/64/3.5.0
module load samtools/1.9
module load sratoolkit/gcc/64/2.8.2-1

cd /project/nal_genomics/leo/NAL_RNA_seq_annotation_pipeline 


ulimit -s unlimited # https://3.basecamp.com/3625179/buckets/5538276/messages/1243658335
RNAseq_annotate.py -i ./example/1049336.tsv -g ../Edan07162013.scaffolds.fa.gz -d
