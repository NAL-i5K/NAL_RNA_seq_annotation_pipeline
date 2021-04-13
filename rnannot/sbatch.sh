#!/bin/sh
#SBATCH -p medium
#SBATCH -N 1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-type=fail         # send email if job fails
#SBATCH --mail-user=<YOUREMAIL@email.com>
echo $(date)
echo "load module java, sratoolkit, samtools, rsem, kentutils"
module load java sratoolkit samtools rsem kentutils
echo "unload module perl"
module unload perl # We need to use the Perl that installed in our env
echo "module list"
module list
echo "start processing..."
#Example 1
download_sra_metadata.py -t 1049336 -o 1049336.tsv
RNAseq_annotate.py -i 1049336.tsv -g Edan07162013.scaffolds.fa.gz -a Edan_2.0 -t

#Example 2
download_sra_metadata.py -t 7460 -o 7460.tsv
RNAseq_annotate.py -i 7460.tsv -g GCF_000002195.4_Amel_4.5_genomic.fna.gz -a Amel_4.5 -t

#Example 3 - process 2 SRR files at most and remove intermediate files
download_sra_metadata.py -t 7070 -o 7070.tsv
RNAseq_annotate.py -i 7070.tsv -g GCF_000002335.3_Tcas5.2_genomic_RefSeqIDs.fna.gz -a Tcas5.2 -m 2

#Example 4
download_sra_metadata.py -t 486640 -o 486640.tsv
RNAseq_annotate.py -i 486640.tsv -g GCF_010583005.1_Obru_v1_genomic.fna.gz -a Obru_v1 -t

#Example 5 
download_sra_metadata.py -t 92692 -o 92692.tsv
RNAseq_annotate.py -i 92692.tsv -g _____RdoDt3_Drdd8_decomES.fasta -a RdoDt3 -t
echo $(date)
