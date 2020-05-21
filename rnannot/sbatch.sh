#!/bin/sh
#SBATCH -p medium
#SBATCH -N 1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-type=fail         # send email if job fails
#SBATCH --mail-user=k2025242322@gmail.com
echo $(date)
echo "load moudle java, sratoolkit, samtools and rsem"
module load java sratoolkit samtools rsem
echo "unload moudle python3 and perl"
module unload perl # We need to use the Perl that installed in our env
echo "module list"
module list
echo "start processing..."
#Example 1
python3 download_sra_metadata.py -t 1049336 -o 1049336.tsv
python3 RNAseq_annotate.py -i 1049336.tsv -g Edan07162013.scaffolds.fa.gz -a Edan_2.0 -t

#Example 2
python3 download_sra_metadata.py -t 7460 -o 7460.tsv
python3 RNAseq_annotate.py -i 7460.tsv -g GCF_000002195.4_Amel_4.5_genomic.fna.gz -a Amel_4.5 -t

#Example 3
python3 download_sra_metadata.py -t 7070 -o 7070.tsv
python3 RNAseq_annotate.py -i 7070.tsv -g GCF_000002335.3_Tcas5.2_genomic_RefSeqIDs.fna.gz -a Tcas5.2 -t
echo $(date)