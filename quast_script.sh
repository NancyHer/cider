#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=0-24:15:00     # day-hr:min:sec
#SBATCH --output=abc.out
#SBATCH --mail-user=nher002@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="test_Run"
#SBATCH -p batch # This is the default partition

# Print current date
date

#Quasting 

#module load java/7u40
module load QUAST

#curl https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/007/245/GCF_000007245.1_ASM724v1/GCF_000007245.1_ASM724v1_genomic.fna.gz --output xf_temecula1_zipped.fna.gz

#curl https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/019/765/GCF_000019765.1_ASM1976v1/GCF_000019765.1_ASM1976v1_genomic.fna.gz --output xf_m23_zipped.fna.gz

#curl https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/013/388/375/GCF_013388375.1_ASM1338837v1/GCF_013388375.1_ASM1338837v1_genomic.fna.gz --output xan_zipped.gz

#quast.py -o xf_xan_output_dir -t 12 -R xf_temecula1_zipped.fna.gz -l Xylella fastidiosa Temecula1, Xylella fastidiosa M23, Xanthomonas campestris pv. raphani xf_temecula1_zipped.fna.gz xf_m23_zipped.fna.gz xan_zipped.gz

quast.py -o xf_xan_output_dir -l Temecula1_Xf,M23_Xf,Xanthomonas_campestris_pv._raphani GCF_000007245.1_ASM724v1_genomic.fasta GCF_000019765.1_ASM1976v1_genomic.fasta GCF_013388375.1_ASM1338837v1_genomic.fasta



