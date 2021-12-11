#!/usr/bin/env python3

import os,gzip
from Bio import SeqIO
from Bio.Seq import Seq

#Xanthamonas genome
x1="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/013/388/375/GCF_013388375.1_ASM1338837v1/GCF_013388375.1_ASM1338837v1_genomic.gff.gz"
#'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/013/388/375/GCF_013388375.1_ASM1338837v1/GCF_013388375.1_ASM1338837v1_genomic.fna.gz'
fasta=os.path.basename(x1)
if not os.path.exists(fasta):
    os.system("curl -O {}".format(x1))
    
aa_line_count=0
protein_coding=0
g_base=0
c_base=0
count = 0
with gzip.open(fasta,"rt") as infh:
    for line in infh:
        line_clean = line.strip().split("\t") #get rid of any white space ends of the line (.strip()) and .split() splits where ever you sepcify the character
        if line[0] != "#":
           #print (len(line_clean)) #ok so there is 9 items in each line
           #print(line_clean[8])
           aa_line_clean= line_clean[8].strip().split(";")
           aa_line_count += 1
           aa_line_count_accounted_for = aa_line_count/2 #there are repeats of the same info
           if aa_line_clean[4] == "gene_biotype=protein_coding":
               #print (aa_line_clean)
               protein_coding +=1
non_protein_coding_genes= aa_line_count_accounted_for-protein_coding
print("There is {} genes in the #Xanthomonas campestris pv. raphani genome.".format(aa_line_count_accounted_for))
print("{} are protein coding genes.".format(protein_coding))
print("{} are non-protein coding genes.".format(non_protein_coding_genes))
#print(aa_line_count)
#print(ninth_column)

proportion= (float(protein_coding)/float(aa_line_count_accounted_for))*100    
decimal_place="{:.2f}".format(proportion)
print("{}% of the genome is protein coding.".format(decimal_place))

##########3
#History of Junk Code

#!/usr/bin/env python3

#import os,gzip
#from Bio import SeqIO
#from Bio.Seq import Seq

#Xanthomonas campestris pv. raphani
#x1="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/013/388/375/GCF_013388375.1_ASM1338837v1/GCF_013388375.1_ASM1338837v1_genomic.gff.gz"
#'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/013/388/375/GCF_013388375.1_ASM1338837v1/GCF_013388375.1_ASM1338837v1_genomic.fna.gz'
#fasta=os.path.basename(x1)
#if not os.path.exists(fasta):
#    os.system("curl -O {}".format(x1))

#xan_cds_total=0
#xan_whole_length=4942039

#with gzip.open(fasta,"rt") as xcr: 
#    for seq_record in SeqIO.parse( xcr , "fasta"):
#        #print(len(seq_record))
#        xan_cds_total = (len(seq_record))+xan_cds_total # gives the total; length of the cds of the genome, need to define the variable first for it to work
#xan_intron_length= xan_whole_length-int(xan_cds_total)            
               
#print('The whole genome of Xanthomonas campestris pv. raphani is {}bp.'.format(xan_whole_length))
#print('The total length calculated with code is {}bp.'.format(xan_cds_total))
#print('The difference should be zero, {}bp. If not the assembly may have discrepancies.'.format(xan_intron_length))
       
