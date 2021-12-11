#!/usr/bin/env python3

import os,gzip
from Bio import SeqIO
from Bio.Seq import Seq

#m23 Xf
m23='https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/019/765/GCA_000019765.1_ASM1976v1/GCA_000019765.1_ASM1976v1_genomic.gff.gz'
fasta=os.path.basename(m23)
if not os.path.exists(fasta):
    os.system("curl -O {}".format(m23))

aa_line_count=0
protein_coding=0
g_base=0
c_base=0
count = 0
ninth_column= []
with gzip.open(fasta,"rt") as infh:
    for line in infh:
        line_clean = line.strip().split("\t") #get rid of any white space ends of the line (.strip()) and .split() splits where ever you sepcify the character
        if line[0] != "#":
           #print (len(line_clean)) #ok so there is 9 items in each line
           #print(line_clean[8])
           #ninth_column.append(line)
           aa_line_clean= line_clean[8].strip().split(";")
           #print (aa_line_clean)
           aa_line_count += 1
           aa_line_count_accounted_for = aa_line_count/2 #there are repeats of the same info
           if aa_line_clean[3] == "gene_biotype=protein_coding":
               #print (aa_line_clean)
               protein_coding +=1
non_protein_coding_genes= aa_line_count_accounted_for-protein_coding
print("There is {} genes in the Xylella fastidiosa M23 genome.".format(aa_line_count_accounted_for))
print("{} are protein coding genes.".format(protein_coding))
print("{} are non-protein coding genes.".format(non_protein_coding_genes))
#print(aa_line_count)
#print(ninth_column)
           
proportion= (float(protein_coding)/float(aa_line_count_accounted_for))*100    
decimal_place="{:.2f}".format(proportion)
print("{}% of the genome is protein coding.".format(decimal_place))
    


       


#####################3
#History of Junk Code
#!/usr/bin/env python3

#import os,gzip
#from Bio import SeqIO
#from Bio.Seq import Seq

#Xylella fastidiosa M23
#m2='https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/019/765/GCA_000019765.1_ASM1976v1/GCA_000019765.1_ASM1976v1_genomic.gff.gz'
#fasta=os.path.basename(m2)
#if not os.path.exists(fasta):
#    os.system("curl -O {}".format(m2))
#a=0
vt1_total_length=2535690

#with gzip.open(fasta,"rt") as xfm:
#    for seq_record in SeqIO.parse( xfm , "fasta"):
#           a = (len(seq_record))+a # gives the total; length of the cds of the genome, need to define the variable first for it to work
#intron_length= t1_total_length-int(a)            
               
#print('The whole genome of Xylella fastidiosa M23 is {}bp.'.format(t1_total_length))
#print('The total length calculated with code is {}bp.'.format(a))
#print('The total length of intons is {}bp.'.format(intron_length))
