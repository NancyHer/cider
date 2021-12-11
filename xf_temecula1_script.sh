#!/usr/bin/env python3

import os,gzip
from Bio import SeqIO
from Bio.Seq import Seq

#temecula1 Xf genome
t1_cds="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/007/245/GCF_000007245.1_ASM724v1/GCF_000007245.1_ASM724v1_genomic.gff.gz"
#"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/007/245/GCF_000007245.1_ASM724v1/GCF_000007245.1_ASM724v1_genomic.fna.gz"
#"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/007/245/GCF_000007245.1_ASM724v1/GCF_000007245.1_ASM724v1_genomic.gff.gz"
#"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/007/245/GCF_000007245.1_ASM724v1/GCF_000007245.1_ASM724v1_genomic.fna.gz"
#'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/007/245/GCF_000007245.1_ASM724v1/GCF_000007245.1_ASM724v1_cds_from_genomic.fna.gz' #this file is only the cds
fasta=os.path.basename(t1_cds)
if not os.path.exists(fasta):
    os.system("curl -O {}".format(t1_cds))

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
           aa_line_count += 1
           aa_line_count_accounted_for = aa_line_count/2 #there are repeats of the same info
           if aa_line_clean[4] == "gene_biotype=protein_coding":
               #print (aa_line_clean)
               protein_coding +=1
non_protein_coding_genes= aa_line_count_accounted_for-protein_coding
print("There is {} genes in the Xylella fastidiosa Temecula1 genome.".format(aa_line_count_accounted_for))
print("{} are protein coding genes.".format(protein_coding))
print("{} are non-protein coding genes.".format(non_protein_coding_genes))
#print(aa_line_count)
#print(ninth_column)
           
proportion= (float(protein_coding)/float(aa_line_count_accounted_for))*100    
decimal_place="{:.2f}".format(proportion)
print("{}% of the genome is protein coding.".format(decimal_place))

###########
#history of junk code
#!/usr/bin/env python3

#import os,gzip
#from Bio import SeqIO
#from Bio.Seq import Seq

#t1_cds=
#"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/007/245/GCF_000007245.1_ASM724v1/GCF_000007245.1_ASM724v1_genomic.gff.gz"
#"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/007/245/GCF_000007245.1_ASM724v1/GCF_000007245.1_ASM724v1_genomic.fna.gz"
#'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/007/245/GCF_000007245.1_ASM724v1/GCF_000007245.1_ASM724v1_cds_from_genomic.fna.gz' #this file is only the cds

#gzip -d GCF_000007245.1_ASM724v1_genomic.gff.gz # for u

#fasta=os.path.basename(t1_cds)
#if not os.path.exists(fasta):
#    os.system("curl -O {}".format(t1_cds))

#a=0
#t1_total_length=2519802

#aa=0
#ab=0
#abcd=0
#total_length=0
#with gzip.open(fasta,"rt") as infh:
   
   #for first in range (20):
        #lines_20= infh.readline()
        #print(lines_20) #to see what the first 20 lines look like
    
#    for line in infh:
        #print(line)
       
#        if line[0] != "#":
            #if line [1] == "RefSeq":
#            print ((line{[9])) #checkpoint
        #line = line.strip('') #This should get rid of any white space from both ends of the line
        #linearray = line.strip().split("[")
        #ab += 1 #so it seems like there is 586166 lines here, where as the filtering below has 232
       
       #if linearray[0] == ">":
            #print (linearray) #checkpoint
        #print (len(linearray))
       
       #if len(linearray) > 2:
            #aa=linearray
            #b=['>']
            #c=any(b in aa for b in aa) #this is some code trickery, anyways it is something to get your head around each time
            #length = linearray[3]
            #length_values= int(length.replace("length=","")) #remembr that int() makes numbers in strings into intergers, xx.replace means to literally replace it with the second value
            #total_length= total_length + length_values
            #print (length) # this prints out a list of the fourth item in the linearry list, which all have length=xxx, time to strip 'length=
            #print(length_values)
            #abcd += 1 # there is 232 lines in here that probably have '>'
            #print(a) #this checks if all of the lines in linearray starts with >, which it does
   
   #print (c) #This prints if statement is true or false
    #print (ab) #This prints how many lines in the whole file
    #print (abcd) #This prints how many lines have >        
#coding_percentage= (total_protein_length/total_length)*100
#decimal_place="{:.2f}".format(coding_percentage)
   
   
   #for seq_record in SeqIO.parse( infh , "fasta"):
            #a = (len(seq_record))+a # gives the tota; length of the cds of the genome, need to define the variable first for it to work
#length_difference= t1_total_length-int(a)


               
#print('The whole genome of Xylella fastidiosa Temecula1 is {}bp.'.format(t1_total_length))
#print('The total length calculated with code is {}bp.'.format(a))
#print('The difference should be zero, {}bp. If not the assembly may have discrepancies.'.format(length_difference))

#prot_file.close()



#!/usr/bin/env python3

#import os,gzip
#from Bio import SeqIO
#from Bio.Seq import Seq

#t1_cds="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/007/245/GCF_000007245.1_ASM724v1/GCF_000007245.1_ASM724v1_genomic.fna.gz"
#"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/007/245/GCF_000007245.1_ASM724v1/GCF_000007245.1_ASM724v1_genomic.gff.gz"
#"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/007/245/GCF_000007245.1_ASM724v1/GCF_000007245.1_ASM724v1_genomic.fna.gz"
#'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/007/245/GCF_000007245.1_ASM724v1/GCF_000007245.1_ASM724v1_cds_from_genomic.fna.gz' #this file is only the cds
#fasta=os.path.basename(t1_cds)
#if not os.path.exists(fasta):
#    os.system("curl -O {}".format(t1_cds))

#a_base=0
#t_base=0
#g_base=0
#c_base=0
#with gzip.open(fasta,"rt") as infh:
#    for line in infh:
#        if line[0] != ">":
#            def countA(text):
#        count = 0
#        for c in text:
#            if c == 'a':
#                count = count + 1
#        return count

#    print(countA(line))
            #print(str_line) #so it does print the fasta nucleotide sequence, so yes it works
#            for a_n in line:
 #               if a_n == "A":
 #                   a_base = a_base+1
  #          for t_n in line:
   #             if t_n == "T":
    #                t_base = a_base+1
     #       for g_n in line:
      #          if g_n == "G":
       #             g_base = a_base+1
        #    for c_n in line:
         #       if c_n == "C":
#                    c_base = a_base+1
 #   print(a_base,t_base,g_base,c_base)
  #          
#!/usr/bin/env python3

#import os,gzip
#from Bio import SeqIO
#from Bio.Seq import Seq

#t1_cds="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/007/245/GCF_000007245.1_ASM724v1/GCF_000007245.1_ASM724v1_genomic.fna.gz"
#"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/007/245/GCF_000007245.1_ASM724v1/GCF_000007245.1_ASM724v1_genomic.gff.gz"
#"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/007/245/GCF_000007245.1_ASM724v1/GCF_000007245.1_ASM724v1_genomic.fna.gz"
#'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/007/245/GCF_000007245.1_ASM724v1/GCF_000007245.1_ASM724v1_cds_from_genomic.fna.gz' #this file is only the cds
#fasta=os.path.basename(t1_cds)
#if not os.path.exists(fasta):
#    os.system("curl -O {}".format(t1_cds))

#a_base=0
#t_base=0
#g_base=0
#c_base=0
#count = 0
#with gzip.open(fasta,"rt") as infh:
#    for line in infh:
#        if line[0] != ">":
#            clean_line = line
#def countA(text):
#    for c in text:
#        if c == 'a':
#            count +=1
#    return count

#print(countA(clean_line))
