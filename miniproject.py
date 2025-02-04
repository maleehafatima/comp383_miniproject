import os
import shutil
#1 get files
os.system('mkdir miniproject_Maleeha_Fatima') #creates directory for project
os.chdir('miniproject_Maleeha_Fatima') #sets that new directory as woring directory

typeRun = input("What data would you like to run? Enter either 'full' or 'test'\n") #asks user what type of run they want to do

if typeRun == "full":  #uses the full data
    os.system('wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660030/SRR5660030.1') #downloads SRR files
    os.system('wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660033/SRR5660033.1')
    os.system('wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660044/SRR5660044.1')
    os.system('wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660045/SRR5660045.1')

    os.system('fastq-dump -I --split-files SRR5660030.1') #splits into forward and reverse files
    os.system('fastq-dump -I --split-files SRR5660033.1')
    os.system('fastq-dump -I --split-files SRR5660044.1')
    os.system('fastq-dump -I --split-files SRR5660045.1')

elif typeRun == "test": #uses the shortened test data files provided and renames and moves them into project directory
    os.chdir("..") 
    shutil.copyfile('testdata/SRR5660030.1_1test.fastq', 'miniproject_Maleeha_Fatima/SRR5660030.1_1.fastq')
    shutil.copyfile('testdata/SRR5660030.1_2test.fastq', 'miniproject_Maleeha_Fatima/SRR5660030.1_2.fastq')
    shutil.copyfile('testdata/SRR5660033.1_1test.fastq', 'miniproject_Maleeha_Fatima/SRR5660033.1_1.fastq')
    shutil.copyfile('testdata/SRR5660033.1_2test.fastq', 'miniproject_Maleeha_Fatima/SRR5660033.1_2.fastq')
    shutil.copyfile('testdata/SRR5660044.1_1test.fastq', 'miniproject_Maleeha_Fatima/SRR5660044.1_1.fastq')
    shutil.copyfile('testdata/SRR5660044.1_2test.fastq', 'miniproject_Maleeha_Fatima/SRR5660044.1_2.fastq')
    shutil.copyfile('testdata/SRR5660045.1_1test.fastq', 'miniproject_Maleeha_Fatima/SRR5660045.1_1.fastq')
    shutil.copyfile('testdata/SRR5660045.1_2test.fastq', 'miniproject_Maleeha_Fatima/SRR5660045.1_2.fastq')
    os.chdir('miniproject_Maleeha_Fatima')
    
else:
    print("Input error: please enter either 'full' or 'test'")


#2 finds number of CDS and stores them
outfile = open('miniproject.log','w')
cdsoutput = open('HCMVcds.fasta','w')

from Bio import SeqIO 
from Bio import Entrez
Entrez.email = "mfatima2@luc.edu"

handle = Entrez.efetch(db = "nucleotide", id = "EF999921", rettype = "gbwithparts", retmode = "text") #uses Entrez to search the nucleotide database using the id number 
entry = SeqIO.read(handle, format = "gb") #reads the handle

cdscount = 0 #initialize counter for number of CDS

for feat in entry.features: #for each feature in the entry
    if feat.type == "CDS": #checks if CDs
        cdscount += 1 #adds one to counter
        featName = feat.qualifiers["protein_id"][0] #saves protein name for the specific CDS
        featSeq = feat.extract(entry.seq) #saves its sequence
        cdsoutput.write(">" + str(featName) + "\n" + str(featSeq) + "\n" ) #writes it out in FASTA format
cdsoutput.close()

outfile.write("The HCMV genome (EF999921) has " + str(cdscount) + " CDS.\n")

os.system('time kallisto index -i HCMVcds.idx HCMVcds.fasta') #builds kallisto index

#3 quantify TPM of each CDS with kallisto
os.system('mkdir results') #creates folder to store results of kallisto quant
os.system('time kallisto quant -i HCMVcds.idx -o results/SRR5660030.1 -b 30 -t 2 SRR5660030.1_1.fastq SRR5660030.1_2.fastq') 
os.system('time kallisto quant -i HCMVcds.idx -o results/SRR5660033.1 -b 30 -t 2 SRR5660033.1_1.fastq SRR5660033.1_2.fastq')
os.system('time kallisto quant -i HCMVcds.idx -o results/SRR5660044.1 -b 30 -t 2 SRR5660044.1_1.fastq SRR5660044.1_2.fastq')
os.system('time kallisto quant -i HCMVcds.idx -o results/SRR5660045.1 -b 30 -t 2 SRR5660045.1_1.fastq SRR5660045.1_2.fastq')

#create table with the 4 samples and their paths to the results
kalltable = open('kallistoOutputTable.txt', 'w')
kalltable.write('sample condition path\n')
kalltable.write('SRR5660030.1 2dpi results/SRR5660030.1\n')
kalltable.write('SRR5660033.1 6dpi results/SRR5660033.1\n')
kalltable.write('SRR5660044.1 2dpi results/SRR5660044.1\n')
kalltable.write('SRR5660045.1 6dpi results/SRR5660045.1\n')
kalltable.close()

#creates R file for performing sleuth
sleuth = open('sleuth.R', 'w')
sleuth.write('library(sleuth)\n')
sleuth.write('library(dplyr)\n')
sleuth.write('ktab <- read.table("kallistoOutputTable.txt",header=TRUE,stringsAsFactors=FALSE)\n')
sleuth.write('so <- sleuth_prep(ktab)\n')
sleuth.write("so <- sleuth_fit(so, ~condition, 'full') #fit a model comparing the two conditions\n")
sleuth.write("so <- sleuth_fit(so, ~1, 'reduced') #fit the reduced model to compare in the likelihood ratio test\n")
sleuth.write("so <- sleuth_lrt(so, 'reduced', 'full') #perform the likelihood ratio test for diffferential expression between conditions\n")
sleuth.write("sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE) #extract the test results from the sleuth object\n")
sleuth.write('sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05) %>% dplyr::arrange(pval)  #filter most sig results (FDR/qval < 0.05) and sort by pval\n')
sleuth.write("selected <- sleuth_significant[,c(1,4,2,3)] #picks the specific columns\n")
sleuth.write("write.table(selected, file = 'fdr05_results.txt', sep='\t', quote = FALSE, row.names = FALSE) #write FDR <0.05 transcripts to file\n")
sleuth.close()
os.system('Rscript sleuth.R')

fdrResults = open('fdr05_results.txt', 'r')
for line in fdrResults.readlines(): #writes the results of sleuth into log file
    outfile.write(line + '\n')

fdrResults.close()

#4 build bowtie index

# get the full HCMV genome sequence to build the bowtie index
completeHCMV = open('completeHCMV.fasta', 'w') #will write complete seq to this file

handle = Entrez.efetch(db = "nucleotide", id = "EF999921", rettype = "fasta", retmode = "text")
entry = SeqIO.read(handle, "fasta")

completeHCMV.write(">" + entry.description + "\n" + str(entry.seq))
completeHCMV.close()

os.system('bowtie2-build completeHCMV.fasta HCMVindex') #uses bowtie2 to build index from completeHCMV

#perform bowtie2 mapping
#--al-conc will output fastq file of bowtie filtered reads
os.system("bowtie2 -x HCMVindex -1 SRR5660030.1_1.fastq -2 SRR5660030.1_2.fastq -S HCMVmap.sam --al-conc SRR5660030.1_mapped.fq")
os.system("bowtie2 -x HCMVindex -1 SRR5660033.1_1.fastq -2 SRR5660033.1_2.fastq -S HCMVmap.sam --al-conc SRR5660033.1_mapped.fq")
os.system("bowtie2 -x HCMVindex -1 SRR5660044.1_1.fastq -2 SRR5660044.1_2.fastq -S HCMVmap.sam --al-conc SRR5660044.1_mapped.fq")
os.system("bowtie2 -x HCMVindex -1 SRR5660045.1_1.fastq -2 SRR5660045.1_2.fastq -S HCMVmap.sam --al-conc SRR5660045.1_mapped.fq")

#counts number of reads before and after bowtie filtering for each donor
beforeFile = open('SRR5660030.1_1.fastq', 'r')
beforeCount = 0
for line in beforeFile.readlines():
    beforeCount += 1
beforeCount = int(beforeCount / 4)

afterFile = open('SRR5660030.1_mapped.1.fq', 'r')
afterCount = 0
for line in afterFile.readlines():
    afterCount += 1
afterCount = int(afterCount / 4)

outfile.write('Donor 1 (2dpi) had ' + str(beforeCount) + ' read pairs before Bowtie2 filtering and ' + str(afterCount) + ' read pairs after.\n')

beforeFile = open('SRR5660033.1_1.fastq', 'r')
beforeCount = 0
for line in beforeFile.readlines():
    beforeCount += 1
beforeCount = int(beforeCount / 4)

afterFile = open('SRR5660033.1_mapped.1.fq', 'r')
afterCount = 0
for line in afterFile.readlines():
    afterCount += 1
afterCount = int(afterCount / 4)

outfile.write('Donor 1 (6dpi) had ' + str(beforeCount) + ' read pairs before Bowtie2 filtering and ' + str(afterCount) + ' read pairs after.\n')

beforeFile = open('SRR5660044.1_1.fastq', 'r')
beforeCount = 0
for line in beforeFile.readlines():
    beforeCount += 1
beforeCount = int(beforeCount / 4)

afterFile = open('SRR5660044.1_mapped.1.fq', 'r')
afterCount = 0
for line in afterFile.readlines():
    afterCount += 1
afterCount = int(afterCount / 4)

outfile.write('Donor 3 (2dpi) had ' + str(beforeCount) + ' read pairs before Bowtie2 filtering and ' + str(afterCount) + ' read pairs after.\n')

beforeFile = open('SRR5660045.1_1.fastq', 'r')
beforeCount = 0
for line in beforeFile.readlines():
    beforeCount += 1
beforeCount = int(beforeCount / 4)

afterFile = open('SRR5660045.1_mapped.1.fq', 'r')
afterCount = 0
for line in afterFile.readlines():
    afterCount += 1
afterCount = int(afterCount / 4)

outfile.write('Donor 3 (6dpi) had ' + str(beforeCount) + ' read pairs before Bowtie2 filtering and ' + str(afterCount) + ' read pairs after.\n')

#5 run spades and write command to log file

os.system("spades -k 127 -t 2 --only-assembler --pe1-1 SRR5660030.1_mapped.1.fq --pe1-2 SRR5660030.1_mapped.2.fq --pe2-1 SRR5660033.1_mapped.1.fq --pe2-2 SRR5660033.1_mapped.2.fq --pe3-1 SRR5660044.1_mapped.1.fq  --pe3-2 SRR5660044.1_mapped.2.fq --pe4-1 SRR5660045.1_mapped.1.fq --pe4-2 SRR5660045.1_mapped.2.fq -o HCMVassembly/")
outfile.write("\nspades -k 127 -t 2 --only-assembler --pe1-1 SRR5660030.1_mapped.1.fq --pe1-2 SRR5660030.1_mapped.2.fq --pe2-1 SRR5660033.1_mapped.1.fq --pe2-2 SRR5660033.1_mapped.2.fq --pe3-1 SRR5660044.1_mapped.1.fq  --pe3-2 SRR5660044.1_mapped.2.fq --pe4-1 SRR5660045.1_mapped.1.fq --pe4-2 SRR5660045.1_mapped.2.fq -o HCMVassembly/\n")

# 6 & 7 count contigs with length >1000 and add up their lengths

contigfile = open('HCMVassembly/contigs.fasta', 'r')

contigList = [] #list for storing sequences with more than 1000 bp

contigNum = 0 #counts number of contigs with length >1000
totalLength = 0 #sums up lengths of those contigs
for record in SeqIO.parse(contigfile, "fasta"):
    if len(record.seq) > 1000: #checks length of sequence
        contigNum += 1 #increases count of number of contigs longer than 1000 bp
        totalLength += len(record.seq) #adds to sum total of bp
        contigList.append(record.seq) #adds sequence to list

        
outfile.write("There are " + str(contigNum) + " contigs > 1000 bp in the assembly.\n")
outfile.write("There are " + str(totalLength) + " bp in the assembly.\n")

#8 run blast

sortedSeq = sorted(contigList, key = len, reverse = True) #sorts the sequence list by length longest to shortest

LC = open('largestContig.fasta', 'w')
LC.write(str(sortedSeq[0])) #stores the first index which is the longest sequence
LC.close()

handle = Entrez.esearch(db = "nucleotide", term = "txid10357[Organism:exp]", retmax = 20000) #search for betaherps seqs but only refseq ones
entry = Entrez.read(handle)

ids = entry.get('IdList') #saves ids of the search results

#retrieve fasta sequences for those ids
handle = Entrez.efetch(db = "nucleotide", id = ids, rettype = 'fasta')
entries = list(SeqIO.parse(handle, 'fasta'))

#writes seqs in fasta format 
betaherps = open('betaherpsSeqs.fasta', 'w')
for entry in entries:
    betaherps.write(entry.format('fasta'))
betaherps.close()

os.system("makeblastdb -in betaherpsSeqs.fasta -out HCMVblast -title HCMVblast -dbtype nucl") #makes db from betaherp sequences

#blasts the largest contig against the betaherps db
os.system('blastn -query largestContig.fasta -db HCMVblast -out HCMVblastResults.txt ' + ' -outfmt "6 sacc pident length qstart qend sstart send bitscore evalue stitle"') 

blast = open('HCMVblastResults.txt', 'r')

outfile.write('sacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tbitscore\tevalue\tstitle\n') #writes out headers with tab spacing

for line in range(10): #writes out top 10 blast results to log file
    outfile.write(str(blast.readline()) + '\n')

outfile.close()
