#!/usr/bin/env python
# coding: utf-8

# In[112]:

##conda activate igib
### Whole Exome Analysis (WES)Pipeline
### To Analysis Illumina WES data from fastq to Annotation
### Requirements:
### Conda env
### FastQC
### Trimgalore/fastp
### BWA
### samtools
### bcftools
### picardtools (outside from env)
### GATK (build)
### ANNOVAR

import os, re, sys
import shutil
import gzip
from collections import defaultdict
from collections import Counter
import glob

print("All modules imported")


# In[113]:

"""
print(os.getcwd())
print(os.chdir("/media/hps2/Pipeline_2023/client/WESII/"))
print(os.getcwd())
print(os.listdir())


# In[ ]:

"""
### Creating directorie's
dirs = ["1.RawData", "2.TrimmedData", "3.Alignment", "4.Variants", "5.Annotation"]
print("Creating directories......")
try:
    for i in dirs:
        os.mkdir(i)
    #print(os.listdir())
except:
    print("Directory already Exist")
print(os.listdir())
print(os.mkdir("1.RawData/QC"))
print(os.listdir("1.RawData"))


# In[ ]:


### Get data into 1.RawData dir
source_folder = "/media/bioinfo/sda1/basesolev/"
destination_folder = "/media/bioinfo/sda1/basesolev/1.RawData/"

# fetch all files
for file_name in os.listdir(source_folder):
    # construct full file path
    source = source_folder + file_name
    destination = destination_folder + file_name
    # move only files
    if os.path.isfile(source):
        shutil.move(source, destination)
        print('Moved:', file_name)

print("Data fetched")
# In[ ]:


### Quality checking of raw data ## Worked while run from terminal
print(os.chdir("1.RawData"))
print(os.getcwd())
os.system(f"fastqc *.gz -t 4 -o QC")

print("QC Done....")
# In[ ]:


### MultiQC for Raw data
print(os.getcwd())
print(os.chdir("QC/"))
print(os.getcwd())
print(os.listdir())
os.system(f"multiqc . -p")

print("MultiQC Done....")
# In[ ]:


### Trimming of raw data using TrimGalore ###Fastp (optional)
print(os.getcwd())
print(os.chdir("../../1.RawData/"))
print(os.listdir())
#print(os.getcwd())
print("Triming of raw data...............")
os.system("parallel --ungroup --link --xapply trim_galore --illumina --trim-n --paired -o ../2.TrimmedData/ ::: *_R1.fastq.gz ::: *_R2.fastq.gz") 

print("Trimminig Done....")
# In[ ]:


### Indexing for Reference genome hg19

print(os.getcwd())
print(os.chdir("../reference/"))
print(os.getcwd())
print(os.listdir())
print("Indexing for reference genome hg19..........")
os.system("bwa index hg19.fa hg19")
print("Done indexing")

##### Samtools create .fai file for reference genome
print("indexing using samtools for reference geneome.........")
os.system(f"samtools faidx hg19.fa")

##### Picard tool create .dict file for reference genome
print("########################### making dict file ##############################################")
os.system(f"java -jar /home/bioinfo/Downloads/picard.jar CreateSequenceDictionary \R=hg19.fa \O=hg19.dict")

print("########################### Done dict file ##############################################")
# In[ ]:


### Alignment of raw trimmed fastq files againest reference genome
print(os.getcwd())
print(os.chdir("../2.TrimmedData/"))
#print(os.chdir("../2.TrimmedData/")) for temporary
print(os.getcwd())
files = os.listdir()
print(files)
refere = "/media/bioinfo/sda1/basesolev/reference/hg19.fa"
#files = os.listdir()
threads = 4

### worked for multiple samples
files = glob.glob("*.gz")
files = [i.rstrip("_val_{1,2}.fq.gz") for i in files]
files = [i.rstrip("_R") for i in files]

#print(*Counter(files))
files1 = Counter(files)
for f in files1:
    print("Alignment for sample ID: ",f)
    print("Aligning....................")
    f1 = f + "_R1_val_1.fq.gz"
    f2 = f + "_R2_val_2.fq.gz"
    os.system(f"bwa mem -t {threads} -M -R '@RG\\tID:{f}\\tPL:ILLUMINA\\tPM:HISEQ\\tSM:{f}' {refere}  {f1} {f2} | samtools view -@ 4 -b - | samtools sort -@ 4 > ../3.Alignment/{f}.sorted.bam")
    #cmd = "bwa mem -t 4 -M /media/hps2/Pipeline_2023/WES/reference/hg19.fa " + f1 + " " + f2+ " > " + f + ".sam"
    #os.system("bwa mem -t 4 -M /media/hps2/Pipeline_2023/WES/reference/hg19.fa " + f1 + " " + f2+ " > " + f + ".sam")
    #Worked# os.system("bwa mem -t 4 -M /media/hps2/Pipeline_2023/WES/reference/hg19.fa " + f1 + " " + f2+ " > " + f + ".sam")
    #Worked# os.system("bwa mem -t 4 -M /media/hps2/Pipeline_2023/WES/reference/hg19.fa " + f1 + " " + f2 + " | samtools view -@ 4 -b - | samtools sort -@ 4 > " + f + ".sorted.bam")
    #os.system("bwa mem -t {threads} -M /media/hps2/Pipeline_2023/WES/reference/hg19.fa 41247351_R1_val_1.fq.gz 41247351_R2_val_2.fq.gz > ../3.Alignment/41247351.sam") ### for one sample its working

print("########################### Alignment DONE ##############################################")
# In[ ]:


### MarkDuplicate

### NOTE: should include REMOVE_DUPLICATES=true option

print(os.getcwd())
print(os.chdir("../3.Alignment/"))
print(os.getcwd())
files = os.listdir()
files = glob.glob("*.sorted.bam")
print(files)
files = [i.rstrip(".sorted.bam") for i in files]
for i in files:
    print("Running for sample ID: ",i)
    f1 = i + ".sorted.bam"
    print("Running for sample ID: ",f1)
    os.system(f"java -jar /home/bioinfo/Downloads/picard.jar MarkDuplicates \I={f1} \O={i}.marked_duplicates.bam \M={i}.marked_dup_metrics.txt")

print("############################### MarkDuplicate Done ###############################################....")
# In[ ]:


### Indexing for final BAM files
print(os.getcwd())
files = glob.glob("*.sorted.bam")
for i in files:
    print("Sorting sample ID: ",i)
    os.system(f"samtools index {i}")

print("####################################### Indexing for BAM files Done ##################################....")
# In[ ]:



### Filtered reads from BAM files using sambamba
### in proccess


# In[ ]:


### HaplotypeCaller
### Call variants from BAM files
print("################################# start variant calling #########################################################")
print(os.getcwd())
#refe = "/media/hps2/Pipeline_2023/WES/reference/hg19.fa" ## fatching from BWA part
BEDIllumina = "/media/bioinfo/sda1/basesolev/BED/hg19_Twist_ILMN_Exome_2.0_Plus_Panel_annotated.BED"
#dbsnp = "/media/hps2/Pipeline_2023/WES/dbsnp/common_all_20180423.vcf"
#print(os.chdir("3.Alignment/"))
files = os.listdir()
files = glob.glob("*.sorted.bam")
print(files)
files = [i.rstrip(".sorted.bam") for i in files]
files
for i in files:
    print("Running for sample: ",i)
    f1 = i + ".sorted.bam"
    print("Running for sample: ",f1)
    ##os.system(f"java -jar /media/hps2/Download/gatk-4.2.5.0/gatk-package-4.2.5.0-local.jar HaplotypeCaller \-R {refe} \-I {f1} \-L {BEDIllumina} \-D {dbsnp} \-O ../4.Variants/{f1}.vcf.gz \-bamout {f1}.bamout.bam")
    os.system(f"/home/bioinfo/Downloads/gatk/./gatk HaplotypeCaller \-R {refere} \-I {f1} \-L {BEDIllumina} \-O ../4.Variants/{f1}.vcf.gz \-bamout ../4.Variants/{f1}.bamout.bam")

print("##################################### HaplotypeCaller Done ###################################################....")
# In[ ]:


### Hard-filterning 1. Select variants 2. Variants filteration 3. merge Indel+SNP > final filterd vcf
### select variants
print("##################################### Select Variant start ###################################################....")
print(os.getcwd())
print(os.chdir("../4.Variants/"))
#print(os.mkdir("SelectVariants"))
print(os.listdir())
files = os.listdir()
files = glob.glob("*.gz")
files = [i.rstrip(".sorted.bam.vcf.gz") for i in files]
files = [i.rstrip(".sorted.bam.vcf.gz.tbi") for i in files]

print(files)
files1 = Counter(files)
#files1
### for SNP
for i in files1:
    print("running for: ",i)
    f1 = i + ".sorted.bam.vcf.gz"
    print("running for: ",f1)
    os.system(f"/home/bioinfo/Downloads/gatk/./gatk SelectVariants  \-V {f1} \-select-type SNP \-O {i}.vcf.snp.gz")
    

### for INDEL
for i in files1:
    print("running for: ",i)
    f1 = i + ".sorted.bam.vcf.gz"
    print("running for: ",f1)
    os.system(f"/home/bioinfo/Downloads/gatk/./gatk SelectVariants  \-V {f1} \-select-type INDEL \-O {i}.vcf.indel.gz")

print("##################################### Select Variant Done ###################################################....")
# In[ ]:


### Hard-filterning 1. Select variants 2. Variants filteration 3. merge Indel+SNP > final filterd vcf
### Filter variants
print("##################################### Filter Variant Start ###################################################....")
print(os.getcwd())
files = os.listdir()
files = glob.glob("*.vcf.snp.gz")
filesin = glob.glob("*.vcf.indel.gz")

#files
#files1 = glob.glob("*.vcf.snp.gz")
#files1
files1 = [i.rstrip(".vcf.snp.gz") for i in files]
filesin1 = [i.rstrip(".vcf.indel.gz") for i in filesin]

#files1
for i in files1:
    f1 = i + ".vcf.snp.gz"
    #print(f1)
    os.system(f"/home/bioinfo/Downloads/gatk/./gatk VariantFiltration \-V {f1} \-filter 'QD < 2.0' --filter-name 'QD2' \-filter 'QUAL < 30.0' --filter-name 'QUAL30' \-filter 'SOR > 3.0' --filter-name 'SOR3' \-filter 'FS > 60.0' --filter-name 'FS60' \-filter 'MQ < 40.0' --filter-name 'MQ40' \-filter 'MQRankSum < -12.5' --filter-name 'MQRankSum-12.5' \-filter 'ReadPosRankSum < -8.0' --filter-name 'ReadPosRankSum-8' \-O {i}.snps_filtered.vcf.gz")
for i in filesin1:
    f2 = i + ".vcf.indel.gz"
    os.system(f"/home/bioinfo/Downloads/gatk/./gatk VariantFiltration \-V {f2} \-filter 'QD < 2.0' --filter-name 'QD2' \-filter 'QUAL < 30.0' --filter-name 'QUAL30' \-filter 'FS > 200.0' --filter-name 'FS200' \-filter 'ReadPosRankSum < -20.0' --filter-name 'ReadPosRankSum-20' \-O {i}.indels_filtered.vcf.gz")

print("##################################### Filter Variant Done ###################################################....")
# In[ ]:


### Hard-filterning 1. Select variants 2. Variants filteration 3. merge Indel+SNP > final filterd vcf
### Merge variants
print("##################################### Merge Variant Start ###################################################....")
print(os.getcwd())
#print(os.chdir("4.Variants/"))
files = os.listdir()
files
files = glob.glob("*s_filtered.vcf.gz")
files
files = [i.rstrip(".indels_filtered.vcf.gz") for i in files]
files = [i.rstrip(".snp") for i in files]
files
files1 = Counter(files)

for i in files1:
    print("sample id: ",i)
    snv = i + ".snps_filtered.vcf.gz"
    ind = i + ".indels_filtered.vcf.gz"
    os.system(f"java -jar /home/bioinfo/Downloads/picard.jar MergeVcfs \I={snv} \I={ind} \O={i}.final.filtered.vcf.gz")
print("##################################### Merge Variant Done ###################################################....")

