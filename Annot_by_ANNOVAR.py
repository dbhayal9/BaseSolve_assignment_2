import os, re, sys
import shutil
import gzip
from collections import defaultdict
from collections import Counter
import glob

print("All modules imported")
print(os.getcwd())
print(os.chdir("4.Variants/"))

f = os.listdir()
f
files = glob.glob("*.final.filtered.vcf.gz")
files
files = [i.rstrip(".final.filtered.vcf.gz") for i in files]
files
humndbA = "/media/hps2/Download/annovar/humandb/"
for i in files:
    print("Annoation by ANNOVAR for sample ID: ",i)
    f = i + ".final.filtered.vcf.gz"
    print("Annoation by ANNOVAR for sample ID: ",f)
    os.system(f"perl /media/hps2/Download/annovar/./table_annovar.pl {f} {humndbA} -buildver hg19 -out ../5.Annotation/{i}.txt -remove -protocol refGene,esp6500siv2_all,1000g2015aug_all,avsnp147,dbnsfp35a,clinvar_20190305,gnomad_genome,dbscsnv11,rmsk,ensGene,knownGene,cytoBand,exac03,cosmic91 -operation g,f,f,f,f,f,f,f,r,g,g,r,f,f -nastring . -vcfinput")

