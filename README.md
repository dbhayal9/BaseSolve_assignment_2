## This is the pipeline to analyse WES data from fastq to vcf annotation

### Getting raw data from NCBI SRA DB
##### sample-1 SRR25234909 run accession ID
##### sample-2 SRR16896259 run accession ID
### NOTE: download raw data script name: rawdata.sh (uploaded in main)
##### wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR252/009/SRR25234909/SRR25234909_1.fastq.gz
##### wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR252/009/SRR25234909/SRR25234909_2.fastq.gz
##### wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR168/059/SRR16896259/SRR16896259_1.fastq.gz
##### wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR168/059/SRR16896259/SRR16896259_2.fastq.gz

### Getting Reference Human Genome HG19 from UCSC
##### wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
### To calculate BAM file alignment summery command is below:
##### samtools flagstat sample1.sorted.bam > sample1_align_summery.txt
##### samtools flagstat sample2.sorted.bam > sample2_align_summery.txt
### NOTE: To process from fastq to VCF please see python script name : fastq_vcf_WES.py (uploaded in main)

### To identify common varinats seperately SNVs and INDELs command below:
#####     bedtools intersect -a sample1.indels_filtered.vcf.gz -b sample2.indels_filtered.vcf.gz > common.indel.vcf.gz
#####     bedtools intersect -a sample1.snps_filtered.vcf.gz -b sample2.snps_filtered.vcf.gz > common.snps.vcf.gz

### For Annotation or final variants by ANNOVAR tools script name : Annot_by_ANNOVAR.py (uploaded in main)

