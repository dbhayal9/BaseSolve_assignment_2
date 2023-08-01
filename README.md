## This is the pipeline to analyse WES data from fastq to vcf annotation

### Getting raw data from NCBI SRA DB
##### sample-1 SRR25234909 run accession ID
##### sample-2 SRR16896259 run accession ID
##### wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR252/009/SRR25234909/SRR25234909_1.fastq.gz
##### wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR252/009/SRR25234909/SRR25234909_2.fastq.gz
##### wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR168/059/SRR16896259/SRR16896259_1.fastq.gz
##### wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR168/059/SRR16896259/SRR16896259_2.fastq.gz

### Getting Reference Human Genome HG19 from UCSC
##### wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
