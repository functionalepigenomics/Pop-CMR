## Download GEO/SRA
prefetch --max-size 70G SRRxxx

## Convert SRA file to fastq
#single-end
fastq-dump  .../SRR***.sra -O out_path

#pair-end
fastq-dump --split-e  .../SRR***.sra -O out_path 