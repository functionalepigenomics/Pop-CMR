## Single-end alignments
bismark --bowtie2 --multicore 32 --gzip -n 1 /home/zdong/Referencedata/hg19_NCBI SRRxxx_trimmed.fq --output_dir /users/zhengd/scratch/Work2/Population


## Paired-end alignments
bismark --bowtie2 --multicore 32 --gzip -n 1 /home/zdong/Referencedata/hg19_NCBI  -1 file_1.fq -2 file_2.fq --output_dir /users/zhengd/scratch/Work2/Population