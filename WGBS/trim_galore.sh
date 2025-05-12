# Paired-end sequencing
trim_galore --paired --gzip SRRxxx_1.fastq SRRxxx_2.fastq --output_dir ./Trim_galore/

# Single-Read Sequencing
trim_galore --gzip SRRxxx.fastq --output_dir ./Trim_galore/