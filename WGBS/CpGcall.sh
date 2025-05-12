## Deduplication
deduplicate_bismark --bam SRRxxx_1_val_1_bismark_bt2_pe.bam --output_dir ./

## Filtering out reads with non-CG methylation 
# single-end file
filter_non_conversion -s SRRxxx_1_val_1_bismark_bt2_pe.deduplicated.bam
# paired-end file
filter_non_conversion -p SRRxxx_1_val_1_bismark_bt2_pe.deduplicated.bam

## Methylation extraction
# single-end file:
bismark_methylation_extractor --multicore 10 -s --comprehensive --merge_non_CpG  --bedGraph --counts --buffer_size 10G --gzip SRRxxxcopy_trimmed_bismark_bt2.deduplicated.nonCG_filtered.bam --output /users/zhengd/scratch/Work2/Population  
# paired-end file:
bismark_methylation_extractor --multicore 10 -p --comprehensive --merge_non_CpG --bedGraph --counts --buffer_size 10G --no_overlap --gzip file_1.fq_bismark_pe.bam --output /users/zhengd/scratch/Work2/Population

coverage2cytosine --merge_CpG --genome_folder /home/zdong/Referencedata/hg19_NCBI -o SRR3274347 SRRxxx_1_val_1_bismark_bt2_pe.deduplicated.nonCG_filtered.bismark.cov.gz