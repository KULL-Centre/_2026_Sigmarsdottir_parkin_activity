
module load python/3.9.9
module load eautils

ls ../fastq_250627/*_R1_001.fastq.gz  |  parallel -P 16 ./call_zerotol_paired.sh {};

Rscript merge_and_map.r samples.csv > merge_and_map.out
