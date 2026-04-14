#!/bin/bash
#PBS -l select=1:ncpus=8:mem=200gb
#PBS -l walltime=24:00:00
#PBS -o output.log
#PBS -e error.log

# Specify directory
WORK_DIR=""
sumsta="/gwas/imt_eur.txt"
out="/twas"


# Munge first
module load anaconda3/personal
cd ~/ldsc
source activate ldsc

./munge_sumstats.py \
  --sumstats $sumsta \
  --N-col N  \
  --out "${out}/munge" \
  --a1 ALLELE1 \
  --a2 ALLELE2 \
  --p P-value \
  --frq Freq1 \
  --snp MarkerName \
  --chunksize 500000

gunzip "${out}/munge.sumstats.gz"

# Use for loop to run from chromosome 1 to 22 for each tissue-specific expression weights
cd $WORK_DIR
source activate r422

for i in {1..22}; do

  for file in ./weight_imt/*.pos; do

    name=$(basename "$file" .pos)

    Rscript FUSION.assoc_test.R \
    --sumstats "${out}/munge.sumstats" \
    --weights $file \
    --weights_dir ./weight_imt/ \
    --ref_ld_chr /UKB/UKB_ref_panel \
    --chr $i \
    --coloc_P 0.05 \
    --GWASN 11737 \
    --out "${out}/res/${name}.${i}"

  done

done