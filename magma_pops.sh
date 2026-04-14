#!/bin/bash
#PBS -l select=1:ncpus=24:mem=64gb
#PBS -l walltime=24:00:00
#PBS -o output.log
#PBS -e error.log

# Specify directory
WORK_DIR=""
sumsta="/gwas/imt_eur.txt"
out="/pops"

module load anaconda3/personal
source activate py39

# MAGMA
cd $WORK_DIR

./magma \
--bfile /UKB/UKB_ref_panel \
--gene-annot hg37.genes.annot \
--pval $sumsta use=MarkerName,9 ncol=N \
--gene-model snp-wise=mean \
--out "${out}/magma.out"

# POPS
python pops.py \
--gene_annot_path example/data/utils/gene_annot_jun10.txt \
--feature_mat_prefix pops_features \
--num_feature_chunks 44 \
--magma_prefix "${out}/magma.out" \
--out_prefix "${out}/pops.out"