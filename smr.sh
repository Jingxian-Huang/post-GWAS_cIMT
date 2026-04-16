#!/bin/bash
#PBS -l select=1:ncpus=32:mem=64gb
#PBS -l walltime=24:00:00
#PBS -o output.log
#PBS -e error.log

smr --bfile mydata --gwas-summary mygwas.ma --beqtl-summary myeqtl --out mysmr --thread-num 10