#!/bin/bash

source paths.sh

if [ ! -d ${OUTPUT_FOLDER} ]; then
  mkdir ${OUTPUT_FOLDER}
fi

if [ $# -eq 0 ]; then
  echo "P-value threshold must be supplied as an argument"
  exit
fi

PVALTHR=$1
OUTPUT_PREFIX=${eQTL%*txt.gz}
OUTPUT_SUFFIX="signif_SNPs_pval_lt_${PVALTHR}.db"
OUTPUT_FILE=${OUTPUT_FOLDER}/${OUTPUT_PREFIX}_${OUTPUT_SUFFIX}
MAX_N_OF_SNPS=$1

python $PY --input_file ${eQTL} \
           --snp_column variant_id \
           --snp_dictionary ${SNP_DICTIONARY} \
           --gene_column gene_id \
           --pval_threshold ${PVALTHR} \
           --max_n_of_snps ${MAX_N_OF_SNPS} \
           --beta_column slope \
           --pvalue_column pval_nominal \
           --output_file ${OUTPUT_FILE} 
