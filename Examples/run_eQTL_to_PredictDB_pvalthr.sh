#!/bin/bash

source paths.sh

if [ ! -d ${OUTPUT_FOLDER} ]; then
  mkdir ${OUTPUT_FOLDER}
fi

if [ $# -eq 0 ]; then
  echo "P-value threshold must be supplied as an argument"
  exit
fi

eQTL=$1
PVALTHR=$2
OUTPUT_PREFIX=${eQTL%*txt.gz}
OUTPUT_PREFIX=${OUTPUT_PREFIX##*/}
OUTPUT_SUFFIX="signif_SNPs_pval_lt_${PVALTHR}.db"
OUTPUT_FILE=${OUTPUT_FOLDER}/${OUTPUT_PREFIX}_${OUTPUT_SUFFIX}

echo $OUTPUT_FILE

python $PY --input_file ${eQTL} \
           --snp_column variant_id \
           --snp_dictionary ${SNP_DICTIONARY} \
           --gene_column gene_id \
           --pval_threshold ${PVALTHR} \
           --beta_column slope \
           --pvalue_column pval_nominal \
           --output_file ${OUTPUT_FILE} 
