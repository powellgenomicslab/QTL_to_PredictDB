#!/bin/bash

export PATH_PREFIX=".."
export PY="Main.py"
export PY=${PATH_PREFIX}/${PY}

export OUTPUT_FOLDER="${PATH_PREFIX}/OutputModels"

export SNP_DICTIONARY="${PATH_PREFIX}/snp_data/gtex_v7_hapmapceu_dbsnp150_snp_annot.txt.gz"
export eQTL="${PATH_PREFIX}/eQTL_SampleData/Whole_Blood_allpairs_truncated.txt.gz"
