#from wdb_main_v2 import main
from wdb_main import main

main("eQTL_SampleData/Whole_Blood_allpairs_truncated.txt", \
     "snps_data/gtex_v7_hapmapceu_dbsnp150_snp_annot.txt.gz", \
     "variant_id", "gene_id", "slope", 1e-2, "pval_nominal")
