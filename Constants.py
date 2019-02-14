GENE = "gene"
GENENAME = "gene_name"
SNP = "snp"
REF_ALLELE = "ref_allele"
ALT_ALLELE = "alt_allele"  # dosage allele
VAR_ID = "var_id"
CHROMOSOME = "chromosome"
POSITION = "position"
DOSAGES = "dosages"
FREQUENCY = "frequency"
BETA = "beta"
SE = "se"
PVALUE = "pvalue"


GENENAME_idx = 0
SNP_idx = 1
CHROMOSOME_idx = 2
POSITION_idx = 3
REF_ALLELE_idx = 4
ALT_ALLELE_idx = 5 # dosage allele
PVALUE_idx = 6
BETA_idx = 7
SE_idx = 8

CHROMOSOME_IDX_FROM_SNP_ANNOTATION = 0
POSITION_IDX_FROM_SNP_ANNOTATION = 1
REF_ALLELE_IDX_FROM_SNP_ANNOTATION = 3
ALT_ALLELE_IDX_FROM_SNP_ANNOTATION = 4
RSID_IDX_FROM_SNP_ANNOTATION = 6

alleles = ["A", "C", "T", "G", "a", "c", "t", "g"]

SNPID_1 = [SNP]
SNPID_2 = [CHROMOSOME, POSITION]
SNPID_3 = [CHROMOSOME, POSITION, REF_ALLELE, ALT_ALLELE]
