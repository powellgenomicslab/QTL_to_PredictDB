import logging
import gzip
import pandas
import numpy
import Exceptions

from Constants import SNP
from Constants import REF_ALLELE
from Constants import ALT_ALLELE
from Constants import VAR_ID
from Constants import CHROMOSOME
from Constants import GENENAME
from Constants import POSITION
from Constants import BETA
from Constants import PVALUE

COLUMN_SNP="column_snp"
COLUMN_REF_ALLELE="column_effect_allele"
COLUMN_ALT_ALLELE="column_non_effect_allele"
COLUMN_CHROMOSOME="column_chromosome"
COLUMN_POSITION="column_position"
COLUMN_FREQ="column_freq"
COLUMN_BETA="column_beta"
COLUMN_PVALUE="column_pvalue"
COLUMN_SE="column_se"


def add_column_arguments_to_parser(parser):
    parser.add_argument("--gtex_format", help="Boolean indicating whether SNP column has the format chr_pos_ref_alt_genomebuild", action="store_true", default=False)
    parser.add_argument("--snp_dictionary", help="Path to file containing chr:position in the first column, and rsID in the second.")
    parser.add_argument("--ignore_if_not_in_dictionary", help="If true (default), discard SNP when it's not in the SNP dictionary provided.", default=True)
    parser.add_argument("--ignore_indels", help="If true (default), insertions and deletions are discarded.", default=True)
    parser.add_argument("--snp_column", help="Name of -snp column- in eQTL input file", default="SNP")
    parser.add_argument("--ref_allele_column", help="Name of -reference allele column- in eQTL input file", default=None)
    parser.add_argument("--alt_allele_column", help="Name of -alternative, or dosage, allele column- in eQTL input file", default=None)
    parser.add_argument("--gene_column", help="Name of -gene column- in eQTL input file", default=None)
    parser.add_argument("--chromosome_column", help="Name of -chromosome column- in eQTL input file", default=None)
    parser.add_argument("--position_column", help="Name of -base position column- in eQTL input file", default=None)
    parser.add_argument("--freq_column", help="Name of -frequency column- in eQTL input file", default=None)
    parser.add_argument("--beta_column", help="Name of snp association's -beta column- in eQTL input file", default=None)
    parser.add_argument("--se_column", help="Name of snp association's -beta standard error- column in eQTL input file", default=None)
    parser.add_argument("--pvalue_column", help="Name of snp association's -pvalue column- in eQTL input file", default=None)


#def format_dict_from_params(dict, parameters):

    #dict["snpID_in_gtex_format"] = parameters.gtex_format
    #if parameters.snp_column: dict[GWAS.COLUMN_SNP] = parameters.snp_column
    #if not parameters.gtex_format:
        #if hasattr(parameters, 'chromosome_column') and parameters.chromosome_column: dict[GWAS.COLUMN_CHROMOSOME] = parameters.chromosome_column
        #if hasattr(parameters, 'position_column') and parameters.position_column: dict[GWAS.COLUMN_POSITION] = parameters.position_column
        #if parameters.ref_column: dict[GWAS.COLUMN_REF_ALLELE] = parameters.ref_column
        #if parameters.alt_allele_column: dict[GWAS.COLUMN_ALT_ALLELE] = parameters.alt_allele_column
    #if hasattr(parameters, 'freq_column') and parameters.freq_column: dict[GWAS.COLUMN_FREQ] = parameters.freq_column
    #if parameters.beta_column: dict[GWAS.COLUMN_BETA] = parameters.beta_column
    #if parameters.se_column: dict[GWAS.COLUMN_SE] = parameters.se_column
    #if parameters.pvalue_column: dict[GWAS.COLUMN_PVALUE] = parameters.pvalue_column


def _f_snp(format): return format[COLUMN_SNP] if COLUMN_SNP in format else None
def _f_ref_allele_column(format): return format[COLUMN_REF_ALLELE] if COLUMN_REF_ALLELE in format else None
def _f_alt_allele_column(format): return format[COLUMN_ALT_ALLELE] if COLUMN_ALT_ALLELE in format else None
def _f_pvalue_column(format): return format[COLUMN_PVALUE] if COLUMN_PVALUE in format else None
def _f_beta_column(format): return format[COLUMN_BETA] if COLUMN_BETA in format else None
def _f_se_column(format): return format[COLUMN_SE] if COLUMN_SE in format else None


# snp, p-value and beta columns are mandatory.
# reference and alternate allele are also mandatory, unless the SNP ID has this information embedded (as is the case of GTEx's varID's.)

def prevalidate_format(params, header_comps):
    valid_format = True
    if params.pvalue_column not in header_comps:
        logging.error("P-value column name is not correct.")
        valid_format = False
    if params.beta_column not in header_comps:
        logging.error("Beta column name is not correct.")
        valid_format = False
    if params.snp_column not in header_comps:
        print params.snp_column, header_comps
        logging.error("SNP column name is not correct.")
        valid_format = False
    if params.gene_column not in header_comps:
        logging.error("Gene column name is not correct.")
        valid_format = False
    if not valid_format:
        exit()


def validate_input_parameters(params, info_from_snpid):
    valid_format = True
    if info_from_snpid == 1 or info_from_snpid == 2:
        if ref_allele_column is None:
            logging.error("You have to provide a reference allele column.")
            valid_format = False
        if alt_allele_column is None:
            logging.error("You have to provide an alternate allele column.")
            valid_format = False
    if not valid_format:
        exit("Aborting...")


  #def validate_format(format):
    #if not format.gtex_format:
        #if not _f_ref_allele_column(format): raise Exceptions.InvalidArguments("Need to provide a -reference allele- column")
        #if not _f_alt_allele_column(format): raise Exceptions.InvalidArguments("Need to provide an -alternative (or dosage) allele- column")
    #if not _f_snp(format): raise Exceptions.InvalidArguments("Need to provide a SNP column")
    #if not _f_snp(format): raise Exceptions.InvalidArguments("Need to provide a SNP column")
    #if not _f_pvalue_column(format): raise Exceptions.InvalidArguments("Need to provide a -p-value- column")
    #if not _f_beta_column(format): raise Exceptions.InvalidArguments("Need to provide a -beta- column")


def write_row_into_db(row, db):
    db("INSERT INTO weights VALUES(\"%s\", \"%s\", \"%s\", \"%s\", \"%s\", \"%s\", \"%s\", NULL, NULL)" % \
        (row[SNP], row[GENENAME], row[BETA], row[VAR_ID], row[REF_ALLELE], row[ALT_ALLELE], row[PVALUE]))


def write_df_into_db(df, db):
    for index, row in df.iterrows():
        write_row_into_db(row, db)


# Replace decimal separator "," by "." if necessary, and NA and "." for None
import re
non_en_number = re.compile("^[-\+]?[0-9]*,{1}[0-9]+([eE]{1}[-\+]?[0-9]+)?$")
def sanitize_component(c):
    if non_en_number.match(c): c = c.replace(",",".")
    if c == "NA": c = None
    if c == ".": c = None
    return c
