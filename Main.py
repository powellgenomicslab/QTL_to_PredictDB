#!/usr/bin/env python

'''
  This tool allows to convert the result of multiple genome-wide univariate regressions
  on a molecular trait (QTLs) [typically gene expression (eQTL), splicing (sQTL),
  protein levels (pQTL)] to a SQLite database containing prediction weights, suitable
  for use with the PrediXcan/S-PrediXcan tools.

  To be implemented: produce a genotype covariance matrix associated to that model, a 
  necessary input for S-PrediXcan.
'''

__author__ = 'rbonazzola'
__version__ = '0.0'

import os, sys, time
import SQLiteDB

import DataFrameStreamer as DFS
import Utilities
import Methods
import logging
import Constants
import pandas
import re

DEBUG = False
N_OF_GENES = 200

########################################################################################################################
#                                                |--------------------|                                                #
#------------------------------------------------|    MAIN FUNCTION   |------------------------------------------------#
#                                                |--------------------|                                                #
########################################################################################################################


def check_if_dbfile_exists(args, batch, ask_user, remove_existent_db):
    # Ask the user how to proceed in the presence of a previously existent homonym DB.
    if os.path.exists(args.output_file):
        if not remove_existent_db:
            if not batch and ask_user:
                overwrite = "." # some initial value different from Y/y/N/n
                while overwrite not in "YyNn":
                    sys.stdout.write("A database called %s already exists. Should I overwrite it? (Y/y/N/n): " % os.path.basename(args.output_file))
                    overwrite = raw_input()
                    sys.stdout.flush()
                if overwrite.lower() == "y":
                    os.remove(args.output_file) # FIXME: or should I copy it into in a temporal file until the new DB is successfully created?
                else:
                    exit("Exiting...")
            else:   # if batch is True, we don't want the user to write input in run time
                exit("A database called %s already exists. Change output DB's name or remove old DB." % os.path.basename(args.output_file))
        else:
           print "Removing old %s file" % os.path.basename(args.output_file)
           os.remove(args.output_file)


def build_white_and_black_gene_lists(args):
    
    gene_white_list, gene_black_list = (None, None)

    if args.gene_white_list is not None:
        # FIXME: maybe I should remove the restriction that the file has a "txt" extension
        if args.gene_white_list.endswith(".txt") and args.gene_white_list in os.listdir("."):
            gene_white_list = []
            with open(args.gene_white_list) as gene_wl:
                for line in gene_wl:
                    gene_white_list.extend([gene for gene in line.strip().split()])

    if args.gene_black_list is not None:
        if args.gene_white_list is not None:
            logging.error("You cannot have both a white list and a black list for genes.")
            exit("Aborting...")
        if args.gene_black_list.endswith(".txt") and args.gene_black_list in os.listdir("."):
            gene_black_list = []
            with open(args.gene_black_list) as gene_bl:
                for line in gene_bl:
                    gene_black_list.extend([gene for gene in line.strip().split()])

    return gene_white_list, gene_black_list



def run(args, batch = False, ask_user = False, remove_existent_db=True):

    start = time.time()
   
    check_if_dbfile_exists(args, batch, ask_user, remove_existent_db)
    
    gene_white_list, gene_black_list = build_white_and_black_gene_lists(args)

    db = SQLiteDB.DB(args.output_file)

    genes_total = 0; genes_passed = 0

    # files to iterate over
    if args.input_file:
        eqtl_files = [args.input_file]
    elif args.input_folder:
        eqtl_files = [os.path.join(args.input_folder, x) for x in os.listdir(args.input_folder)]
    else:
        logging.error("You must provide exactly one of the followings arguments: --input_folder or --input_file.")
        exit("Aborting...")

    if args.snpid_format != "rsid":
       snpid_regex = re.compile(args.snpid_format)
    else:
       snpid_regex = None

    # Each chunk corresponds to all SNPs associated to the expression of a given gene
    for eqtl_file in eqtl_files:
        logging.info("Processing %s" % eqtl_file)
        for i, chunk in enumerate(DFS.read_eQTL_file(eqtl_file, args, snpid_regex)):
            chunk = Methods.generate_weights(chunk, args) # convert to a Pandas data frame
            if len(chunk.index) != 0:
                genes_passed += 1
                gene_name = chunk.iloc[0][Constants.GENENAME] # extract gene name
                db("INSERT INTO extra VALUES(\"%s\", \"%s\",  NULL, NULL, NULL)" % (gene_name, gene_name))
            Utilities.write_df_into_db(chunk, db)
            genes_total = i + 1
            if DEBUG and i == N_OF_GENES:
                break

    print "DB %s created, with %d/%d genes." % (args.output_file, genes_passed, genes_total)
    print "--- %s seconds ---" % (time.time() - start)

    db.create_indexes()
    db.close()


########################################################################################################################


if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(description='Create a PrediXcan-compatible SQLite database with gene expression prediction weights computed from eQTL input files.')

    parser.add_argument("--input_file", help="Path to input eQTL file.")
    parser.add_argument("--input_folder", help="Folder containing eQTL files.")

    # add to parses arguments concerning file column names
    Utilities.add_column_arguments_to_parser(parser)

    parser.add_argument("--output_file", "-o", default="output.db",
                        help="SQLite DB output file")

    # FIXME: change names of methods.
    # I want to call "naive" the method that takes the betas as weights, regardless of the number of SNPs or the p-value threshold.
    # univar_to_multivar -> multivariate
    weights_method_lst = ["top_eQTL", "signif_eQTL", "univar_to_multivar", "BLUP", "BSLMM", "LDPred", "LASSO"]

    parser.add_argument("--method", choices = weights_method_lst, default="naive",
                        help="Method for obtaining weights from eQTL's betas")

    # FIXME: validate that it's greater than zero and less than or equal to 1.
    parser.add_argument("--pval_threshold", default=None,
                        help="P-value above which SNPs are discarded to compute the prediction weights.")

    parser.add_argument("--max_n_of_snps", default=None,
                        help="Maximum number of SNPs to be considered for the model of a given gene.")

    # FIXME: they may not be organized by chromosome.
    parser.add_argument("--covariance_matrix",
                        help="Folder containing SNP covariance matrices, organized by chromosome. (Not used yet)")

    parser.add_argument("--genotype_LD_reference",
                        help = """Folder containing dosage files for a reference panel, in PrediXcan format.
                                These files are used to compute the covariance matrices for the SNPs involved in each gene's prediction model.
                                This parameter and --covariance_matrix are mutually exclusive. (Not used yet)""",
                        default = None)

    parser.add_argument("--gene_white_list",
                        help="White list of genes to generate models for. It can be the path for a txt file with the names of the genes, or the names themselves.", default = None)

    parser.add_argument("--gene_black_list",
                        help="Idem as white list, but for genes to exclude.", default = None)

    parser.add_argument("--logs_file", default="log.txt",
                        help="Name of the file with the error output")

    args = parser.parse_args()

    run(args)
