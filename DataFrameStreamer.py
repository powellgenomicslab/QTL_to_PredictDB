import gzip
import pandas
import logging
import Utilities
import Constants
import os
import re

logging.basicConfig(level=logging.DEBUG, format="%(asctime)s: %(message)s")

def to_rsid(snp_dict, chr, pos):

    # logging.info("Loading SNPs definition.")
    if isinstance(rsid_dict, str):
        rsid_dict = load_snp_dictionary(rsid_dict)  # If rsid_dict is a string, it's assumed that it's the path to a file containing two columns (chr:pos and rsid).
    elif not isinstance(rsid_dict, dict):
        exit("Argument rsid_dict must be a path to a dictionary file or a python dictionary.")

    #logging.info("Dictionary loaded.")
    key = (chr, pos)

    try:
        rsid = rsid_dict[key]
        return rsid
    except KeyError:
        return None


def get_col_indexes(header_comps, params):
    # TODO: Assuming that "params" is correct, i.e. parameters passed as arguments do appear in the file header.
    # TODO: consider the case where the header is not the first line in the file, adding a skip_until_header argument)
    # header_comps = header.strip().split()  # column names
    col_indexes = {
        Constants.SNP: header_comps.index(params.snp_column),
        Constants.GENENAME: header_comps.index(params.gene_column),
        Constants.PVALUE: header_comps.index(params.pvalue_column),
        Constants.BETA: header_comps.index(params.beta_column),
    }
    if params.se_column is not None:
        col_indexes[Constants.SE] = header_comps.index(params.se_column)
    if params.ref_allele_column is not None:
        col_indexes[Constants.REF_ALLELE] = header_comps.index(params.ref_allele_column)
    if params.alt_allele_column is not None:
        col_indexes[Constants.ALT_ALLELE] = header_comps.index(params.alt_allele_column)
    if params.chromosome_column is not None:
        col_indexes[Constants.CHROMOSOME] = header_comps.index(params.chromosome_column)
    if params.position_column is not None:
        col_indexes[Constants.POSITION] = header_comps.index(params.position_column)

    return col_indexes


def to_dataframe(data, columns, to_numeric=None, fill_na=None):

    data = zip(*data)
    if to_numeric:
        data = [pandas.to_numeric(x, errors=to_numeric) for x in data]
    if len(data) == 0:
        data = [[] for i in range(0, len(columns))]
    data = { columns[i]: data[i] for i in range(0, len(columns)) }
    data = pandas.DataFrame(data)
    data = data[columns]
    if fill_na:
        data = data.fillna(fill_na)
    return data


def get_chr_pos_from_snpid(snpid, sep = "_"):
    chromosome, position = snpid.split(sep)
    return {Constants.CHROMOSOME: chromosome, Constants.POSITION: position}


def get_rsid(snpid):
    return {Constants.SNP: snpid}


def get_info_from_var_id(varID):
    chr, pos, ref_allele, alt_allele, build_version = varID.split("_")
    return {Constants.CHROMOSOME: chr, Constants.POSITION: pos, \
            Constants.REF_ALLELE: ref_allele, Constants.ALT_ALLELE: alt_allele}


def build_var_id(chromosome, position, ref_allele, alt_allele):
    return "_".join([str(chromosome), str(position), ref_allele, alt_allele, "b37"])


def build_dictionary(filepath, rsid_as_key = False):

    def chr_from_snps_file(comps):           chr = comps[0]; return chr
    def pos_from_snps_file(comps):           pos = comps[1]; return pos
    def ref_allele_from_snps_file(comps):    ref = comps[3]; return ref
    def alt_allele_from_snps_file(comps):    alt = comps[4]; return alt
    def rsid_from_snps_file(comps):         rsid = comps[7]; return rsid
    def generate_key(chr, pos):                    return (chr, pos)

    dictionary = {}

    import os
    print os.getcwd()
    print filepath

    _open = gzip.open if filepath.endswith(".gz") else open
    with _open(filepath) as d:
        print "Loading HapMap CEU SNP list from %s..." % filepath
        for i, line in enumerate(d):
            if i == 0:
                continue
            comps = line.strip().split()
            chromosome = chr_from_snps_file(comps)
            position = pos_from_snps_file(comps)
            rsid = rsid_from_snps_file(comps)
            if rsid_as_key:
                dictionary[rsid] = (chromosome, position)
            else:
                dictionary[(chromosome, position)] = rsid
    return dictionary


def decide_snpid_parsing_function(snpid):

    if re.compile("^rs[0-9]+$").match(snpid):
        return 1, get_rsid
    elif re.compile("^[0-9]{1,2}_[0-9]+$").match(snpid):
        def with_underscore(snpid):
            return get_chr_pos_from_snpid(snpid, sep = "_")
        return 2, with_underscore
    elif re.compile("^[0-9]{1,2}:[0-9]+$").match(snpid):
        def with_colon(snpid):
            return get_chr_pos_from_snpid(snpid, sep = ":")
        return 2, with_colon
    elif re.compile("^[0-9]{1,2}_[0-9]+_[ACGTacgt]+_[ACGTacgt]+").match(snpid):
        return 3, get_info_from_var_id
    else:
        exit("SNP ID format not supported.")


# Estimate the number of lines from the file size and the number of bytes of one line.
def total_n_of_snps_in_eQTL(filename):

    if filename.endswith("gz"):
        # if the file is gzipped we cannot do a reasonable estimation of the number of lines
        return None
    else:
        with open(filename) as eqtl:
            eqtl.readline() # skip header  #TODO: consider headers with several lines
            first_line = eqtl.readline()
            filesize = os.path.getsize(filename)
            est_n_lines = filesize / len(first_line)
            est_n_lines = 100000 * (est_n_lines / 100000)
            return est_n_lines


def show_progress(progress, current_line, total_lines):
    if total_lines is not None:  # if we have an estimation for the number of lines
         k = int(100*current_line/total_lines)
         if progress != k:
             progress = k
             print "{}\% completed ({}/~{})".format(k, current_line, total_lines)
         return progress
    else:
        pass


def get_df_colnames():
    df_colnames = [Constants.CHROMOSOME, \
                   Constants.POSITION, \
                   Constants.REF_ALLELE, \
                   Constants.ALT_ALLELE, \
                   Constants.SNP, \
                   Constants.VAR_ID, \
                   Constants.GENENAME, \
                   Constants.BETA, \
                   Constants.PVALUE]
    return df_colnames


# get the relevant information from the rows of a eQTL file from the GTEx consortium
# the SNP ID for these files is chr_pos_alt_ref_build, e.g. 1_123456789_A_G_b37
# ignore_if_not_in_dictionary == True: ignore rows with SNPs absent from dictionary***
def process_row(comps, pval_thr, indexes, \
                snp_dictionary, \
                parse_snpid, \
                info_from_snp_id, \
                ignore_if_not_in_dictionary = True, \
                ignore_indels = True):

    pvalue = comps[indexes[Constants.PVALUE]]

    if float(pvalue) > pval_thr:
        return None

    snp_id = comps[indexes[Constants.SNP]]

    # extract info from snp ID.
    if info_from_snp_id == 1: # rsid
        info_from_snp = parse_snpid(snp_id)
        snp_id = info_from_snp[Constants.SNP]
        if snp_id not in snp_dictionary:
            if ignore_if_not_in_dictionary:
                return None
        ref_allele = comps[indexes[Constants.REF_ALLELE]]
        alt_allele = comps[indexes[Constants.ALT_ALLELE]]
        var_id = snp_id
        chromosome, position = snp_dictionary[var_id]

    elif info_from_snp_id == 2: # chr, pos
        info_from_snp = parse_snpid(snp_id)
        chromosome, position = info_from_snp[Constants.CHROMOSOME], info_from_snp[Constants.POSITION]
        ref_allele = comps[indexes[Constants.REF_ALLELE]]
        alt_allele = comps[indexes[Constants.ALT_ALLELE]]
        key = (chromosome, position)
        if key not in snp_dictionary:
            if ignore_if_not_in_dictionary:
                return None
        else:
            snp_id = snp_dictionary[key]
        var_id = build_var_id(chromosome, position, ref_allele, alt_allele)

    elif info_from_snp_id == 3: # chr, pos, ref, alt
        info_from_snp = parse_snpid(snp_id)
        chromosome, position, ref_allele, alt_allele = \
            info_from_snp[Constants.CHROMOSOME], info_from_snp[Constants.POSITION], info_from_snp[Constants.REF_ALLELE], info_from_snp[Constants.ALT_ALLELE]
        key = (chromosome, position)
        if key not in snp_dictionary and ignore_if_not_in_dictionary:
            #print "%s is not in dictionary" % (key,)
            return None
        else:
            #print "%s is in dictionary" % (key,)
            snp_id = snp_dictionary[key]
            var_id = build_var_id(chromosome, position, ref_allele, alt_allele)

    # if we are ignoring indels and the variant happens to be one...
    if ignore_indels and (ref_allele not in Constants.alleles or alt_allele not in Constants.alleles):
        return None

    genename = comps[indexes[Constants.GENENAME]]
    beta = comps[indexes[Constants.BETA]]

    row = [chromosome, position, ref_allele, alt_allele, var_id, \
             snp_id, genename, beta, pvalue ]

    return row


def include_gene_(gene, white_list, black_list):
    white_list_condition = (white_list is None or gene in white_list)
    black_list_condition = (black_list is not None and gene in black_list)
    include_gene = white_list_condition and not black_list_condition
    return include_gene


# FIXME: this function is too long
def read_eQTL_file(params):

    found = set()
    gene_column = params.gene_column  # sentinel_column
    data = pandas.DataFrame(columns=get_df_colnames())

    if params.pval_threshold is None:
        params.pval_threshold = 1
    else:
        params.pval_threshold = float(params.pval_threshold)

    # estimate number of SNPs in QTL file, to display progress of model generation
    total_snps = total_n_of_snps_in_eQTL(params.input_file)

    total_for_gene = 0
    passed_for_gene = 0

    _open = gzip.open if ".gz" in params.input_file else open
    with _open(params.input_file) as ifile:

        sentinel = None
        buf = []
        progress = 0
        number_of_genes_passed = 0

        for i, line in enumerate(ifile):
            comps = line.strip().split()

            if i == 0:  # HEADER # TODO: consider the possibility of headers with multiple lines
                # header = ifile.readline().strip().split()
                Utilities.prevalidate_format(params, comps)
                indexes = get_col_indexes(comps, params)
                gene_index = indexes[Constants.GENENAME]
                continue

            # FIRST LINE OF DATA
            # FIXME: put this into a separate function
            if i == 1:
                snpid = comps[indexes[Constants.SNP]]
                snpid_type, parse_snpid = decide_snpid_parsing_function(snpid)
                if snpid_type == 1:
                    snp_dictionary = build_dictionary(params.snp_dictionary, rsid_as_key = True)
                elif snpid_type == 2:
                    snp_dictionary = build_dictionary(params.snp_dictionary, rsid_as_key = False)
                elif snpid_type == 3:
                    snp_dictionary = build_dictionary(params.snp_dictionary, rsid_as_key = False)
                Utilities.validate_input_parameters(params, snpid_type)

            progress = show_progress(progress, i, total_snps)

            row_sentinel = comps[gene_index]

            if sentinel is None:
                sentinel = row_sentinel
                include_gene = include_gene_(sentinel, params.gene_white_list, params.gene_black_list)

            # if the current row has still the same gene as the previous
            if row_sentinel == sentinel:
                if include_gene:
                    total_for_gene += 1
                    new_row = process_row(comps, params.pval_threshold, indexes, \
                                          snp_dictionary, parse_snpid, \
                                          snpid_type, params.ignore_if_not_in_dictionary, \
                                          params.ignore_indels)

                    if new_row is not None:
                        passed_for_gene += 1
                        buf.append(new_row)
                continue

            else: # if the gene just changed, we have to decide whether we have to include it or not
                include_gene = include_gene_(sentinel, params.gene_white_list, params.gene_black_list)
                if include_gene:
                    number_of_genes_passed += 1

            # if we're here, this is the last line for the current gene (because the gene just changed),
            # therefore we write the data into a pandas dataframe
            logging.info("Done with gene {} ({}/{} SNPs passed)".format(sentinel, passed_for_gene, total_for_gene))
            data = to_dataframe(buf, get_df_colnames(), to_numeric="ignore")
            yield data

            total_for_gene = 0
            passed_for_gene = 0

            if params.gene_white_list is not None:
                found.add(sentinel)
                if len(found) == len(params.gene_white_list):
                    logging.info("Found everything in the whitelist")
                    return

            new_row = process_row(comps, params.pval_threshold, indexes, \
                                  snp_dictionary, parse_snpid, \
                                  snpid_type, params.ignore_if_not_in_dictionary, \
                                  params.ignore_indels)

            buf = [new_row] if new_row else []
            sentinel = comps[gene_index]

        # If we reach the end of the file and buf still contains something...
        if len(buf):
            data = to_dataframe(buf, get_df_colnames(), to_numeric="ignore")
            yield data
