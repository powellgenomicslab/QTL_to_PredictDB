def main(input_file, snp_dictionary, snp_column, gene_column, beta_column, pval_threshold, pvalue_column):

  from subprocess import check_output #call
  import tempfile
  import os

  command = """python Main.py --input_file {input_file} \
           --snp_dictionary {snp_dictionary} \
           --snp_column {snp_column} \
           --gene_column {gene_column} \
           --beta_column {beta_column} \
           --pval_threshold {pval_threshold} \
           --pvalue_column {pvalue_column} \
           --output_file {output_filepath}"""

  # print snp_dictionary
  
  dirpath = tempfile.mkdtemp()
  input_file_dir, input_file_basename = os.path.split(input_file)
  output_filename = "".join(input_file_basename.split(".")[:-2]) if input_file_basename.endswith("gz") else "".join(input_file_basename.split(".")[:-1])
  output_filename += ".db"
  #output_filepath = os.path.join(dirpath, output_basename)

  command = command.format(input_file = input_file, \
                           snp_dictionary = snp_dictionary, \
                           snp_column = snp_column, \
                           gene_column = gene_column, \
                           beta_column = beta_column, \
                           pval_threshold = pval_threshold, \
                           pvalue_column = pvalue_column, \
                           output_filepath = output_filename)

  print command

  check_call(command.split()) 

  return output_filename
