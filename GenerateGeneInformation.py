#/usr/bin/env python

import pickle
import gzip

GENE_ID = "gene_id"
GENE_NAME = "gene_name"
GENE_TYPE = "gene_type"
CHROMOSOME = "chromosome"
FEATURE = "feature"
START = "start"
END = "end"
STRAND = "strand"

CHROMOSOME_idx = 0
FEATURE_idx = 2
START_idx = 3
END_idx = 4
STRAND_idx = 6
KEY_VALUE_PAIRS = 8


def extract_gene_info(dictionary, key_value_pairs):

  for x in key_value_pairs:
    if GENE_ID in x:
      idx = x.index(GENE_ID)
      gene_id = x[idx+1]
      gene_id = gene_id.replace('"', '')
      gene_id = gene_id[0:15]
    if GENE_NAME in x:
      idx = x.index(GENE_NAME)
      gene_name = x[idx+1]
      gene_name = gene_name.replace('"', '')
    if GENE_TYPE in x:
      idx = x.index(GENE_TYPE)
      gene_type = x[idx+1]
      gene_type = gene_type.replace('"', '')

  dictionary[GENE_ID] = gene_id
  dictionary[GENE_NAME] = gene_name
  dictionary[GENE_TYPE] = gene_type
  return dictionary
 

def extract_geneID_and_genename(dictionary, key_value_pairs):
  entry = None
  for x in key_value_pairs:
    if GENE_ID in x:
      idx = x.index(GENE_ID)
      key = x[idx+1]
      key = key.replace('"', '')
      key = key[0:15]
    if GENE_NAME in x:
      idx = x.index(GENE_NAME)
      value = x[idx+1]
      value = value.replace('"', '')
  if key is not None and value is not None:
    dictionary[key] = value
  

def generate_dictionary(gencode_file, key):

  d = {}
  _open = gzip.open if gencode_file.endswith(".gz") else open

  with _open(gencode_file) as gencode:
    for i, line in enumerate(gencode):
      if "##" in line:
        continue
      comps = line.strip().split("\t")
      key_value_pairs = comps[KEY_VALUE_PAIRS].split(";")
      key_value_pairs = [x.split(" ") for x in key_value_pairs]
      
      feature = comps[FEATURE_idx]
      if feature == "gene":
        dd = {}
        dd = extract_gene_info(dd, key_value_pairs)
        dd[CHROMOSOME] = comps[CHROMOSOME_idx]
        dd[START] = comps[START_idx]
        dd[END] = comps[END_idx]
        dd[STRAND] = comps[STRAND_idx]

        gene_id = dd[GENE_ID]
        d[gene_id] = dd

        if i == 1000:
          print "Processed %s lines" % str(i)
          break

  return d


def save_pkl(obj, pkl_file):
  with open(pkl_file, "wb") as ff: 
    pickle.dump(obj, ff) 


def load_pkl(pkl_file):
  with open(pkl_file, "rb") as ff: 
    d = pickle.load(ff)
  return d


def run(args):

  d = generate_dictionary(args.gencode_file, args.dictionary_key)
  
  if args.save_as_binary:
    if not args.output_file.endswith(".pkl"):
      args.output_file = args.output_file + ".pkl"
    save_pkl(d, args.output_file) 
  else:
    with open(args.output_file, "w") as of: 
      header = [GENE_ID, GENE_NAME, GENE_TYPE, CHROMOSOME, START, END, STRAND]
      header = "%s\n" % "\t".join(header)
      of.write(header)
      for key in d:
        dd = d[key]
        of.write("%s\n" % "\t".join([dd[GENE_ID], dd[GENE_NAME], dd[GENE_TYPE], dd[CHROMOSOME], dd[START], dd[END], dd[STRAND]]))
    
        
if __name__ == "__main__":
  
  import argparse

  parser = argparse.ArgumentParser()

  parser.add_argument("--gencode_file")
  parser.add_argument("--output_file")
  parser.add_argument("--dictionary_key", default="gene_id")
  parser.add_argument("--save_as_binary", action="store_true")
 
  args = parser.parse_args()

  run(args)
