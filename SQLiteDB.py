import sqlite3

from Constants import SNP
from Constants import REF_ALLELE
from Constants import ALT_ALLELE
from Constants import CHROMOSOME
from Constants import GENENAME
from Constants import POSITION
from Constants import BETA
from Constants import SE
from Constants import PVALUE

from Constants import SNP_idx
from Constants import REF_ALLELE_idx
from Constants import ALT_ALLELE_idx
from Constants import CHROMOSOME_idx
from Constants import GENENAME_idx
from Constants import POSITION_idx
from Constants import BETA_idx
from Constants import SE_idx
from Constants import PVALUE_idx


class DB:
    def __init__(self, db_file):
        ''' Opens connection to database '''
        self.connection = sqlite3.connect(db_file)
        self.c = self.connection.cursor()
        self("CREATE TABLE weights (rsid TEXT, gene TEXT, weight DOUBLE, varID TEXT, ref_allele CHARACTER, eff_allele CHARACTER, pval DOUBLE, N INTEGER, cis INTEGER)")
        self("CREATE TABLE extra (gene TEXT, genename TEXT, [pred.perf.R2] DOUBLE, [pred.perf.pval] DOUBLE, [pred.perf.qval] DOUBLE)")

    def __call__(self, sql, args=None):
        c = self.connection.cursor()
        if args:  c.execute(sql, args) # print sql # DEBUGGING
        else:  c.execute(sql)  # print sql # DEBUGGING

    def insert_row(self, row):
        self("INSERT INTO weights VALUES(\"%s\", \"%s\", \"%s\", \"%s\", \"%s\", \"%s\", NULL, NULL)" % \
              (row[SNP_idx], row[GENENAME_idx], row[BETA_idx], row[REF_ALLELE_idx], row[ALT_ALLELE_idx], row[PVALUE_idx]))
        """alt allele is the dosage/effect allele in GTEx data"""

    def insert_row_extra(self, row):
        self("INSERT INTO extra VALUES(\"%s\", \"%s\", NULL, NULL, NULL)" % (row[GENENAME_idx], row[GENENAME_idx]))

    def create_indexes(self):
        self("CREATE INDEX weights_rsid ON weights (rsid)")
        self("CREATE INDEX weights_gene ON weights (gene)")
        # OPTIMIZE: what about reversing the order in the double index?
        # I think it should improve the performance when reading. But I'm not sure it'll be an appreciable gain.
        self("CREATE INDEX weights_rsid_gene ON weights (rsid, gene)")

    def drop_indexes(self):
        self("DROP weights_rsid")
        self("DROP weights_gene")
        self("DROP weights_rsid_gene")

    def close(self):
        self.connection.commit()
