Code for generating prediction models for molecular traits (gene expression, methylation rate, intron excision, etc.) as SQLite databases, in the format required by both PrediXcan and S-PrediXcan, starting from the corresponding QTL summary statistics.

------------
Example on how to run from Rodrigo's email

Let's say that the eQTL summary statistics files looks like this:

> Beta    SE      Pvalue  gene_id CHROM:POS:NON_EFF:EFF   EAF\
> 0.426   0.16    0.00858 ENSG00000006555 1:54965556:G:A  0.148\
> 0.341   0.148   0.0225  ENSG00000006555 1:54968447:T:A  0.188\
> 0.404   0.166   0.016   ENSG00000006555 1:54969716:C:A  0.141\
> 0.373   0.174   0.034   ENSG00000006555 1:54975683:G:A  0.129

and there is a a set of such files (one for each chromosome) in a folder called sum_stats.

Notice that you have no separate rsid/EffAllele/NonEffAllele columns, but a single column condensing chromosome, position and alleles.
In this case, you will need to provide a SNP annotation file that includes chromosome, position and rsID, like this one which is based on hg19.
(If you already have rsid's it's simpler)

The command to run the script for this case would be teh following (look at how the regular expression for the SNP ID format is defined):

> python Main.py --input_folder sum_stats \
>               --snp_column CHROM:POS:NON_EFF:EFF \
>               --snpid_format "(?P<chromosome>.*):(?P<position>.*):(?P<non_effect_allele>.*):(?P<effect_allele>.*)" \
>               --gene_column gene_id \
>               --beta_column Beta \
>               --pvalue_column Pvalue \
>               --pval_threshold 0.01 \
>               --snp_annot_file gtex_v7_hapmapceu_dbsnp150_snp_annot.txt.gz \
>               --snp_annot_snpid_column rsid \
>              --snp_annot_chr_column chromosome \
>               --snp_annot_position_column pos \
>               --output_file output.db

In this command, the snp_annot_*_column arguments correspond to column names in the SNP annotation file. The rest of the *_column arguments are column names of the input eQTL file. I'm using a p-value threshold of 1%. If you don't want to use a p-value threshold to filter the SNPs, you can just remove the --pval_threshold argument.
You can check some more options by running python Main.py --help.
Let me know if I can help with anything else.
