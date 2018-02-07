import Constants

# FIXME: check default value for n
def top_n_eQTL(df, n):
    df = df.nsmallest(int(n), Constants.PVALUE)
    return df


def generate_weights(df, args):
    if args.method == "naive":
        if args.max_n_of_snps is None:
            return df
        else:
            return top_n_eQTL(df, args.max_n_of_snps)
    elif args.method == "BSLMM":
        return generate_weights_top_eQTL(df, cov)
    elif args.method == "LDPred":
        return generate_weights_LDPred(df, cov)
    elif args.method == "LASSO":
        return generate_weights_LASSO(df, cov)
    elif args.method == "BLUP":
        return generate_weights_BLUP(df, cov)
    return
