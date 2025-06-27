import pandas as pd

# function to generate positive and negative set fro training ML model
def pos_neg(pos_size,  omics_file):
    """create a positive and negative set from your omics input file.
    Positives will be the top x in our omics data file, which is sorted so that most important genes at top of list.
    The top x is determined by checking a benchmark plot for the disease and it's gold standard.
    Negatives are anything that are not in the top X * 3, so if top 1000, you do not sample negatives from top 3000.
    Negative and positive lists will be saved in folder where code is run, for safety.
    """

    pos_size = pos_size
    pos_omics = omics_file[:pos_size]

    neg_start_idx = (pos_size*3)+1
    neg_size = pos_size*5
    neg_omics_sampler = omics_file.iloc[neg_start_idx:]

    neg_omics = neg_omics_sampler.sample(n=neg_size, random_state=42)
    # Assign classes to the pos and neg data
    pos_omics['class'] = 1
    neg_omics['class'] = 0
    
    pos_neg = pd.concat([pos_omics, neg_omics])
    ###pos_neg.to_csv("pos_neg_set.tsv", index=False, sep="\t")
    return pos_neg


