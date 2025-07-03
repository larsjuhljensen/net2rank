import h5py
import numpy as np
from typing import List, Tuple, Iterable
import pandas as pd
import os

# load h5 embeddings
class H5Loader:
    """Class to load protein embeddings from an h5 file."""
    def __init__(self, h5_file:str):
        self.h5_file = h5_file
        self.proteins, self.embeddings = self.load_embeddings(h5_file)
        self.prot2ind = {protein: i for i, protein in enumerate(self.proteins)}
    
    def load_embeddings(self,embeddings_file:str) -> Tuple[List[str], np.ndarray]:
        """Load embeddings from an h5 file."""
        if not embeddings_file.endswith('.h5'):
            raise ValueError("The provided file is not an h5 file.")
        with h5py.File(embeddings_file, 'r') as f:
            proteins = list(f.keys())
            embeddings = np.array([f[protein][:] for protein in proteins])
            
        return proteins, embeddings
    
    def get_embeddings(self,proteins):
        """Get embeddings for a list of proteins."""
        # check if any proteins are in the loaded embeddings
        if len(self.proteins) == 0 or len(self.embeddings) == 0:
            raise ValueError("No embeddings loaded. Please check the h5 file.")
        for protein in proteins:
            if protein not in self.prot2ind:
                raise ValueError(f"Protein {protein} not found in the embeddings.")
        
        indices = [self.prot2ind[protein] for protein in proteins if protein in self.prot2ind]
        return self.embeddings[indices] if indices else np.array([])
    

def pos_neg_threshold(pos_size,  omics_file):
    """create a positive and negative set from your omics input file.
    Positives will be the top x in our omics data file, which is sorted so that most important genes at top of list.
    The top x is determined by checking a benchmark plot for the disease and it's gold standard.
    Negatives are anything that are not in the top X * 3, so if top 1000, you do not sample negatives from top 3000.
    Negative and positive lists will be saved in folder where code is run, for safety.
    """

    pos_size = pos_size
    pos_omics = omics_file[:pos_size]

    neg_start_idx = (pos_size*3)+1
    # neg_size = pos_size*5
    neg_omics = omics_file.iloc[neg_start_idx:]

    proteins = pos_omics.iloc[:, 0].tolist() + neg_omics.iloc[:, 0].tolist()
    classes = [1] * len(pos_omics) + [0] * len(neg_omics)
    ###pos_neg.to_csv("pos_neg_set.tsv", index=False, sep="\t")
    return pd.DataFrame({
        'protein': proteins,
        'class': classes
    })


def process_train_file(train_file: str, file_type: str, 
                       protein_space: Iterable[str],
                       pos_size: int = 1000,
                       ):
    
    if file_type == 'pvalue':
        df_train = pd.read_csv(train_file, sep='\t')
        df_train = df_train.sort_values(by=df_train.columns[1], ascending=True)
        df_train = pos_neg_threshold(pos_size, df_train)
        # filter out proteins not in the protein space
        df_train = df_train[df_train['protein'].isin(protein_space)]
        df_train.columns = ['protein', 'class']
        return df_train
    
    elif file_type == 'list':
        # do negative sampling
        pos = open(train_file, 'r').read().splitlines()
        pos = [p.strip() for p in pos if p.strip() in protein_space]
        neg = list(set(protein_space) - set(pos))
        return pd.DataFrame({
            'protein': pos + neg,
            'class': [1] * len(pos) + [0] * len(neg)
        })
    
    elif file_type == 'label':
        # just load the file
        df = pd.read_csv(train_file, sep='\t', header=None)
        df.columns = ['protein', 'class']
        return df
    
    
def save_predictions(y_test, y_pred, test_set, save_dir):
            
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    
    # save the predictions
    predictions_df = pd.DataFrame({
        'protein': test_set.iloc[:, 0],
        'true_label': y_test,
        'predicted_score': y_pred
    })
    
    # sort the predictions by predicted score
    predictions_df = predictions_df.sort_values(by='predicted_score', ascending=False)
    
    return predictions_df
