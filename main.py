import fire 
import pandas as pd
import random  
import numpy as np
from net2rank.utils import H5Loader, process_train_file, save_predictions

from sklearn.model_selection import KFold
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.linear_model import LogisticRegression

import matplotlib.pyplot as plt
import os

class Net2rank:
    """
    A machine learning framework for protein-disease association prediction using network embeddings.
    
    This class provides methods for cross-validation and train-test evaluation of protein embeddings
    for disease prediction tasks using logistic regression models.
    
    Usage:
        python main.py cross_validation --help
        python main.py train_test --help
    """
    
    def cross_validation(self, train_file:str, file_type:str, 
                         embedding_file: str = 'data/9606.node2vec64.h5',
                         pos_size: int = 1000,
                         k: int = 5, 
                         roc_plot: bool = True,
                         save_dir: str = 'results/cv_results',
                         random_state: int = 42):
        """
        Perform k-fold cross-validation on protein-disease association data.
        
        This method splits the training data into k folds, trains a logistic regression
        model on k-1 folds, and evaluates on the remaining fold. This process is repeated
        k times to provide robust performance estimates.
        
        Args:
            train_file (str): Path to the training data file containing protein-disease associations.
            file_type (str): Type of the training file format (e.g., 'pvalue', 'list', 'label').
            embedding_file (str): Path to the H5 file containing protein embeddings.
            pos_size (int, optional): Maximum number of positive samples to use. Defaults to 1000.
            k (int, optional): Number of folds for cross-validation. Defaults to 5.
            roc_plot (bool, optional): Whether to generate and save ROC curve plots. Defaults to True.
            save_dir (str, optional): Directory to save results and plots. Defaults to 'results/cv_results'.
            random_state (int, optional): Random seed for reproducibility. Defaults to 42.
            
        Returns:
            None: Results are saved to files and printed to console.
            
        Raises:
            ValueError: If roc_plot is True but save_dir is None.
            
        Note:
            - The method maintains class balance across folds
            - Results include prediction scores, true labels, and AUC values for each fold
            - Combined results from all folds are saved as a TSV file
            - ROC curves are plotted for each fold if roc_plot is enabled
       """
        
        
        random.seed(random_state)
        np.random.seed(random_state)
        
        if roc_plot and save_dir is None:
            raise ValueError("If roc_plot is True, save_dir must be specified to save the ROC curve.")
        
        human_embeddings = H5Loader(embedding_file)
        protein_space = human_embeddings.proteins
        df_train = process_train_file(train_file, file_type, protein_space,pos_size)
        disease_name = os.path.basename(train_file).split('.')[0]
        
        # make sure that we keep the balance of positive and negative samples
        if save_dir is not None:
            save_dir = os.path.join(save_dir, disease_name)
            save_results = []
            if not os.path.exists(save_dir):
                os.makedirs(save_dir)
        for i, (train_index, test_index) in enumerate(KFold(n_splits=k, shuffle=True, random_state=random_state).split(df_train)):
            train_set = df_train.iloc[train_index]
            test_set = df_train.iloc[test_index]
            
            X_train = human_embeddings.get_embeddings(train_set.iloc[:,0].tolist())
            y_train = train_set.iloc[:,1].values
            
            X_test = human_embeddings.get_embeddings(test_set.iloc[:,0].tolist())
            y_test = test_set.iloc[:,1].values
            
            model = LogisticRegression(max_iter=1000, random_state=random_state)

            model.fit(X_train, y_train)
            
            y_pred = model.predict_proba(X_test)[:, 1]
            auc = roc_auc_score(y_test, y_pred)
            print(f"Fold {i+1}/{k}, AUC: {auc:.4f}")
            
            if save_dir is not None:
                fold_result = save_predictions(y_test, y_pred, test_set, save_dir)
                fold_result['fold'] = i + 1
                save_results.append(fold_result)
        
        if save_dir is not None:
            # save the combined results of all folds
            combined_results = pd.concat(save_results)
            combined_results.to_csv(os.path.join(save_dir, f'{disease_name}_cv_prediction_results.tsv'), sep='\t', index=False)
            

            if roc_plot:  
                for fold_index,fold_result in enumerate(save_results):
                    y_test = fold_result['true_label']
                    y_pred = fold_result['predicted_score']
                    fpr, tpr, _ = roc_curve(y_test, y_pred)
                    plt.plot(fpr, tpr, label=f'Fold {fold_index+1}: AUC = {roc_auc_score(y_test, y_pred):.2f}')
                plt.plot([0, 1], [0, 1], 'k--')
                plt.xlabel('False Positive Rate')
                plt.ylabel('True Positive Rate')
                plt.title(f'ROC Curve for {disease_name}')
                plt.legend()
                plt.savefig(os.path.join(save_dir, f'{disease_name}_roc_curve.png'),dpi=300)
                plt.close()
        print(f"Cross-validation completed for {disease_name}. Results saved to {save_dir if save_dir else 'current directory'}.")
        return None
            
    def train_test(self, train_file: str, 
                   file_type: str,
                   test_file: str,
                   embedding_file: str='data/9606.node2vec64.h5',
                   pos_size: int = 1000,
                   roc_plot: bool = True,
                   save_dir: str = 'results/test_results',
                   random_state: int = 42,
                   ):
        """
        Train a model on training data and evaluate on separate test data.
        
        This method trains a logistic regression model on the provided training data
        and evaluates its performance on a separate test dataset. This is useful for
        final model evaluation or when you have predefined train/test splits.
        
        Args:
            train_file (str): Path to the training data file containing protein-disease associations.
            file_type (str): Type of the training file format (e.g., 'pvalue', 'list', 'label').
            embedding_file (str): Path to the H5 file containing protein embeddings.
            test_file (str): Path to the test data file for evaluation.
            pos_size (int, optional): Maximum number of positive samples to use. Defaults to 1000.
            roc_plot (bool, optional): Whether to generate and save ROC curve plot. Defaults to True.
            save_dir (str, optional): Directory to save results and plots. Defaults to 'results/test_results'.
            random_state (int, optional): Random seed for reproducibility. Defaults to 42.
            
        Returns:
            None: Results are saved to files and printed to console.
            
        Raises:
            ValueError: If roc_plot is True but save_dir is None.
            
        Note:
            - The test file is processed with 'label' format regardless of train file format
            - Results include prediction scores, true labels, and AUC value
            - ROC curve is plotted and saved if roc_plot is enabled
            - All results are saved in a subdirectory named after the disease
        """
        
        if roc_plot and save_dir is None:
            raise ValueError("If roc_plot is True, save_dir must be specified to save the ROC curve.")
        
        human_embeddings = H5Loader(embedding_file)
        protein_space = human_embeddings.proteins
        df_train = process_train_file(train_file, file_type, protein_space,pos_size,)
        df_test = process_train_file(test_file, 'label', protein_space, pos_size,)
        
        disease_name = os.path.basename(train_file).split('.')[0]
        
        if save_dir is not None:
            save_dir = os.path.join(save_dir, disease_name)
            if not os.path.exists(save_dir):
                os.makedirs(save_dir)
                
        X_train = human_embeddings.get_embeddings(df_train.iloc[:,0].tolist())
        y_train = df_train.iloc[:,1].values
        
        X_test = human_embeddings.get_embeddings(df_test.iloc[:,0].tolist())
        y_test = df_test.iloc[:,1].values
        
        model = LogisticRegression(max_iter=1000, random_state=random_state)
        model.fit(X_train, y_train)
        
        y_pred = model.predict_proba(X_test)[:, 1]
        auc = roc_auc_score(y_test, y_pred)
        print(f"Test AUC: {auc:.4f}")
        
        if save_dir is not None:
            fold_result = save_predictions(y_test, y_pred, df_test, save_dir)
            fold_result.to_csv(os.path.join(save_dir, 
                                            f'{disease_name}_test_prediction_results.tsv'), 
                               sep='\t', index=False)
        
            if roc_plot:  
                fpr, tpr, _ = roc_curve(y_test, y_pred)
                plt.plot(fpr, tpr, label=f'AUC = {auc:.2f}')
                plt.plot([0, 1], [0, 1], 'k--')
                plt.xlabel('False Positive Rate')
                plt.ylabel('True Positive Rate')
                plt.title(f'ROC Curve for {disease_name}')
                plt.legend()
                plt.savefig(os.path.join(save_dir, f'{disease_name}_roc_curve.png'), dpi=300)
                plt.close()
                
        print(f"Test completed for {disease_name}. Results saved to {save_dir if save_dir else 'current directory'}.")
        return None


if __name__ == "__main__":
    fire.Fire(Net2rank)
