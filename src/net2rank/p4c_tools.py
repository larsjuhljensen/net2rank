import py4cytoscape as p4c
import pandas as pd
from sklearn.metrics import roc_curve
from typing import List, Dict, Iterable

import json
import requests
import math
from net2rank.utils import process_train_file, H5Loader
import os

def make_network(proteins:Iterable[str],
                 network_name:str,
                 cutoff:float=0.7,
                 species:str='Homo sapiens',
                 networktype='full STRING network') -> str:            
    protein_string = ','.join(proteins)
    command = (f'string protein query '
            f'query="{protein_string}" '
            f'species="{species}" '
            f'cutoff="{cutoff}" '
            f'networkType="{networktype}" '
            f'limit="0" '
            f'newNetName="{network_name}" ')
    suid = p4c.commands_post(command)
    return suid

def network_cluster(inflation_parameter=4,showUI=False):
    cluster_cmd = ('cluster mcl '
            'network=current '
            f'inflation_parameter={inflation_parameter} '
            f'showUI={showUI}'
            )
    cluster_result = p4c.commands_post(cluster_cmd) 
    return cluster_result

def function_enrichment(proteins):
    string_api_url = "https://version-12-0.string-db.org/api"
    output_format = "json"
    method = "enrichment"

    request_url = "/".join([string_api_url, output_format, method])

    params = {

        "identifiers" : "%0d".join(proteins), # your protein
        "species" : 9606, # NCBI/STRING taxon identifier 
        "caller_identity" : "net2rank" # your app name
    }

    response = requests.post(request_url, data=params)

    data = json.loads(response.text)
    
    return data


def parse_cluster_result(cluster_result):
    
    clusters = cluster_result['clusters']
    
    node_list = []
    for cluster in clusters:
        nodes = [node['name'] for node in cluster['nodes']]
        
        if len(nodes) < 10:
            break
        node_list.append(nodes)
    return node_list


def function_enrichment_from_clusters(cluster_result:List) -> pd.DataFrame:
    """ Perform function enrichment analysis for each cluster in the cluster result.
    Parameters
    ----------
    cluster_result : List
        List of cluster results, where each cluster contains a list of proteins.
        From the output of `network_cluster()` function.
    Returns
    -------
    pd.DataFrame
        DataFrame containing enrichment results for each cluster.
    """
    # enrichment analysis per cluster
    result = list()
    cluster_proteins = list()
    for idx,proteins in enumerate(parse_cluster_result(cluster_result)):
        enrichment = function_enrichment(proteins=proteins)
            # select cat
        use_cats = {'KEGG','Process','RCTM','WikiPathways','Process'}

        for row in enrichment:
            if float(row["fdr"]) > 0.05:
                continue
            
            # if category is not None and row["category"] != category:
            #     continue
            
            if row["category"] not in use_cats:
                continue

            result.append({
                'term': row['term'],
                'description': row["description"],
                'fdr': float(row["fdr"]),
                'num_genes': row['number_of_genes'],
                'num_genes_background': row['number_of_genes_in_background'],
                'category': row["category"],
                'cluter_number':idx
            }) 
        cluster_proteins.append(proteins)
    return pd.DataFrame(result),cluster_proteins


def parse_enrichment_results(enrichment_results:List,category=None) -> pd.DataFrame:
    result = list()
    
    for row in enrichment_results:
        if float(row["fdr"]) > 0.05:
            continue
        
        if category is not None and row["category"] != category:
            continue

        result.append({
            'term': row['term'],
            'description': row["description"],
            'fdr': float(row["fdr"]),
            'num_genes': row['number_of_genes'],
            'num_genes_background': row['number_of_genes_in_background'],
            
            'category': row["category"],
        })
    df = pd.DataFrame(result)
    df = df.sort_values(by='fdr')
    
    return df

def select_proteins_by_enrichment(df:pd.DataFrame, term:str, category:str='Process') -> List[str]:
    """
    Select proteins by enrichment term and category.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing enrichment results.
    term : str
        The term to filter by.
    category : str
        The category to filter by (default is 'Process').
        
    Returns
    -------
    List[str]
        List of proteins associated with the specified term and category.
    """
    filtered_df = df[(df['term'] == term) & (df['category'] == category)]
    
    if filtered_df.empty:
        return []
    
    return filtered_df['num_genes'].tolist()

def load_protein_data(predictions_file:str,
                      train_file:str,
                      train_file_type:str,
                      pos_size:int=1000,
                      human_network_emb:str='../data/9606.node2vec64.h5',
                      threshold:float=0.05,
                      ) -> pd.DataFrame:

    df_predictions = pd.read_csv(predictions_file, sep='\t')
    fpr, tpr, thresholds = roc_curve(df_predictions['true_label'], df_predictions['predicted_score'])
    threshold = thresholds[fpr <= 0.05][-1]
    df_predictions = df_predictions[df_predictions['predicted_score'] >= threshold]
    # protein_list = df_predictions[df_predictions['predicted_score'] >= threshold]['protein'].unique()

    protein_space = H5Loader(human_network_emb).proteins

    # train_file = '../data/train/atopic_dermatitis.integrated.tsv'
    # train_file_type = 'pvalue'
    # pos_size = 1000
    df_train = process_train_file(train_file, train_file_type, protein_space,pos_size)

    df_protein_category = list()

    for idx, row in df_predictions.iterrows():
        
        protein = row['protein']
        text_mining = True if row['true_label'] == 1 else False
        
        protein_in_training = df_train[df_train['protein'] == protein] ## all the proteins in the training set
        # Check if the protein is in the training set and has a class of 1
        if not protein_in_training.empty and protein_in_training['class'].values[0] == 1:
            training = True
        else:
            training = False
        
        df_protein_category.append({
            'protein': protein,
            'text_mining': text_mining,
            'training': training
        })
    df_protein_category = pd.DataFrame(df_protein_category)

    categories = list()
    for idx, row in df_protein_category.iterrows():
        if row['text_mining'] and row['training']:
            categories.append('both')
        elif row['text_mining'] and not row['training']:
            categories.append('text_mining')
        elif not row['text_mining'] and row['training']:
            categories.append('omics')
        else:
            categories.append('novel')
    df_protein_category['category'] = categories
    df_protein_category = df_protein_category.merge(df_predictions[['protein', 'predicted_score', 'true_label']], on='protein', how='left')
    
    return df_protein_category

def make_network_plot(df_protein_category:pd.DataFrame, enrichment_results:List[Dict],selected_term:str,
                      disease_name:str,cutoff:float=0.7,networktype:str='full STRING network') -> str:

    # start making the module network plot
    df_enrichment = parse_enrichment_results(enrichment_results)
    process_name = df_enrichment[df_enrichment['term'] == selected_term]['description'].values[0]

    for row in enrichment_results:
        if row['category'] == 'Process' and row['term'] == selected_term:
            selected_proteins = row['inputGenes']
            break
        
    suid = make_network(
        proteins=selected_proteins,
        network_name=f'{disease_name}_module_network_{process_name}',
        cutoff=cutoff,
        networktype=networktype)



def make_network_plot_with_pie_charts(df_protein_category:pd.DataFrame, 
                                      selected_proteins,suid, map_node_size=False,
                                      ukb_proteins=None) -> str:

    # Get the proteins that are in the 'both' category for pie chart logic
    both_proteins = df_protein_category[df_protein_category['category'] == 'both']['protein'].values

    if map_node_size:
        if ukb_proteins is None:
            raise ValueError("ukb_proteins must be provided when map_node_size is True")
        ukb_proteins = open(ukb_proteins).read().splitlines()
    
        df_protein_category['ukb_measured'] = 'not_measured'
        df_protein_category.loc[df_protein_category['protein'].isin(ukb_proteins), 'ukb_measured'] = 'measured'
    
    pie_data = []
    for node in selected_proteins:
        if node in both_proteins:
            # For novel candidates, create pie chart data (e.g., 50% each category)
            pie_data.append({'name': node, 'value1': 50, 'value2': 50})
        else:
            # For other nodes, set to 0 so they don't show pie charts
            pie_data.append({'name': node, 'value1': 0, 'value2': 0})

    pie_df = pd.DataFrame(pie_data)

    p4c.set_visual_style('Revelen',network=suid)

    # Pass a copy of the DataFrame to avoid SettingWithCopyWarning
    p4c.load_table_data(df_protein_category.copy(), data_key_column='protein', table='node')
    
    # Set node colors
    p4c.set_node_color_mapping(
        table_column='category',
        table_column_values=['text_mining', 'training', 'novel'],
        colors=['#fa6600', '#27a59b', '#d8d8d8'],
        mapping_type='d',
        style_name='Revelen',
        network=suid
    )
    
    if map_node_size:
        # Set node sizes - small for 'training' and 'both', big for 'text_mining'
        p4c.set_node_size_mapping(
            table_column='ukb_measured',
            table_column_values=['not_measured', 'measured'],
            sizes=[30, 30*math.sqrt(2)],  # Small for not_measured, big for measured
            mapping_type='discrete',
            style_name='Revelen',
            network=suid
        )


    # Set up the pie chart with start angle from top (270 degrees or -90 degrees)
    p4c.load_table_data(pie_df.copy(), data_key_column='name', table='node')
    p4c.set_node_custom_pie_chart(
        columns=['value1', 'value2'], 
        colors=['#fa6600', '#27a59b'],  # Orange and teal
        slot=1, 
        style_name='Revelen',
        start_angle=270  # This sets the start position to the top (12 o'clock)
    )
    return suid


def pick_cluster_proteins(clusters, idx):
    """
    Select proteins from a specific cluster index.
    """
    return [node_dict['name'] for node_dict in clusters['clusters'][idx]['nodes']]

def cluster_enrichment_pipeline(disease_name,
                                predictions,
                                cutoff=0.7,networktype='full STRING network',topk=None):
    """
    Main pipeline to load predictions, create network, cluster, and perform enrichment analysis.
    Args:
        disease_name (str): Name of the disease to load predictions for.
        cutoff (float): Cutoff value for network creation.
        networktype (str): Type of network to create. Defaults to 'full STRING network'; or 'physical subnetwork'.
    """
    # df_cat = load_predictions(disease_name,file_pattern)
    df_cat = pd.read_csv(predictions, sep='\t')
    if topk is not None:
        df_cat = df_cat.sort_values(by='predicted_score', ascending=False).head(topk)
    
    suid = make_network(
        proteins=df_cat['protein'].tolist(),
        network_name=f'{disease_name}_main_network',
        cutoff=cutoff,
        networktype=networktype)
    
    cluster_result = network_cluster()
    
    df_clusters_enrichment, cluster_proteins = function_enrichment_from_clusters(cluster_result)
    
    return df_cat,df_clusters_enrichment, cluster_proteins

def plot_cluster_network(cluster_proteins, idx, df,disease_name):
    """
    Plot the network for a specific cluster index with pie charts.
    """
    suid = make_network(
        proteins=cluster_proteins[idx],
        network_name=f'{disease_name}_cluster_{idx}_network',
        cutoff=0.7)
    
    suid = make_network_plot_with_pie_charts(
        df_protein_category=df,
        selected_proteins=cluster_proteins[idx],
        suid=suid,)

    p4c.commands.commands_run(f'layout apply preferred networkSelected={suid}')
    p4c.notebook_export_show_image()
    
    return suid

def save_clusters_and_enrichment(disease_name,
                                 save_to_dir,
                                 cutoff=0.7,
                                 networktype='full STRING network',topk=None):
    if not os.path.exists(save_to_dir):
        os.makedirs(save_to_dir)
    
    df_cat, df_clusters_enrichment,cluster_proteins = cluster_enrichment_pipeline(disease_name,
                                                                                  save_to_dir,
                                                                                  cutoff,
                                                                                  networktype,topk) 
    
    # save the cluster enrichment results and cluster proteins
    networktype = networktype.replace(' ', '_').lower()
    file_name = f'{disease_name}_clusters_enrichment_{cutoff}_networktype_{networktype}.tsv'
    df_clusters_enrichment.to_csv(os.path.join(save_to_dir, file_name), sep='\t', index=False)
    
    # save the cluster proteins
    file_name = f'{disease_name}_clusters_proteins_{cutoff}_networktype_{networktype}.tsv'
    with open(os.path.join(save_to_dir, file_name), 'w') as f:
        f.write('protein\tcluster_idx\n')
        for idx, cluster in enumerate(cluster_proteins):
            for protein in cluster:
                f.write(f'{protein}\t{idx}\n')
    print(f'Saved {disease_name} cluster enrichment and proteins to {save_to_dir}')      
    return None

def load_records(disease_name, 
                 predictions,
                 file_pattern,cut_off=0.7,networktype='full STRING network',):
    df_cat = pd.read_csv(predictions, sep='\t')
    df_clusters_enrichment = pd.read_csv(f'{file_pattern}/{disease_name}_clusters_enrichment_{cut_off}_networktype_{networktype.replace(" ", "_").lower()}.tsv', sep='\t')
    cluster_proteins = pd.read_csv(f'{file_pattern}/{disease_name}_clusters_proteins_{cut_off}_networktype_{networktype.replace(" ", "_").lower()}.tsv', sep='\t')
    cluster_proteins = cluster_proteins.groupby('cluster_idx')['protein'].apply(list).tolist()
    return df_cat,df_clusters_enrichment, cluster_proteins