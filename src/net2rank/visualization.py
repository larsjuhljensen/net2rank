"""
Visualization module for network analysis and protein clustering.

This module provides functions for loading cluster records, plotting cluster networks,
and saving cluster enrichment results for protein network analysis.
"""

import os
import pandas as pd
import py4cytoscape as p4c
from typing import List, Tuple, Optional
from .p4c_tools import make_network, make_network_plot_with_pie_charts, cluster_enrichment_pipeline


def load_records(
    disease_name: str,
    predictions: str,
    file_pattern: str,
    cut_off: float = 0.7,
    networktype: str = 'full STRING network'
) -> Tuple[pd.DataFrame, pd.DataFrame, List[List[str]]]:
    """
    Load cluster records including predictions, enrichment results, and cluster proteins.
    
    Parameters
    ----------
    disease_name : str
        Name of the disease for file naming conventions
    predictions : str
        Path to the predictions file
    file_pattern : str
        Base directory path for cluster files
    cut_off : float, optional
        Network cutoff value, by default 0.7
    networktype : str, optional
        Type of network used, by default 'full STRING network'
        
    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame, List[List[str]]]
        - df_cat: DataFrame with prediction results
        - df_clusters_enrichment: DataFrame with cluster enrichment results
        - cluster_proteins: List of protein lists for each cluster
    """
    # Load predictions
    df_cat = pd.read_csv(predictions, sep='\t')
    
    # Format network type for filename
    networktype_formatted = networktype.replace(" ", "_").lower()
    
    # Load cluster enrichment results
    enrichment_filename = (
        f'{disease_name}_clusters_enrichment_{cut_off}_networktype_{networktype_formatted}.tsv'
    )
    enrichment_path = os.path.join(file_pattern, enrichment_filename)
    df_clusters_enrichment = pd.read_csv(enrichment_path, sep='\t')
    
    # Load cluster proteins
    proteins_filename = (
        f'{disease_name}_clusters_proteins_{cut_off}_networktype_{networktype_formatted}.tsv'
    )
    proteins_path = os.path.join(file_pattern, proteins_filename)
    cluster_proteins_df = pd.read_csv(proteins_path, sep='\t')
    
    # Convert to list of lists grouped by cluster index
    cluster_proteins = cluster_proteins_df.groupby('cluster_idx')['protein'].apply(list).tolist()
    
    return df_cat, df_clusters_enrichment, cluster_proteins


def plot_cluster_network(
    cluster_proteins: List[List[str]],
    idx: int,
    df: pd.DataFrame,
    disease_name: str,
    cutoff: float = 0.7
) -> str:
    """
    Plot the network for a specific cluster index with pie charts.
    
    Parameters
    ----------
    cluster_proteins : List[List[str]]
        List of protein lists for each cluster
    idx : int
        Index of the cluster to plot
    df : pd.DataFrame
        DataFrame containing protein category information
    disease_name : str
        Name of the disease for network naming
    cutoff : float, optional
        Network cutoff value, by default 0.7
        
    Returns
    -------
    str
        Network SUID (Session Unique Identifier)
    """
    # Validate cluster index
    if idx >= len(cluster_proteins):
        raise IndexError(f"Cluster index {idx} out of range. Available clusters: 0-{len(cluster_proteins)-1}")
    
    # Create network for the specific cluster
    network_name = f'{disease_name}_cluster_{idx}_network'
    suid = make_network(
        proteins=cluster_proteins[idx],
        network_name=network_name,
        cutoff=cutoff
    )
    
    if 'aortic_aneurysm' in disease_name.lower():
        map_node_size = True
    else:
        map_node_size = False
    
    # Add pie charts to the network
    suid = make_network_plot_with_pie_charts(
        df_protein_category=df,
        selected_proteins=cluster_proteins[idx],
        suid=suid,
        map_node_size=map_node_size
    )
    
    # Apply layout and display
    p4c.commands.commands_run(f'layout apply preferred networkSelected={suid}')
    p4c.notebook_export_show_image()
    
    return suid


def save_clusters_and_enrichment(
    disease_name: str,
    predictions: str,
    save_to_dir: str,
    cutoff: float = 0.7,
    networktype: str = 'full STRING network',
    topk: Optional[int] = None
) -> None:
    """
    Save cluster enrichment results and cluster proteins to files.
    
    Parameters
    ----------
    disease_name : str
        Name of the disease for file naming
    predictions : str
        Path to the predictions file
    save_to_dir : str
        Directory to save the results
    cutoff : float, optional
        Network cutoff value, by default 0.7
    networktype : str, optional
        Type of network to create, by default 'full STRING network'
    topk : Optional[int], optional
        Number of top predictions to use, by default None (use all)
        
    Returns
    -------
    None
    """
    # Create output directory if it doesn't exist
    if not os.path.exists(save_to_dir):
        os.makedirs(save_to_dir)
    
    # Run the cluster enrichment pipeline
    df_cat, df_clusters_enrichment, cluster_proteins = cluster_enrichment_pipeline(
        disease_name=disease_name,
        predictions=predictions,
        cutoff=cutoff,
        networktype=networktype,
        topk=topk
    )
    
    # Format network type for filename
    networktype_formatted = networktype.replace(' ', '_').lower()
    
    # Save cluster enrichment results
    enrichment_filename = (
        f'{disease_name}_clusters_enrichment_{cutoff}_networktype_{networktype_formatted}.tsv'
    )
    enrichment_path = os.path.join(save_to_dir, enrichment_filename)
    df_clusters_enrichment.to_csv(enrichment_path, sep='\t', index=False)
    
    # Save cluster proteins
    proteins_filename = (
        f'{disease_name}_clusters_proteins_{cutoff}_networktype_{networktype_formatted}.tsv'
    )
    proteins_path = os.path.join(save_to_dir, proteins_filename)
    
    with open(proteins_path, 'w') as f:
        f.write('protein\tcluster_idx\n')
        for idx, cluster in enumerate(cluster_proteins):
            for protein in cluster:
                f.write(f'{protein}\t{idx}\n')
    
    print(f'Saved {disease_name} cluster enrichment and proteins to {save_to_dir}')