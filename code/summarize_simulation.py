"""
Summarizes PhysiCell ABM output with spatial statistics from the Squidpy package.
"""

import numpy as np
import squidpy as sq
import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import pcdl, os, pickle, itertools, glob

def calc_perc_celltype(x):
    cnames = list(x.obs['cell_type'].value_counts(normalize=True).index)
    values = (x.obs['cell_type'].value_counts(normalize=True) * 100).values
    df = pd.DataFrame([values], columns=cnames)
    return(df)

def calc_im(x):
    values = sq.gr.interaction_matrix(x, cluster_key = "cell_type", copy=True).ravel()
    cell_types = list(x.obs['cell_type'].cat.categories)
    cnames = ["im_"+i for i in ["_".join(x) for x in itertools.product(cell_types, cell_types)]]
    df = pd.DataFrame([values], columns = cnames)
    return(df)

def calc_nes(x):
    values = sq.gr.nhood_enrichment(x, cluster_key = "cell_type", copy = True, seed=123, n_jobs=8, n_perms=500)[0].ravel()
    cell_types = list(x.obs['cell_type'].cat.categories)
    cnames = ["nes_"+i for i in ["_".join(x) for x in itertools.product(cell_types, cell_types)]]
    df = pd.DataFrame([values], columns = cnames)
    return(df)

def calc_centrality(x):
    values = sq.gr.centrality_scores(x, cluster_key = "cell_type", copy = True).to_numpy().ravel()
    cell_types = list(sq.gr.centrality_scores(x, cluster_key = "cell_type", copy = True).index)
    measures = list(sq.gr.centrality_scores(x, cluster_key = "cell_type", copy = True).columns)
    cnames = ["_".join(item) for item in itertools.product(cell_types, measures)]
    df = pd.DataFrame([values], columns = cnames)
    return(df)

def append_dataframes_with_zeros(dataframes):
    # Step 1: Identify all unique columns across all dataframes
    all_columns = set()
    for df in dataframes:
        all_columns.update(df.columns)

    # Convert the set of all columns to a sorted list
    all_columns = sorted(all_columns)

    # Step 2: Create a list to store modified dataframes
    modified_dfs = []

    for df in dataframes:
        # Step 3: Create a DataFrame with all columns filled with zeros
        df_with_all_columns = pd.DataFrame(0, index=df.index, columns=all_columns)
        
        # Step 4: Update the DataFrame with the original data
        for col in df.columns:
            df_with_all_columns[col] = df[col]
        
        # Step 5: Add the modified DataFrame to the list
        modified_dfs.append(df_with_all_columns)
    
    # Step 6: Concatenate all modified dataframes into a single dataframe
    result_df = pd.concat(modified_dfs, ignore_index=True)

    return result_df

def summarize_simulation(sim_dir, **kwargs):
    all_timestep_output_fnames = sorted(glob.glob(f"{sim_dir}/output*.xml"))
    config_file = kwargs.get('config.xml', 'config.xml')
    mcds_all = [pcdl.TimeStep(f, settingxml = config_file) for f in all_timestep_output_fnames]

    # get the anndata for each of the timesteps
    adata_all = [x.get_anndata(values=2, scale='maxabs') for x in mcds_all]

    # calculate the connectivity graph for each of the anndata objects
    for i, ad in enumerate(adata_all):
        ad.obs['cell_type'] = ad.obs['cell_type'].astype("category")
        sq.gr.spatial_neighbors(ad)

    centralities = [calc_centrality(x) for x in adata_all]
    centralities_df = append_dataframes_with_zeros(centralities)
    nes = [calc_nes(x) for x in adata_all]
    nes_df = append_dataframes_with_zeros(nes)
    imd = [calc_im(x) for x in adata_all]
    imd_df = append_dataframes_with_zeros(imd)
    perc_ct = [calc_perc_celltype(x) for x in adata_all]
    perc_ct_df = append_dataframes_with_zeros(perc_ct)

    feature_matrix = pd.concat([centralities_df, nes_df, imd_df, perc_ct_df], axis=1)
    return(feature_matrix)