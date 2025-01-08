import numpy as np
import pandas as pd
import os
import sys

# run from the 'code' directory
DATA_DIR = "../output/tables/parameter_sampling"
WRITE_DIR = "../output/tables"
VERSION = "1"

def unravel_dataframe(df):
    """
    Unravels a Pandas DataFrame with time series data into a single-row 
    DataFrame with column names prepended with "t{i}_".

    Args:
        df: The input DataFrame with timepoints as rows and features as columns.

    Returns:
        A single-row DataFrame with modified column names, or None if the input
        is not a DataFrame.
    """
    if not isinstance(df, pd.DataFrame):
        print("Error: Input must be a Pandas DataFrame.")
        return None

    n_timepoints = df.shape[0]
    feature_names = df.columns

    # Create new column names
    new_columns = []
    for timepoint in range(n_timepoints):
        for feature in feature_names:
            new_columns.append(f"t{timepoint}_{feature}")

    # Create a new DataFrame with a single row and the new column names
    unraveled_df = pd.DataFrame(index=[0], columns=new_columns)

    # Fill the new DataFrame
    k=0
    for timepoint in range(n_timepoints):
        for feature in feature_names:
            unraveled_df.iloc[0, k] = df.loc[timepoint, feature] # Assign values from the original df
            k+=1


    return unraveled_df

if __name__ == "__main__":
    # for converting to a single dataframe
    csvs = [f for f in os.listdir(DATA_DIR) if ".csv" in f]
    dfs = [pd.read_csv(os.path.join(DATA_DIR, f)) for f in csvs]
    dfs1d = [unravel_dataframe(df) for df in dfs]
    df = pd.concat(dfs1d)
    df = df.drop(df.filter(regex='Unnamed').columns, axis=1)
    # set nans and inf to 0 as this indicates a given cell type in the statistic was not present at the timestep 
    df = df.replace([np.inf, -np.inf], np.nan) 
    df = df.replace(np.nan, 0)
    df = df.assign(sim_id=pd.Series(csvs).values)
    df.to_csv(os.path.join(WRITE_DIR, f"simulations_{VERSION}.csv"), index=False)
