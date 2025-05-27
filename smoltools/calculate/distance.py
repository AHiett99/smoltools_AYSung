"""Functions for calculating distances between atoms in a PDB structure."""

import numpy as np
import pandas as pd
import scipy.spatial.distance as ssd


def _pairwise_distance(df_a: pd.DataFrame, df_b: pd.DataFrame) -> np.ndarray:
    """Return the euclidean distance between all 3D coordinates."""
    return ssd.cdist(df_a, df_b, 'euclidean')


def _tidy_pairwise_distances(df: pd.DataFrame) -> pd.DataFrame:
    """Take a square dataframe of pairwise distances and convert it to tidy format."""
    return df.melt(value_name='distance', ignore_index=False).reset_index()

def pairwise_distances(df_a: pd.DataFrame, df_b: pd.DataFrame = None) -> pd.DataFrame:
    """Given two dataframes with 3D coordinates of each residue, calculate the pairwise
    distance between each residue and return in tidy form.
    """
    if df_b is None:
        df_b = df_a
    return (
        pd.DataFrame(
            _pairwise_distance(df_a, df_b),
            index=df_a.index,
            columns=df_b.index,
        )
        .rename_axis(index='id_1', columns='id_2')
        .pipe(_tidy_pairwise_distances)
    )