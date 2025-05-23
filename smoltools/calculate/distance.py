"""Functions for calculating atomic distances."""

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


def _merge_pairwise_distances(df_a: pd.DataFrame, df_b: pd.DataFrame) -> pd.DataFrame:
    """Merge two DataFrames of pairwise distances (intersection of residues pairs in
    each dataset)
    """
    return pd.merge(
        df_a,
        df_b,
        on=['id_1', 'id_2'],
        suffixes=['_a', '_b'],
        validate='1:1',
    )


def pairwise_distances_between_conformations(
    distances_a: pd.DataFrame, distances_b: pd.DataFrame
) -> pd.DataFrame:
    """Given two DataFrames with the pairwise distances between each residue in two
    conformation, return a merged DataFrame that also contains the difference in
    pairwise distances between the conformations.

    Parameters:
    -----------
    distances_a (DataFrame): Dataframe with the atom IDs (residue number, carbon ID)
        of each atom pair and the distance (in angstroms) between each pair.
    distances_b (DataFrame): Dataframe with the atom IDs (residue number, carbon ID)
        of each atom pair and the distance (in angstroms) between each pair.

    Returns:
    --------
    DataFrame: DataFrame of the pairwise distance between residues for each
        of the two conformations, and the difference in the pairwise distances
        between the conformations. Distances reported in angstroms.
    """
    return _merge_pairwise_distances(
        distances_a,
        distances_b,
    ).assign(delta_distance=lambda x: x.distance_a - x.distance_b)
