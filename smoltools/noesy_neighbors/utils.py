import numpy as np
import pandas as pd


def extract_residue_number(s: pd.Series) -> pd.Series:
    return s.str.partition('-')[0].str[3:].astype(int)


def lower_triangle(df: pd.DataFrame) -> pd.Series:
    return extract_residue_number(df.id_1) < extract_residue_number(df.id_2)


def splice_conformation_tables(
    df_a: pd.DataFrame,
    df_b: pd.DataFrame,
    chain_a_id: str = 'A',
    chain_b_id: str = 'B',
) -> pd.DataFrame:
    """Splice distance tables for two conformations together.

    Parameters:
    -----------
    df_a (DataFrame): Dataframe with the atom IDs (residue number, carbon ID) of each
        atom pair and the distance (in angstroms) between each pair.
    df_b (DataFrame): Dataframe with the atom IDs (residue number, carbon ID) of each
        atom pair and the distance (in angstroms) between each pair.

    Returns:
    --------
    DataFrame: DataFrame with the lower triangle of the DataFrame containing values
        from the first conformation and the upper triangle of the DataFrame containing
        values from the second conformation.
    """
    overlap = set(df_a.id_1) & set(df_b.id_1)
    return (
        pd.concat(
            [
                df_a.loc[
                    lambda x: extract_residue_number(x.id_1)
                    <= extract_residue_number(x.id_2)
                ].assign(subunit=chain_a_id),
                df_b.loc[
                    lambda x: extract_residue_number(x.id_1)
                    > extract_residue_number(x.id_2)
                ].assign(subunit=chain_b_id),
            ]
        )
        .loc[lambda x: x.id_1.isin(overlap) & x.id_2.isin(overlap)]
        .sort_values(['id_1', 'id_2'], key=extract_residue_number)
        .astype({'subunit': 'category'})
    )


def add_noe_bins(df: pd.DataFrame, bins, labels) -> pd.DataFrame:
    """Add NOE strength bins based on user bins."""
    return df.assign(
        noe_strength=pd.cut(
            df.distance,
            bins=bins,
            include_lowest=True,
            labels=labels,
            ordered=True,
        )
    )