import pandas as pd
import numpy as np
import smoltools.pdbtools.pdb_select as pdb_select
from smoltools.pdbtools.load import get_labeled_atoms
from smoltools.pdbtools import coordinate_table
from smoltools.calculate.distance import pairwise_distances

def generate_df_chains(structure, labeled_atoms: dict[str, list[str]]) -> dict[str, pd.DataFrame]:
    """Extract coordinate DataFrames per chain from a structure using selected atoms."""
    df_dict = {}
    for model in structure.get_models():
        for chain in model.get_chains():
            residues = pdb_select.get_residues(chain, residue_filter=set(labeled_atoms.keys()))
            atoms = get_labeled_atoms(residues, labeled_atoms)
            if atoms:
                df = (
                    coordinate_table(atoms)
                    .assign(id=lambda x: x.residue_name + x.residue_number.astype(str) + '-' + x.atom_id)
                    .set_index('id')[['x', 'y', 'z']]
                )
                df_dict[chain.id] = df
    return df_dict

def calculate_distances_by_chain(df_dict: dict[str, pd.DataFrame]) -> tuple[dict, dict]:
    """Calculate intra- and inter-chain distances using your distance module."""
    intra = {}
    inter = {}

    chains = list(df_dict.keys())

    for i, chain_i in enumerate(chains):
        df_i = df_dict[chain_i]

        # ✅ Store one intra-chain distance matrix PER chain
        intra[chain_i] = pairwise_distances(df_i)

        # Inter-chain distances (pairwise)
        for j in range(i + 1, len(chains)):
            chain_j = chains[j]
            df_j = df_dict[chain_j]
            inter[(chain_i, chain_j)] = pairwise_distances(df_i, df_j)

    return intra, inter

def _pairwise_distance(df_a: pd.DataFrame, df_b: pd.DataFrame) -> np.ndarray:
    """Return the euclidean distance between all 3D coordinates."""
    return ssd.cdist(df_a, df_b, 'euclidean')