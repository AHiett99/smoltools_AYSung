"""Convert list of residues to a table of atomic coordinates"""
from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
import smoltools.pdbtools.pdb_select as pdb_select
from smoltools.pdbtools.load import get_labeled_atoms
from smoltools.pdbtools import path_to_chain
import pandas as pd

def coordinates_from_chain(
    chain: Chain, labeled_atoms: dict[str, list[str]]
) -> pd.DataFrame:
    """Calculate pairwise distances of terminal carbons of branched-chain amino acids
    in the given Chain object. Use if a chain object is already loaded.

    Parameters:
    -----------
    chain (Chain): PDB Chain object.
    labeled_atoms (dict): Dictionary mapping three letter residue ID (e.g. 'ILE')
        to list of atoms to select (e.g. ['CD', 'CG2'])

    Returns:
    --------
    DataFrame: Dataframe with the atom IDs (residue number, carbon ID) of each atom pair
        and the distance (in angstroms) between each pair.
    """
    residue_filter = set(labeled_atoms.keys())
    residues = pdb_select.get_residues(chain, residue_filter=residue_filter)
    atoms = get_labeled_atoms(residues, labeled_atoms)
    return (
        coordinate_table(atoms)
        .assign(
            id=lambda x: x.residue_name + x.residue_number.astype(str) + '-' + x.atom_id
        )
        .set_index('id')
        .loc[:, ['x', 'y', 'z']]
    )

def coordinates_from_path(
    path: str,
    labeled_atoms: dict[str, list[str]],
    model: int = 0,
    chain: str = 'A',
) -> pd.DataFrame:
    """Calculate pairwise distances of terminal carbons of branched-chain amino acids
    in the specified chain from a PDB file. Use if starting directly from PDB file.

    Parameters:
    -----------
    path (str): Path to PDB file.
    labeled_atoms (dict): Dictionary mapping three letter residue ID (e.g. 'ILE')
        to list of atoms to select (e.g. ['CD', 'CG2'])
    model (int): Model number of desired chain (default = 0)
    chain (str): Chain ID of desired chain (default = 'A')

    Returns:
    --------
    DataFrame: Dataframe with the atom IDs (residue number, carbon ID) of each atom pair
        and the distance (in angstroms) between each pair.
    """
    chain = path_to_chain(path, model=model, chain=chain)
    return coordinates_from_chain(chain, labeled_atoms)

def coordinate_table(atoms: list[Atom]) -> pd.DataFrame:
    """Extract 3D coordinates from list of atoms into DataFrame.

    Parameters:
    -----------
    atoms (list[Atom]): List of PDB Atom.

    Returns:
    --------
    DataFrame: Dataframe with the atom ID (residue number, carbon ID) as the index
        and the x, y, z coordinate of each atom as the columns.
    """

    def _get_atom_info(atom: Atom) -> tuple:
        parent_residue = atom.get_parent()
        residue_number = parent_residue.get_id()[1]
        residue_name = parent_residue.get_resname()
        atom_id = atom.get_id()
        coordinates = atom.get_coord()

        return residue_name, residue_number, atom_id, *coordinates

    atom_info = [_get_atom_info(atom) for atom in atoms]
    info_columns = ['residue_name', 'residue_number', 'atom_id']

    return pd.DataFrame(
        atom_info,
        columns=[*info_columns, 'x', 'y', 'z'],
    )
