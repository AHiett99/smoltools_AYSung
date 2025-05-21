"""Functions for loading PDB files."""
import io
from pathlib import Path
from collections import defaultdict

from Bio.PDB import PDBParser
from Bio.PDB.Structure import Structure
from Bio.PDB.Atom import Atom
from Bio.PDB.Residue import Residue

import smoltools.pdbtools.pdb_select as pdb_select
import pandas as pd


def convert_to_path(path: str) -> Path:
    if not isinstance(path, Path):
        return Path(path)
    else:
        return path


def read_pdb_from_bytes(id: str, pdb_bytes: bytes) -> Structure:
    """
    Reads pdb file into a Structure object.

    Parameters:
    -----------
    id (str): id of structure object.
    pdb_bytes (bytes): byte object containing file data.

    Returns:
    --------
    Structure: Structure object containing data from the PDB file.
    """
    pdb_stream = io.StringIO(pdb_bytes.decode('utf-8').replace('\r', '\n'))
    return PDBParser().get_structure(id, pdb_stream)


def read_pdb_from_path(pdb_path: Path | str) -> Structure:
    """
    Reads a pdb file into a Structure object.

    Parameters:
    -----------
    pdb_path (Path | str): path to pdb file.

    Returns:
    --------
    Structure: Structure object containing data from the PDB file.
    """
    pdb_path = convert_to_path(pdb_path)
    id = pdb_path.stem
    return PDBParser().get_structure(id, pdb_path)

def get_atom_names_by_residue(structure: Structure) -> dict[str, list[str]]:
    """
    Extract atom names grouped by residue type from a Structure.

    Parameters:
    -----------
    structure (Structure): Biopython Structure object.

    Returns:
    --------
    dict: Dictionary mapping residue name (e.g., 'ILE') to list of atom names (e.g., ['CD1', 'CG2']).
    """
    residue_atoms = defaultdict(set)

    for model in structure:
        for chain in model:
            for residue in chain:
                # Skip heteroatoms and water
                hetfield, resseq, icode = residue.id
                if hetfield != ' ':
                    continue

                resname = residue.get_resname()
                for atom in residue:
                    residue_atoms[resname].add(atom.get_name())

    return {res: sorted(atoms) for res, atoms in residue_atoms.items()}

def get_labeled_atoms(
    residues: list[Residue], labeled_atoms: dict[str, list[str]]
) -> list[Atom]:
    """Retrieve labelled carbons from branched-chain amino acids (VAL, LEU, ILE)
    from a list of residues.

    Parameters:
    -----------
    residues (list[Residue]): List of PDB Residue objects.
    labeled_atoms (dict): Dictionary mapping three letter residue ID (e.g. 'ILE')
        to list of atoms to select (e.g. ['CD', 'CG2'])

    Returns:
    list[Atom]: List of PDB Atom objects.
    """

    return pdb_select.get_atoms(residues, labeled_atoms)


#commenting out this section since it will probably become obselete
# def coordinates_from_path_presets(
#     path: str,
#     mode: str = 'ILV',
#     model: int = 0,
#     chain: str = 'A',
# ) -> pd.DataFrame:
#     """Calculate pairwise distances of terminal carbons of branched-chain amino acids
#     in the specified chain from a PDB file. Use if starting directly from PDB file.

#     Parameters:
#     -----------
#     path (str): Path to PDB file.
#     mode (str): Predefined labeled atom selections (choices are 'ILV', 'ILVA', and 'ILVMAT')
#     model (int): Model number of desired chain (default = 0)
#     chain (str): Chain ID of desired chain (default = 'A')

#     Returns:
#     --------
#     DataFrame: Dataframe with the atom IDs (residue number, carbon ID) of each atom pair
#         and the distance (in angstroms) between each pair.
#     """
#     labeled_atoms = PREDEFINED_LIST_OF_LABELED_ATOMS[mode]
#     chain = path_to_chain(path, model=model, chain=chain)
#     return coordinates_from_chain(chain, labeled_atoms)