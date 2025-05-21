"""Functions for selecting residues and atoms from PDB structure."""

import itertools

from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure

from smoltools.pdbtools.exceptions import ChainNotFound, NoResiduesFound, NoAtomsFound


def get_chain(structure: Structure, model: int, chain: str) -> Chain:
    """Returns a chain from a PDB structure object.

    Parameters:
    -----------
    structure (Structure): PDB structure object.
    model (int): Model number.
    chain (str): Chain identifier.

    Returns:
    --------
    Chain: PDB chain object.
    """
    try:
        return structure[model][chain]
    except KeyError as e:
        raise ChainNotFound(structure.get_id(), model, chain) from e


def get_residues(chain: Chain, residue_filter: set[str] = None) -> list[Residue]:
    """Produces a list of all residues in a PDB chain. Can provide a set of specific
    residues to keep.

    Parameters:
    -----------
    chain (Chain): PDB chain object.
    residue_filter (set[str]): Optional, a set (or other list-like) of three letter
        amino codes for the residues to keep. Default is to return all residues.

    Returns:
    --------
    list[Residue]: List of PDB residue objects in the given entity that meet the
        residue filter.
    """
    if residue_filter is None:
        residues = [
            residue for residue in chain.get_residues() if residue.get_id()[0] == ' '
        ]
    else:
        residues = [
            residue
            for residue in chain.get_residues()
            if residue.get_resname() in residue_filter
        ]

    return _validate_residues(residues)


def _validate_residues(residues: list[Residue]) -> list[Residue]:
    if not residues:
        raise NoResiduesFound
    else:
        return residues

def get_atoms(
    residues: list[Residue], atom_select: dict[str, list[str]]
) -> list[Atom]:
    """
    Returns a list of atoms from a list of residues that match the provided
    atom selection criteria.

    Parameters:
    -----------
    residues (list[Residue]): List of PDB residue objects.
    atom_select (dict): Dictionary of residue names (3-letter codes) to a list of
                        atom names to include. E.g., {'ILE': ['CD1'], 'LEU': ['CD2']}.

    Returns:
    --------
    list[Atom]: List of selected PDB Atom objects.
    """

    def _get_atoms(residue: Residue) -> list[Atom]:
        resname = residue.get_resname()
        if resname not in atom_select:
            return []
        atom_filter = atom_select[resname]
        return [atom for atom in residue.get_atoms() if atom.get_name() in atom_filter]

    atoms = _flatten_list([_get_atoms(residue) for residue in residues])
    return _validate_atoms(atoms)
              
# def get_alpha_carbons(residues: list[Residue]) -> list[Atom]:
#     """Returns a list of alpha carbons for a given list of residues.

#     Parameters:
#     -----------
#     residues (list[Residue]): list of PDB residue objects.

#     Returns:
#     --------
#     list[Atom]: list of alpha carbons as PDB atom objects.
#     """

#     def _get_atoms(residue) -> list[Atom]:
#         return [atom for atom in residue.get_atoms() if atom.get_id() == 'CA']

#     atoms = _flatten_list([_get_atoms(residue) for residue in residues])

#     if not atoms:
#         raise NoAtomsFound
#     else:
#         return atoms


# def get_carbons(
#     residues: list[Residue], atom_select: dict[str : list[str]]
# ) -> list[Atom]:
#     """Returns a list of atoms from a list of residues that meet the atom selection
#     criteria. Requires a dictionary of the names of the atoms to retrieve for each
#     amino acid.
#     """

#     def _get_atoms(residue: Residue):
#         atom_filter = atom_select[residue.get_resname()]
#         return [atom for atom in residue.get_atoms() if atom.get_name() in atom_filter]

#     atoms = _flatten_list([_get_atoms(residue) for residue in residues])
#     return _validate_atoms(atoms)


# def filter_by_b_factor(atoms: list[Atom], cutoff) -> list[Atom]:
#     """Returns a list of atoms with a b factor that meets the provided cutoff."""
#     atoms = [atom for atom in atoms if atom.get_bfactor() > cutoff]
#     return _validate_atoms(atoms)


def _validate_atoms(atoms=list[Atom]) -> list[Atom]:
    if not atoms:
        raise NoAtomsFound
    else:
        return atoms


def _flatten_list(nested_list: list[list]) -> list:
    return list(itertools.chain(*nested_list))
