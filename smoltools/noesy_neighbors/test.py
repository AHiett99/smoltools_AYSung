#can the pdb be read in?
from Bio.PDB import PDBParser
from Bio.PDB.PDBIO import PDBIO

from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure
import itertools
import pandas as pd

parser = PDBParser()
test = parser.get_structure('test1', '/Users/ahiett/Desktop/KH34 AF3 FUBP3 Predictions/kh34_model0'
'/kh34_model0.pdb')

# io = PDBIO()
# output_test = io.set_structure(test)
# io.save('test.pdb')


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
        print('nope not today')
    else:
        return residues
    
blblbl = get_residues(test)
#print(blblbl)

def _validate_atoms(atoms=list[Atom]) -> list[Atom]:
    if not atoms:
        print('ugly')
    else:
        return atoms


def _flatten_list(nested_list: list[list]) -> list:
    return list(itertools.chain(*nested_list))

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

    return get_atoms(residues, labeled_atoms)

def coordinate_table(atoms: list[Atom]) -> pd.DataFrame:
    """Extract 3D coordinates from list of atoms into DataFrame.

    Parameters:
    -----------
    atoms (list[Atom]): List of PDB Atom.

    Returns:
    --------
    DataFrame: Dataframe with the atom ID (residue number, atom ID) as the index
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

def generate_df_chains(structure, labeled_atoms: dict[str, list[str]]) -> dict[str, pd.DataFrame]:
    """Extract coordinate DataFrames per chain from a structure using selected atoms."""
    df_dict = {}
    for model in structure.get_models():
        for chain in model.get_chains():
            residues = get_residues(chain, residue_filter=set(labeled_atoms.keys()))
            atoms = get_labeled_atoms(residues, labeled_atoms)
            if atoms:
                df = (
                    coordinate_table(atoms)
                    .assign(id=lambda x: x.residue_name + x.residue_number.astype(str) + '-' + x.atom_id)
                    .set_index('id')[['x', 'y', 'z']]
                )
                df_dict[chain.id] = df
    return df_dict

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
        print('donuts')

ahh = get_chain(test, 0, 'A')
weee = get_residues(test)
beans = generate_df_chains(weee, labeled_atoms={'N': 'tester'})
print(beans)