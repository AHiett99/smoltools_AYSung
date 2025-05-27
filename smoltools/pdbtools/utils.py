from Bio.PDB.Chain import Chain
from smoltools.pdbtools import load, pdb_select

def path_to_chain(path: str, model: int = 0, chain: str = 'A') -> Chain:
    """Extract the specified chain from a PDB file.

    Parameters:
    -----------
    path (str): Path to PDB file.
    model (int): Model number of desired chain (default = 0)
    chain (str): Chain ID of desired chain (default = 'A')

    Returns:
    --------
    Chain: PDB Chain object.
    """
    structure = load.read_pdb_from_path(path)
    return pdb_select.get_chain(structure, model=model, chain=chain)