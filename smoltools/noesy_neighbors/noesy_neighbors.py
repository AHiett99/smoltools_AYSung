'''Find expected NOEs after importing a PDB structure'''
from PyQt5.QtWidgets import (QApplication, QWidget, QVBoxLayout, QPushButton, QFileDialog, QListWidget, QListWidgetItem, QLabel, QGroupBox, QHBoxLayout, QScrollArea, QDoubleSpinBox, QMessageBox)
import pandas as pd
from pathlib import Path
import panel as pn
import numpy as np
import sys
import altair as alt
alt.data_transformers.disable_max_rows()

from collections import defaultdict
import itertools
import scipy.spatial.distance as ssd

from Bio.PDB import PDBParser
from Bio.PDB.Structure import Structure
from Bio.PDB.Atom import Atom
from Bio.PDB.Residue import Residue
from Bio.PDB.Chain import Chain

#errors
class ChainNotFound(KeyError):
    def __init__(self, structure_id: str, model_id: str, chain_id: str):
        message = f'Chain {structure_id}/{model_id}/{chain_id} not in structure'
        super().__init__(message)

class NoResiduesFound(ValueError):
    def __init__(self) -> None:
        message = 'No residues matching filter criteria found.'
        super().__init__(message)


class NoAtomsFound(ValueError):
    def __init__(self) -> None:
        message = 'No atoms matching filter criteria found.'
        super().__init__(message)

#load PDB stuff
def convert_to_path(path: str) -> Path:
    if not isinstance(path, Path):
        return Path(path)
    else:
        return path
    
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

def _flatten_list(nested_list: list[list]) -> list:
    return list(itertools.chain(*nested_list))

def _validate_atoms(atoms=list[Atom]) -> list[Atom]:
    if not atoms:
        raise NoAtomsFound
    else:
        return atoms
    
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

    return get_atoms(residues, labeled_atoms)


#Stuff to calculate distances
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
    
#Stuff to discern chains
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

def calculate_distances_by_chain(df_dict: dict[str, pd.DataFrame]) -> tuple[dict, dict]:
    """Calculate intra- and inter-chain distances using your distance module."""
    intra = {}
    inter = {}

    chains = list(df_dict.keys())

    for i, chain_i in enumerate(chains):
        df_i = df_dict[chain_i]

        #Store one intra-chain distance matrix PER chain
        intra[chain_i] = pairwise_distances(df_i)

        #Inter-chain distances (pairwise)
        for j in range(i + 1, len(chains)):
            chain_j = chains[j]
            df_j = df_dict[chain_j]
            inter[(chain_i, chain_j)] = pairwise_distances(df_i, df_j)

    return intra, inter


#Export HTML files
def export_charts(titles, tables, charts, parent_widget=None):
    msg = QMessageBox()
    msg.setIcon(QMessageBox.Information)
    msg.setWindowTitle("Export Charts")
    msg.setText("Please select the folder where the HTML files will be saved. This may take a moment!")
    msg.exec_()
    
    #Ask where to save
    folder = QFileDialog.getExistingDirectory(parent_widget, "Select Folder to Save Charts")

    if folder:
        folder_path = Path(folder)
        
        charts_dir = folder_path / 'charts'
        charts_dir.mkdir(parents=True, exist_ok=True)

        tables_dir = folder_path / 'tables'
        tables_dir.mkdir(parents=True, exist_ok=True)

    for title, chart in zip(titles, charts):
        file_path = charts_dir / f"{title}_chart.html"
        chart.save(str(file_path))

    for title, table in zip(titles, tables):
        file_path = tables_dir / f"{title}_table.html"
        panel_table = pn.widgets.Tabulator(table, pagination=None)
        panel_table.save(str(file_path))
    
    msg = QMessageBox()
    msg.setIcon(QMessageBox.NoIcon)
    msg.setWindowTitle("Export Charts")
    msg.setText(f"Success! Saved {len(charts)} charts and {len(tables)} tables to {folder}.")
    msg.exec_()

class PDBAtomSelector(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.structure = None
        self.residue_atom_dict = {}
        self.atom_checkboxes = {}
        self._build_ui()

    def _build_ui(self):
        self.layout = QVBoxLayout()

        self.load_button = QPushButton("Load PDB File")
        self.load_button.clicked.connect(self.load_pdb)

        self.selector_label = QLabel("No PDB loaded.")
        self.scroll_area = QScrollArea()
        self.scroll_area.setWidgetResizable(True)

        self.layout.addWidget(self.load_button)
        self.layout.addWidget(self.selector_label)
        self.layout.addWidget(self.scroll_area)
        self.setLayout(self.layout)

    def load_pdb(self):
        pdb_file, _ = QFileDialog.getOpenFileName(
            self, "Select PDB File", "", "PDB Files (*.pdb);;All Files (*)"
        )
        if not pdb_file:
            return

        try:
            self.structure = read_pdb_from_path(Path(pdb_file))
            self.residue_atom_dict = get_atom_names_by_residue(self.structure)
            self.selector_label.setText(f"PDB Loaded: {Path(pdb_file).name}. Please select atoms to calculate expected NOEs.")
            self.populate_atom_lists()
        except Exception as e:
            self.selector_label.setText(f"Failed to load PDB: {e}")

    def populate_atom_lists(self):
        content = QWidget()
        atom_layout = QVBoxLayout(content)
        self.atom_lists = {}  # Maps residue name to QListWidget

        for residue_name, atoms in self.residue_atom_dict.items():
            group = QGroupBox(residue_name)
            vbox = QVBoxLayout()
            list_widget = QListWidget()
            list_widget.setSelectionMode(QListWidget.MultiSelection)
            for atom in atoms:
                item = QListWidgetItem(atom)
                list_widget.addItem(item)
            self.atom_lists[residue_name] = list_widget
            vbox.addWidget(list_widget)
            group.setLayout(vbox)
            atom_layout.addWidget(group)

        content.setLayout(atom_layout)
        self.scroll_area.setWidget(content)

    def get_labeled_atoms_dict(self) -> dict[str, list[str]]:
        atoms_dict = {}
        for resname, list_widget in self.atom_lists.items():
            selected_items = list_widget.selectedItems()
            selected_atoms = [item.text() for item in selected_items]
            if selected_atoms:
                atoms_dict[resname] = selected_atoms
        return atoms_dict
    
    def calculate_distances(self):
        labeled_atoms = self.get_labeled_atoms_dict()
        if not labeled_atoms:
            self.selector_label.setText("No atoms selected!")
            return

        structure = self.get_structure()
        if structure is None:
            self.selector_label.setText("No PDB loaded!")
            return

        try:
            df_chains = generate_df_chains(structure, labeled_atoms)
            intra, inter = calculate_distances_by_chain(df_chains)

            self.intra_chain_dfs = intra
            self.inter_chain_dfs = inter

            self.selector_label.setText(
                f"Calculated distances for {len(df_chains)} chains: "
                f"{len(intra)} intra-chain, {len(inter)} inter-chain."
            )

        except Exception as e:
            self.selector_label.setText(f"Error: {str(e)}")
    
    def get_structure(self):
        return self.structure  

def add_noe_bins(df: pd.DataFrame, bins, labels) -> pd.DataFrame:
    """Add NOE strength bins based on user bins."""
    df = df.assign(
        noe_strength=pd.cut(
            df.distance,
            bins=bins,
            include_lowest=True,
            labels=labels,
            ordered=True,
        )
    )
    df['noe_strength'] = df['noe_strength'].astype(str)
    return df

class BinThresholdWidget(QWidget):
    def __init__(self):
        super().__init__()

        layout = QVBoxLayout()

        self.strong_spin = QDoubleSpinBox()
        self.strong_spin.setRange(0.0, 100.0)
        self.strong_spin.setValue(5.0)
        layout.addWidget(QLabel("Strong threshold (Å):"))
        layout.addWidget(self.strong_spin)

        self.medium_spin = QDoubleSpinBox()
        self.medium_spin.setRange(0.0, 100.0)
        self.medium_spin.setValue(8.0)
        layout.addWidget(QLabel("Medium threshold (Å):"))
        layout.addWidget(self.medium_spin)

        self.weak_spin = QDoubleSpinBox()
        self.weak_spin.setRange(0.0, 100.0)
        self.weak_spin.setValue(10.0)
        layout.addWidget(QLabel("Weak threshold (Å):"))
        layout.addWidget(self.weak_spin)

        # Button to trigger plotting or updating bins
        self.plot_btn = QPushButton("Plot NOE Map")
        layout.addWidget(self.plot_btn)

        self.setLayout(layout)

    def get_bins(self):
        strong = self.strong_spin.value()
        medium = self.medium_spin.value()
        weak = self.weak_spin.value()

        # Validate ascending order
        if not (strong < medium < weak):
            raise ValueError("Thresholds must be strictly ascending: strong < medium < weak")

        bins = [0, strong, medium, weak, float('inf')]
        labels = ['strong', 'medium', 'weak', 'none']

        return bins, labels
    
class MainWindow(QWidget):
    def __init__(self):
        super().__init__()

        layout = QHBoxLayout()

        self.selector = PDBAtomSelector()
        self.bins = BinThresholdWidget()

        layout.addWidget(self.selector)
        layout.addWidget(self.bins)

        self.setLayout(layout)

        self.bins.plot_btn.clicked.connect(self.handle_plot)

    def handle_plot(self):
        #calculate distances
        bins, labels = self.bins.get_bins()
        print("Calculating distance data now...")
        self.selector.calculate_distances()
        
        tables = []
        charts = []
        titles = []
        for chain, df in self.selector.intra_chain_dfs.items():
            
            df_binned = add_noe_bins(df, bins, labels)
            n_x = df_binned['id_1'].nunique()
            n_y = df_binned['id_2'].nunique()

            size_per_atom = 20
            max_size = 2000
            width = min(n_x * size_per_atom, max_size)
            height = min(n_y * size_per_atom, max_size)

            title = 'Expected Intra-chain NOEs'

            #Make tables
            table = df_binned.rename(columns={'id_1':f'Atom #1: Chain {chain}', 'id_2':f'Atom #2: Chain {chain}', 'distance':'Distance', 'noe_strength':'NOE strength'})
            table = table[table['NOE strength'].str.contains('none') == False]
            mask = table['Distance'] == 0.0
            table = table[~mask]
            tables.append(table)

            #make charts
            chart = alt.Chart(df_binned).mark_rect().encode(
                x=alt.X('id_1', title=f'Chain {chain}', sort=None, axis=alt.Axis(labelFontSize=10)),
                y=alt.Y('id_2', title=f'Chain {chain}', sort=None, axis=alt.Axis(labelFontSize=10)),
                color=alt.Color(
                    'noe_strength',
                    title='NOE Strength',
                    scale=alt.Scale(domain=labels, scheme='blues', reverse=True)
                ),
                tooltip=[
                    alt.Tooltip('id_1', title=f'Atom #1: Chain {chain}'),
                    alt.Tooltip('id_2', title=f'Atom #2: Chain {chain}'),
                    alt.Tooltip('distance', title='Distance (Å)', format='.2f'),
                    alt.Tooltip('noe_strength', title='NOE Strength'),
                ]
            ).properties(width=width, height=height, title=title)
            chart_interactive = chart.interactive()
            charts.append(chart_interactive)
            titles.append(str(f'Intra-chain_{chain}'))

        for (chain1, chain2), df in self.selector.inter_chain_dfs.items():

            df_binned = add_noe_bins(df, bins, labels)
            table = df_binned.rename(columns={'id_1':f'Atom #1: Chain {chain1}', 'id_2':f'Atom #2: Chain {chain2}', 'distance':'Distance', 'noe_strength':'NOE strength'})
            n_x = df_binned['id_1'].nunique()
            n_y = df_binned['id_2'].nunique()

            size_per_atom = 20
            max_size = 2000
            width = min(n_x * size_per_atom, max_size)
            height = min(n_y * size_per_atom, max_size)

            title = 'Expected Inter-chain NOEs'

            #Make tables
            table = df_binned.rename(columns={'id_1':f'Atom #1: Chain {chain1}', 'id_2':f'Atom #2: Chain {chain2}', 'distance':'Distance', 'noe_strength':'NOE strength'})
            table = table[table['NOE strength'].str.contains('none') == False]
            mask = table['Distance'] == 0.0
            table = table[~mask]
            tables.append(table)

            #make charts
            chart = alt.Chart(df_binned).mark_rect().encode(
                x=alt.X('id_1', title=f'Chain {chain1}', sort=None, axis=alt.Axis(labelFontSize=10)),
                y=alt.Y('id_2', title=f'Chain {chain2}', sort=None, axis=alt.Axis(labelFontSize=10)),
                color=alt.Color(
                    'noe_strength',
                    title='NOE Strength',
                    scale=alt.Scale(domain=labels, scheme='blues', reverse=True)
                ),
                tooltip=[
                    alt.Tooltip('id_1', title=f'Atom #1: Chain {chain1}'),
                    alt.Tooltip('id_2', title=f'Atom #2: Chain {chain2}'),
                    alt.Tooltip('distance', title='Distance (Å)', format='.2f'),
                    alt.Tooltip('noe_strength', title='NOE Strength'),
                ]
            ).properties(width=width, height=height, title=title)
            chart_interactive = chart.interactive()
            charts.append(chart_interactive)
            titles.append(str(f'Inter-chains_{chain1}-{chain2}'))
        
        export_charts(titles, tables, charts, parent_widget=self)
        self.close()
        

#launch GUI
def main():
    app = QApplication(sys.argv)
    window = MainWindow()
    window.setWindowTitle("NOESY Neighbors - created by Andrew Sung")
    window.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()