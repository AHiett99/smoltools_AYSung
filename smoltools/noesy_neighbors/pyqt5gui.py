'''gui logic & plotting'''
from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QPushButton, QFileDialog, QListWidget, QListWidgetItem, QLabel, QGroupBox, QHBoxLayout, QScrollArea, QDoubleSpinBox, QMessageBox)
import pandas as pd
from smoltools.pdbtools.load import read_pdb_from_path, get_atom_names_by_residue
from smoltools.noesy_neighbors.chains import generate_df_chains, calculate_distances_by_chain
from pathlib import Path
import panel as pn
import altair as alt
alt.data_transformers.disable_max_rows()

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