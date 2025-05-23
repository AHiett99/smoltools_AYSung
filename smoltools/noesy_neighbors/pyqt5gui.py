import pandas as pd
from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QPushButton,
    QFileDialog, QListWidget, QListWidgetItem, QLabel, QGroupBox, QCheckBox, QHBoxLayout, QScrollArea, QDoubleSpinBox
)
from smoltools.noesy_neighbors import read_pdb_from_path, get_atom_names_by_residue
from smoltools.noesy_neighbors.chains import generate_df_chains, calculate_distances_by_chain
from smoltools.noesy_neighbors.utils import add_noe_bins
from smoltools.noesy_neighbors.plots import make_chart_title
from pathlib import Path
import altair as alt

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

            # Save for use later
            self.intra_chain_dfs = intra
            self.inter_chain_dfs = inter

            # print("Intra-chain distances:")
            # for chain_id, df in intra.items():
            #     print(f"Chain {chain_id}:")
            #     print(df.head(20))

            # print("Inter-chain distances:")
            # for (chain1, chain2), df in inter.items():
            #     print(f"Between {chain1} and {chain2}:")
            #     print(df.head(20))

            self.selector_label.setText(
                f"Calculated distances for {len(df_chains)} chains: "
                f"{len(intra)} intra-chain, {len(inter)} inter-chain."
            )

        except Exception as e:
            self.selector_label.setText(f"Error: {str(e)}")


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

    def get_structure(self):
        return self.structure   
    
    def get_data(self) -> dict[str, pd.DataFrame]:
        """
        Return a dictionary of all intra- and inter-chain distance DataFrames.
        Keys will indicate the chain or chain pair.
        """
        if not hasattr(self, 'intra_chain_dfs') or not hasattr(self, 'inter_chain_dfs'):
            print("No distances calculated yet.")
            return {}

        combined = {}

        for chain_id, df in self.intra_chain_dfs.items():
            combined[f"Intra-Chain distances: Atom 1 and Atom 2 = Chain {chain_id}"] = df

        for (chain1, chain2), df in self.inter_chain_dfs.items():
            combined[f"Inter-chain distances: Atom 1 = Chain {chain1} and Atom 2 = Chain {chain2}"] = df

        return combined


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

        # Validate ascending order (optional)
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
        try:
            bins, labels = self.bins.get_bins()
        except ValueError as e:
            print("Error:", e)
            return

        #calculate distances
        print("Calculating distance data now...")
        self.selector.calculate_distances()
        
        dfs = self.selector.get_data()

        charts = []
        for chain, df in dfs.items():
            
            df_binned = add_noe_bins(df, bins, labels)
            n_x = df_binned['id_1'].nunique()
            n_y = df_binned['id_2'].nunique()

            size_per_atom = 10
            max_size = 1000
            width = min(n_x * size_per_atom, max_size)
            height = min(n_y * size_per_atom, max_size)

            title = make_chart_title(chain)

            # Your plotting code here, e.g. saving to HTML and opening in browser
            chart = alt.Chart(df_binned).mark_rect().encode(
                x=alt.X('id_1', title='Atom #1', sort=None, axis=alt.Axis(labelFontSize=5)),
                y=alt.Y('id_2', title='Atom #2', sort=None, axis=alt.Axis(labelFontSize=5)),
                color=alt.Color(
                    'noe_strength',
                    title='NOE Strength',
                    scale=alt.Scale(domain=labels, scheme='blues', reverse=True)
                ),
                tooltip=[
                    alt.Tooltip('id_1', title='Atom #1'),
                    alt.Tooltip('id_2', title='Atom #2'),
                    alt.Tooltip('distance', title='Distance (Å)', format='.2f'),
                    alt.Tooltip('noe_strength', title='NOE Strength'),
                ]
            ).properties(width=width, height=height, title=title)

            charts.append(chart)

        combined = alt.vconcat(*charts)
        combined.save('./noe_charts.html')
        print("Saved noe_charts.html")
            
        #open in browser:
        import webbrowser
        webbrowser.open('./noe_charts.html')