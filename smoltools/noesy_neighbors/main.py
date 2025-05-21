#use to launch GUI
import sys
from PyQt5.QtWidgets import QApplication
from smoltools.noesy_neighbors.pyqt5gui import PDBAtomSelector, MainWindow

def main():
    app = QApplication(sys.argv)
    window = MainWindow()
    window.setWindowTitle("PDB Selector and Bin Thresholds")
    window.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()