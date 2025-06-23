
import sys
import os
import copy
from rdkit.Chem import  AllChem, rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from PySide6.QtWidgets import  QApplication, QFileDialog, QDialog, QVBoxLayout, QTextEdit
from PySide6 import QtCore, QtGui, QtWidgets
from PySide6.QtCore import QProcess
from PySide6.QtGui import QIcon,  QKeySequence
from molEditWidget import MolEditWidget
from ptable_widget import PTable
from rdkit import Chem
import shutil
import numpy as np



if getattr(sys, 'frozen', False):
   import pyi_splash

mainpath = os.path.dirname(__file__)
appicon = os.path.join(mainpath, r'Img\MS-Draw.ico')
apppng = os.path.join(mainpath, r'Img\MS-logo.png')
fileopen = os.path.join(mainpath, r'Img\File.png')
drawpng = os.path.join(mainpath, r'Img\draw.png')
cleanstr = os.path.join(mainpath, r'Img\clean.png')
cursormode = os.path.join(mainpath, r'Img\cursor_mode.png')
removeatom = os.path.join(mainpath, r'Img\Remove.png')
Stereopng = os.path.join(mainpath, r'Img\Stereo.png')
Stereoezpng = os.path.join(mainpath, r'Img\ez.png')
tablepng = os.path.join(mainpath, r'Img\table.png')
trashpng = os.path.join(mainpath, r'Img\trash.png')
xypng = os.path.join(mainpath, r'Img\xy.png')
pluschargepng = os.path.join(mainpath, r'Img\pluscharge.png')
negativechargepng = os.path.join(mainpath, r'Img\minuscharge.png')
hashmappng = os.path.join(mainpath, r'Img\Hashmap.png')
singlebondpng = os.path.join(mainpath, r'Img\singlebond.png')
doublebondpng = os.path.join(mainpath, r'Img\double.png')
triplebondpng = os.path.join(mainpath, r'Img\Triplebond.png')
benzene = os.path.join(mainpath, r'Img\benzene.png')
copypng = os.path.join(mainpath, r'Img\copy.png')
cyclohexane = os.path.join(mainpath, r'Img\cyclohexane')
cyclopentane = os.path.join(mainpath, r"Img\cyclopentane.png")
cyclobutane = os.path.join(mainpath, r"Img\cyclobutane.png")
cyclopropane = os.path.join(mainpath, r"Img\cyclopropane.png")
output_dir = os.path.join(mainpath, r"OpenBabel\output")
mol2batch = os.path.join(mainpath, r"OpenBabel\mol2prep.bat")
pdbqtbatch = os.path.join(mainpath, r"OpenBabel\pdbqtprep.bat")


class BatchProcessMOL2Dialog(QDialog):
    def __init__(self, batch_file_path, forcefield, protonation_state, steps):
        super().__init__()
        self.setWindowTitle("Energy Minimization")
        self.setWindowIcon(QIcon.fromTheme(appicon))
        self.resize(600, 400)

        # Text box to show logs
        self.text_edit = QTextEdit(self)
        self.text_edit.setStyleSheet("""
    QTextEdit {
        background-color: #1e1e1e;
        color: #dcdcdc;
        font-family: Consolas, Courier New, monospace;
        font-size: 14px;
        border: none;
        padding: 5px;
    }
""")
        self.text_edit.setReadOnly(True)
        self.log_data = ""

        layout = QVBoxLayout(self)
        layout.addWidget(self.text_edit)
        forcefield = str(forcefield)
        protonation_state = str(protonation_state)
        steps = str(steps)

        # QProcess setup
        self.process = QProcess(self)
        self.process.setProgram("cmd.exe")
        self.process.setArguments(["/c", batch_file_path, forcefield, protonation_state, steps])  # /c means run then exit
        self.text_edit.append(f"FORCEFIELD: {forcefield}\n PROTONATION_STATE: {protonation_state} \n STEPS: {steps} \n OUTPUT: MOL2")
        self.process.readyReadStandardOutput.connect(self.handle_stdout)
        self.process.readyReadStandardError.connect(self.handle_stderr)
        self.process.finished.connect(self.on_finished)

        # Start the batch process
        self.append_log(f"Starting: {batch_file_path}")
        self.append_log(f"Arguments:\n  Forcefield: {forcefield}\n  Protonation: {protonation_state}\n  Steps: {steps}")
        self.process.start()


    def append_log(self, text):
        self.log_data += text + "\n"
        self.text_edit.append(text)

    def handle_stdout(self):
        text = self.process.readAllStandardOutput().data().decode()
        self.text_edit.append(text)
        self.append_log(f"{text}")

    def handle_stderr(self):
        text = self.process.readAllStandardError().data().decode()
        self.text_edit.append(f"{text}")
        self.append_log(f"{text}")

    def on_finished(self):
        self.text_edit.append("\nBatch process completed.")
        print(self.log_data)  # Print the log data to console
        self.accept() 


class BatchProcessPDBQTDialog(QDialog):
    def __init__(self, batch_file_path, forcefield, protonation_state, steps):
        super().__init__()
        self.setWindowTitle("Energy Minimization")
        self.setWindowIcon(QIcon.fromTheme(appicon))
        self.resize(600, 400)

        # Text box to show logs
        self.text_edit = QTextEdit(self)
        self.text_edit.setStyleSheet("""
    QTextEdit {
        background-color: #1e1e1e;
        color: #dcdcdc;
        font-family: Consolas, Courier New, monospace;
        font-size: 14px;
        border: none;
        padding: 5px;
    }
""")
        self.text_edit.setReadOnly(True)

        layout = QVBoxLayout(self)
        layout.addWidget(self.text_edit)

        # QProcess setup
        forcefield = str(forcefield)
        protonation_state = str(protonation_state)
        steps = str(steps)
        

        self.process = QProcess(self)
        self.text_edit.append(f"FORCEFIELD: {forcefield}\n PROTONATION_STATE: {protonation_state} \n STEPS: {steps} \n OUTPUT: PDBQT")
        self.process.setProgram("cmd.exe")
        self.process.setArguments(["/c", batch_file_path, forcefield, protonation_state, steps])  # /c means run then exit

        self.process.readyReadStandardOutput.connect(self.handle_stdout)
        self.process.readyReadStandardError.connect(self.handle_stderr)
        self.process.finished.connect(self.on_finished)

        # Start the batch process
        self.process.start()

    def handle_stdout(self):
        text = self.process.readAllStandardOutput().data().decode()
        self.text_edit.append(text)

    def handle_stderr(self):
        text = self.process.readAllStandardError().data().decode()
        self.text_edit.append(f"{text}")

    def on_finished(self):
        self.text_edit.append("\nBatch process completed.")
        self.accept() 
    
class StyledMessageDialog(QtWidgets.QDialog):
    def __init__(self, title, message, parent=None):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.setModal(True)
        self.setFixedSize(350, 150)
        self.setStyleSheet("""
            QDialog {
                background-color: #232629;
                border-radius: 10px;
            }
            QLabel {
                color: #f8f8f2;
                font-size: 14px;
                padding: 10px;
            }
            QPushButton {
                background-color: #6272a4;
                color: #f8f8f2;
                border-radius: 5px;
                padding: 6px 18px;
                font-size: 13px;
            }
            QPushButton:hover {
                background-color: #44475a;
            }
        """)
        layout = QtWidgets.QVBoxLayout(self)
        label = QtWidgets.QLabel(message)
        label.setWordWrap(True)
        layout.addWidget(label)
        btn = QtWidgets.QPushButton("OK")
        btn.clicked.connect(self.accept)
        btn.setDefault(True)
        layout.addWidget(btn, alignment=QtCore.Qt.AlignCenter)



class MainWindow(QtWidgets.QMainWindow):
    # Constructor function
    def __init__(self, fileName=None, loglevel="WARNING"):
        super(MainWindow, self).__init__()
        
        self.pixmappath = os.path.abspath(os.path.dirname(__file__)) + "/pixmaps/"
        QtGui.QIcon.setThemeSearchPaths(
            # QtGui.QIcon.themeSearchPaths() +
            [os.path.abspath(os.path.dirname(__file__)) + "/icon_themes/"]
        )
        self.loglevels = ["Critical", "Error", "Warning", "Info", "Debug", "Notset"]
        # RDKit draw options, tooltip, default value is read from molViewWidget
        self._drawopts_actions = [
            (
                "prepareMolsBeforeDrawing",
                "Prepare molecules before drawing (i.e. fix stereochemistry and annotations)",
            ),
            (
                "addStereoAnnotation",
                "Add stereo annotation (R/S and E/Z)",
            ),
            (
                "unspecifiedStereoIsUnknown",
                "Show wiggly bond at potential undefined chiral stereo centres "
                + "and cross bonds for undefined doublebonds",
            ),
        ]

        self.editor = MolEditWidget()
        self.editor.setStyleSheet("""
                                   QDialog {
                                    background-color: #1e1e2f;
                                    border-radius: 12px;
                                }
                                QLabel {
                                    color: #e0e0e0;
                                    font-size: 16px;
                                    padding: 12px;
                                }
                                QPushButton {
                                    background-color: #3a7ca5;
                                    color: #fff;
                                    border-radius: 7px;
                                    padding: 8px 22px;
                                    font-size: 14px;
                                }
                                QPushButton:hover {
                                    background-color: #28516b;
                                }
                                  """)
        self.chemEntityActionGroup = QtGui.QActionGroup(self, exclusive=True)
        self.ptable = PTable(self.chemEntityActionGroup)
        # Example for one button, apply to all atom buttons
        self.ptable.setStyleSheet("""
                QWidget {
                    background-color: #f4f6fb;
                    border-radius: 12px;
                }
                QToolButton {
                    background-color: #e3e6ee;
                    color: #23272e;
                    border: 2px solid #b0b0b0;
                    border-radius: 8px;
                    font-size: 14px;
                    font-weight: bold;
                    min-width: 24px;
                    min-height: 24px;
                    margin: 2px;
                    padding: 4px;
                }
                QToolButton:hover {
                    background-color: #b3c6f7;
                    color: #23272e;
                    border: 2px solid #6272a4;
                }
                QToolButton:checked, QToolButton:checked:!hover {
                    background-color: #6272a4;
                    color: #fff;
                    border: 2px solid #22223b;
                }
                QLabel {
                    color: #23272e;
                    font-size: 15px;
                    font-weight: bold;
                }
            """)
        icon1 = QtGui.QIcon()
        icon1.addPixmap(QtGui.QPixmap(appicon), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.ptable.setWindowIcon(QtGui.QIcon(icon1))
        def on_atomtype_changed(atomname):
            self.editor.setChemEntity(atomname)
            self.statusbar.showMessage(f"Atomtype {atomname} selected")
        self.ptable.atomtypeChanged.connect(on_atomtype_changed)

        self._fileName = None
        self.initGUI(fileName=fileName)
     
        # self.singleBondAction.trigger()
        self.ptable.atomActions["C"].trigger()

    # Properties
    @property
    def fileName(self):
        return self._fileName

    @fileName.setter
    def fileName(self, filename):
        if filename != self._fileName:
            self._fileName = filename
            

    def initGUI(self, fileName=None):
        self.setWindowTitle("Chem-Parse")
        self.setWindowIcon(QIcon.fromTheme(appicon))
        #self.setGeometry(100, 100, 200, 150)
        self.resize(800, 800)
        container = QtWidgets.QWidget()
        container.setObjectName("mainContainer")
        container.setMinimumWidth(700)
        container.setMaximumWidth(700)
        container.setStyleSheet("""
                QWidget#mainContainer {
                background: #ffffff;
                border: 5px solid #b0b0b0;
                border-radius: 30px;
            }
        """)
        layout = QtWidgets.QHBoxLayout(container)
        layout.setAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignCenter) 
        layout.setContentsMargins(20, 5, 0, 0) # Align left (and top if you want)
        self.editor.setMinimumSize(600, 600)
        self.editor.setMaximumSize(600, 600)
        
        layout.addWidget(self.editor)
        
        
        self.setCentralWidget(container)
        self.fileName = fileName

        self.filters = "MOL Files (*.mol *.mol);;SMILES Files (*.smi *.smi);;Any File (*)"
        self.filters_save = "MOL Files (*.mol *.mol);;SMILES Files (*.smi *.smi);;PNG File (*.png *.png);;SDF File (*.sdf *.sdf);;Any File (*)"
        self.filters_open = "MOL Files (*.mol *.mol);;SMILES Files (*.smi *.smi);;SDF File (*.sdf *.sdf);;Any File (*)"
        

        self.infobar = QtWidgets.QLabel("")

        if self.fileName is not None:
            self.editor.logger.info("Loading molecule from %s" % self.fileName)
           

        self.editor.sanitizeSignal.connect(self.infobar.setText)
        self.menubar = QtWidgets.QMenuBar(self)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1234, 26))
        self.menubar.setObjectName("menubar")
        self.menuFile = QtWidgets.QMenu(self.menubar)
        self.menuFile.setObjectName("menuFile")
        self.menuHelp = QtWidgets.QMenu(self.menubar)
        self.menuHelp.setObjectName("menuHelp")
        self.menuHelp_2 = QtWidgets.QMenu(self.menubar)
        self.menuHelp_2.setObjectName("menuHelp_2")
        self.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(self)
        self.statusbar.setObjectName("statusbar")
        self.setStatusBar(self.statusbar)
        self.logoDock = QtWidgets.QDockWidget(self)
        self.logoDock.setAllowedAreas(QtCore.Qt.RightDockWidgetArea)
        self.logoDock.setFeatures(QtWidgets.QDockWidget.DockWidgetMovable)
        self.logoDock.setTitleBarWidget(QtWidgets.QWidget())  # Hide the title bar
        self.logoDock.setMinimumWidth(0)  # Adjust width as needed
        self.logoDock.setObjectName("logoDock")

        mainWidget = QtWidgets.QWidget()
        mainWidget.setMinimumSize(0, 0)
        mainWidget.setStyleSheet(""" 
    QWidget {
        background: #ffffff;
        
    }
    QTableWidget {
        background: #ffffff;
        border: 1px solid #dcdde1;
        font-size: 14px;
        selection-background-color: #d0e6fa;
        selection-color: #222;
    }
     QTableCornerButton::section {
        background: #dcdde1;
        border: none;
    }
                                 
            QScrollBar:vertical {
        border: none;
        background: #f5f6fa;
        width: 12px;
        margin: 0px 0px 0px 0px;
        border-radius: 6px;
    }
    QScrollBar::handle:vertical {
        background: #b2bec3;
        min-height: 20px;
        border-radius: 6px;
    }
    QScrollBar::handle:vertical:hover {
        background: #4078c0;
    }
    QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {
        height: 0px;
        background: none;
        border: none;
    }
       QScrollBar:horizontal {
        border: none;
        background: #f5f6fa;
        height: 12px;
        margin: 0px 0px 0px 0px;
        border-radius: 6px;
    }
    QScrollBar::handle:horizontal {
        background: #b2bec3;
        min-width: 20px;
        border-radius: 6px;
    }
    QScrollBar::handle:horizontal:hover {
        background: #4078c0;
    }
    QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal {
        width: 0px;
        background: none;
        border: none;
    }
    QHeaderView::section {
        background: #dcdde1;
        font-weight: bold;
        padding: 6px;
        border: none;
        font-size: 14px;
    }
    QLineEdit {
        background: #ffffff;
        border: 1px solid #b2bec3;
        border-radius: 4px;
        padding: 6px;
        font-size: 14px;
        color: #222f3e;
        selection-background-color: #d0e6fa;
    }
    QLineEdit:focus {
        border: 1.5px solid #4078c0;
        background: #f0f6ff;
    }
    QPushButton {
        background-color: #4078c0;
        color: white;
        border-radius: 6px;
        padding: 8px 0;
        margin-bottom: 12px;
        font-size: 14px;
    }
    QPushButton:hover {
        background-color: #305080;
    }
    QComboBox {
        background: #ffffff;
        border: 1px solid #dcdde1;
        border-radius: 4px;
        padding: 6px;
        font-size: 14px;
    }
""")
        mainLayout = QtWidgets.QHBoxLayout(mainWidget)
        mainLayout.setContentsMargins(10, 10, 10, 10)

        # Table widget with headers
        self.tableWidget = QtWidgets.QTableWidget()
        self.tableWidget.setStyleSheet("""
               QScrollBar:vertical {
        border: none;
        background: #f5f6fa;
        width: 12px;
        margin: 0px 0px 0px 0px;
        border-radius: 6px;
    }
    QScrollBar::handle:vertical {
        background: #b2bec3;
        min-height: 20px;
        border-radius: 6px;
    }
    QScrollBar::handle:vertical:hover {
        background: #4078c0;
    }
    QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {
        height: 0px;
        background: none;
        border: none;
    }
       QScrollBar:horizontal {
        border: none;
        background: #f5f6fa;
        height: 12px;
        margin: 0px 0px 0px 0px;
        border-radius: 6px;
    }
    QScrollBar::handle:horizontal {
        background: #b2bec3;
        min-width: 20px;
        border-radius: 6px;
    }
    QScrollBar::handle:horizontal:hover {
        background: #4078c0;
    }
    QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal {
        width: 0px;
        background: none;
        border: none;
    }
                                       """)
        self.tableWidget.setColumnCount(3)
        self.tableWidget.setHorizontalHeaderLabels(["Name","SMILES", "Value"])
        self.tableWidget.horizontalHeader().setStretchLastSection(True)
        self.tableWidget.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.tableWidget.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.tableWidget.setMaximumWidth(300)
        self.tableWidget.setMinimumWidth(0)
        self.tableWidget.setEditTriggers(QtWidgets.QAbstractItemView.AllEditTriggers)

        mainLayout.addWidget(self.tableWidget, stretch=3)

        # Right-side controls
        rightWidget = QtWidgets.QWidget()
        rightWidget.setStyleSheet("""
                    QComboBox {
        background: #ffffff;
        border: 1px solid #b2bec3;
        border-radius: 4px;
        padding: 6px;
        font-size: 14px;
        color: #222f3e;
    }
    QComboBox:focus {
        border: 1.5px solid #4078c0;
        background: #f0f6ff;
    }
    QComboBox QAbstractItemView {
        background: #ffffff;
        selection-background-color: #d0e6fa;
        font-size: 14px;
    } QComboBox {
        background: #ffffff;
        border: 1px solid #b2bec3;
        border-radius: 4px;
        padding: 6px;
        font-size: 14px;
        color: #222f3e;
    }
    QComboBox:focus {
        border: 1.5px solid #4078c0;
        background: #f0f6ff;
    }
    QComboBox QAbstractItemView {
        background: #ffffff;
        selection-background-color: #d0e6fa;
        font-size: 14px;
    }
                                  """)
        rightWidget.setMaximumWidth(400)
        rightWidget.setMinimumWidth(0)
        rightWidget.setContentsMargins(10, 10, 60, 10)
        rightLayout = QtWidgets.QVBoxLayout(rightWidget)
        rightLayout.setAlignment(QtCore.Qt.AlignTop)
        def process_column(df, column_name, method):
                """
                Processes a column in the DataFrame based on the selected method.
                
                Parameters:
                    df (pd.DataFrame): The input DataFrame.
                    column_name (str): The column to process.
                    method (str): One of the keys from ylabel_processors.

                Returns:
                    pd.Series: The processed column.
                """
                ylabel_processors = {
                "None": lambda col: col,
                "pIC50(nm)": lambda col: -np.log10(col * 1e-9),   # Convert nM to M then to pIC50
                "pIC50(mm)": lambda col: -np.log10(col * 1e-3),   # Convert mM to M then to pIC50
                "LOG": lambda col: np.log10(col),
                "Min-Max": lambda col: (col - col.min()) / (col.max() - col.min()) if col.max() != col.min() else col * 0,
                "STD scalar": lambda col: (col - col.mean()) / col.std() if col.std() != 0 else col * 0
            }
                if method not in ylabel_processors:
                    raise ValueError(f"Unknown method: {method}")
                import pandas as pd
                df[column_name] = pd.to_numeric(df[column_name], errors='coerce')
                df[column_name] = ylabel_processors[method](df[column_name])

        def savetable(self):
            # Open save dialog with filters
          
            options = QtWidgets.QFileDialog.Options()
            filename, selected_filter = QtWidgets.QFileDialog.getSaveFileName(
                self,
                "Save Table",
                "",
                "Excel Files (*.xlsx);;CSV Files (*.csv);;SMILES Files (*.smi)",
                options=options
            )
            if not filename:
                return

            # Get table data
            row_count = self.tableWidget.rowCount()
            col_count = self.tableWidget.columnCount()
            headers = [self.tableWidget.horizontalHeaderItem(i).text() for i in range(col_count)]
            data = []
            for row in range(row_count):
                row_data = []
                for col in range(col_count):
                    item = self.tableWidget.item(row, col)
                    row_data.append(item.text() if item else "")
                data.append(row_data)

            # Save as Excel
            if filename.endswith('.xlsx'):
                try:
                    import pandas as pd
                    df = pd.DataFrame(data, columns=headers)
                    method = ylabelcombo.currentText()
                    process_column(df, "Value", method)
                    df.to_excel(filename, index=False)
                    StyledMessageDialog("Success", f"XLSX file saved to {filename}", self).exec()
                except Exception as e:
                    StyledMessageDialog("Error", "Pandas is required to save as Excel.", self).exec()
            # Save as CSV
            elif filename.endswith('.csv'):
                try:
                    import pandas as pd
                    df = pd.DataFrame(data, columns=headers)
                    method = ylabelcombo.currentText()
                    process_column(df, "Value", method)
                    df.to_csv(filename, index=False)
                    StyledMessageDialog("Success", f"CSV file saved to {filename}", self).exec()
                except ImportError:
                    StyledMessageDialog("Error", "Pandas is required to save as CSV.", self).exec()
                    
            # Save as SMILES
            elif filename.endswith('.smi'):
                # Check for SMILES and Name columns
                try:
                    smiles_idx = headers.index("SMILES")
                    name_idx = headers.index("Name")
                except ValueError:
                  
                    StyledMessageDialog("Error", "Table must have 'SMILES' and 'Name' columns to save as .smi.", self).exec()
                    return
                with open(filename, "w", encoding="utf-8") as f:
                    for row in data:
                        smiles = row[smiles_idx].strip()
                        name = row[name_idx].strip()
                        if smiles and name:
                            f.write(f"{smiles}\t{name}\n")
                    StyledMessageDialog("Success", f"SMILES file saved to {filename}", self ).exec()

        saveBtn = QtWidgets.QPushButton("Save")
        saveBtn.clicked.connect(lambda: savetable(self))    
        
        def cleartable():
           try:
            self.tableWidget.clearContents()
           except Exception as e:
               print(str(e))

        clearBtn = QtWidgets.QPushButton("Clear")
        clearBtn.clicked.connect(cleartable)
        def OPT(self):
         
            from rdkit import Chem
            from rdkit.Chem import AllChem
            if os.path.isdir(output_dir):
                for file in os.listdir(output_dir):
                    os.remove(os.path.join(output_dir, file))

            options = QtWidgets.QFileDialog.Options()
            filename = QtWidgets.QFileDialog.getExistingDirectory(
                self,
                "Save 3D Structures",
                "",
                options=options
            )
            if not filename:
                return
            headers = [self.tableWidget.horizontalHeaderItem(i).text() for i in range(self.tableWidget.columnCount())]
            try:
                smiles_idx = headers.index("SMILES")
                name_idx = headers.index("Name")
            except ValueError:
             
                StyledMessageDialog("Error", "Table must have 'SMILES' , 'Name' columns to save as SDF.", self).exec()
                return
           
            for row in range(self.tableWidget.rowCount()):
                smiles_item = self.tableWidget.item(row, smiles_idx)
                name_item = self.tableWidget.item(row, name_idx)
                if smiles_item and name_item:
                    smiles = smiles_item.text().strip()
                    name = name_item.text().strip()
                    text = combofile.currentText()
                    dirx = os.path.dirname(output_dir)
                    fullname =  name + ".sdf"
                    fullpath = os.path.join(dirx, fullname)
                    writer = Chem.SDWriter(fullpath)
                        

                    if smiles:
                        mol = Chem.MolFromSmiles(smiles)
                        if mol:
                            mol = Chem.AddHs(mol)
                            try:
                               
                                AllChem.EmbedMolecule(mol, AllChem.ETKDG())
                                AllChem.UFFOptimizeMolecule(mol)
                                writer.write(mol)
                            except Exception as e:
                                self.editor.logger.warning(f"3D generation failed for {name}: {e}")
                                StyledMessageDialog("Error", f"3D generation failed for {name}: {e}", self).exec()
            writer.close()
            dirx = os.path.dirname(output_dir) 
            os.chdir(dirx)
            if combofile.currentText() == 'MOL2':
                ff = ffcombo.currentText()
                try:
                 prot = prot_line.text()
                 stepsx = stepsline.text()
                 protfloat = float(prot)
                 stepsxy = int(stepsx)
                 if protfloat < 1 :
                     return StyledMessageDialog("Error", "Protonation state should be in range 1-14", self).exec()
                 elif protfloat > 14:
                     return StyledMessageDialog("Error", "Protonation state should be in range 1-14", self).exec()
                 dialog = BatchProcessMOL2Dialog(mol2batch, ff, protfloat, stepsxy)
                 dialog.exec()
                 StyledMessageDialog("Chem-Parser", "Successfully Converted to MOL2", self).exec()
                 files = os.listdir(output_dir)
                 if files:
                     for file in files:
                         path = os.path.join(output_dir, file)
                         shutil.move(path, filename)
                     StyledMessageDialog("Chem-Parser", "Successfully Saved MOL2", self).exec()

                except Exception as e:
                    StyledMessageDialog("Error", "Give Proper Value Type ")
            elif combofile.currentText() == 'PDBQT':

                ff = ffcombo.currentText()
                try:
                 prot = prot_line.text()
                 stepsx = stepsline.text()
                 protfloat = float(prot)
                 stepsxy = int(stepsx)
                 if protfloat < 1 :
                     return StyledMessageDialog("Error", "Protonation state should be in range 1-14", self).exec()
                 elif protfloat > 14:
                     return StyledMessageDialog("Error", "Protonation state should be in range 1-14", self).exec()
                 dialog = BatchProcessPDBQTDialog(pdbqtbatch, ff, protfloat, stepsxy)
                 dialog.exec()
                 StyledMessageDialog("Chem-Parser", "Successfully Converted to PDBQT", self).exec()
                 files = os.listdir(output_dir)
                 if files:
                     for file in files:
                         path = os.path.join(output_dir, file)
                         shutil.move(path, filename)
                     StyledMessageDialog("Chem-Parser", "Successfully Saved PDBQT", self).exec()
                 
                except Exception as e:
                    StyledMessageDialog("Error", "Give Proper Value Type ")


        def save3D(self):
            from rdkit import Chem
            from rdkit.Chem import AllChem

            row_count = self.tableWidget.rowCount()
            col_count = self.tableWidget.columnCount()
            headers = [self.tableWidget.horizontalHeaderItem(i).text() for i in range(col_count)]
            data = []
            for row in range(row_count):
                row_data = []
                for col in range(col_count):
                    item = self.tableWidget.item(row, col)
                    row_data.append(item.text() if item else "")
                data.append(row_data)
            import pandas as pd
            df = pd.DataFrame(data, columns=headers)
            method = ylabelcombo.currentText()
            if "Value" in df.columns:
                ylabel_processors = {
                    "None": lambda col: col,
                    "pIC50(nm)": lambda col: -np.log10(col * 1e-9),   # Convert nM to M then to pIC50
                    "pIC50(mm)": lambda col: -np.log10(col * 1e-3),   # Convert mM to M then to pIC50
                    "LOG": lambda col: np.log10(col),
                    "Min-Max": lambda col: (col - col.min()) / (col.max() - col.min()) if col.max() != col.min() else col * 0,
                    "STD scalar": lambda col: (col - col.mean()) / col.std() if col.std() != 0 else col * 0
                }
                try:
                    import pandas as pd
                    df["Value"] = pd.to_numeric(df["Value"], errors='coerce')
                    df["Value"] = ylabel_processors[method](df["Value"])
                except Exception as e:
                    StyledMessageDialog("Error", f"Failed to process 'Value' column: {e}", self).exec()
                    return
            options = QtWidgets.QFileDialog.Options()
            filename, _ = QtWidgets.QFileDialog.getSaveFileName(
                self,
                "Save 3D Structures",
                "",
                "SDF Files (*.sdf)",
                options=options
            )
            if not filename:
                return

            # Validate necessary columns
            try:
                smiles_idx = headers.index("SMILES")
                name_idx = headers.index("Name")
            except ValueError:
                StyledMessageDialog("Error", "Table must have 'SMILES' and 'Name' columns.", self).exec()
                return

            writer = Chem.SDWriter(filename)
            for row in range(row_count):
                smiles = df.at[row, "SMILES"]
                name = df.at[row, "Name"]
                value = df.at[row, "Value"] if "Value" in df.columns else None

                if smiles and name:
                    mol = Chem.MolFromSmiles(smiles)
                    if mol:
                        mol = Chem.AddHs(mol)
                        try:
                            steps = stepsline.text()
                            int_steps = int(steps)
                            AllChem.EmbedMolecule(mol, AllChem.ETKDG())
                            AllChem.UFFOptimizeMolecule(mol, maxIters=int_steps)
                        except Exception as e:
                            StyledMessageDialog("Chem-Parse", "Steps should be an integer", self).exec()
                            continue

                        mol.SetProp("_Name", name)
                        if value is not None and pd.notna(value):
                            mol.SetProp("_Value", str(value))

                        writer.write(mol)

            writer.close()
            StyledMessageDialog("Success", f"3D structures saved to {filename}", self).exec()
            # Open save dialog for SDF
            
        def choiceformat(self):
            try:
                text = combofile.currentText()
                if text == 'SDF':
                    save3D(self)
                else:
                    OPT(self)
            except Exception as e:
                StyledMessageDialog("Error", str(e), self).exec

        save3DBtn = QtWidgets.QPushButton("3D")
        save3DBtn.clicked.connect(lambda: choiceformat(self))
        filelabel = QtWidgets.QLabel("3D Format:")
        combofile = QtWidgets.QComboBox()
        combofile.addItem("SDF")
        combofile.addItem("MOL2")
        combofile.addItem("PDBQT")
        fflabel = QtWidgets.QLabel("Forcefield:")
        ffcombo = QtWidgets.QComboBox()
        ffcombo.addItem("MMFF94")
        ffcombo.addItem("MMFF94s")
        ffcombo.addItem("GAFF")
        ffcombo.addItem("Ghemical")
        ffcombo.addItem("UFF")
        prot_label = QtWidgets.QLabel("Protonation:")
        prot_line = QtWidgets.QLineEdit()
        prot_line.setText("7.4")
        stepslabel = QtWidgets.QLabel("Steps:")
        stepsline = QtWidgets.QLineEdit()
        stepsline.setText("2500")
        y_label = QtWidgets.QLabel("Y-Label Processing:")
        ylabelcombo = QtWidgets.QComboBox()
        ylabelcombo.addItem("None")
        ylabelcombo.addItem("pIC50(nm)")
        ylabelcombo.addItem("pIC50(mm)")
        ylabelcombo.addItem("LOG")
        ylabelcombo.addItem("Min-Max")
        ylabelcombo.addItem("STD scalar")


        rightLayout.addWidget(saveBtn)
        rightLayout.addWidget(clearBtn)
        rightLayout.addWidget(save3DBtn)
        rightLayout.addWidget(filelabel)
        rightLayout.addWidget(combofile)
        rightLayout.addWidget(fflabel)
        rightLayout.addWidget(ffcombo)
        rightLayout.addWidget(prot_label)
        rightLayout.addWidget(prot_line)
        rightLayout.addWidget(stepslabel)
        rightLayout.addWidget(stepsline)
        rightLayout.addWidget(y_label)
        rightLayout.addWidget(ylabelcombo)

        mainLayout.addWidget(rightWidget, stretch=3)

        self.logoDock.setWidget(mainWidget)
        self.addDockWidget(QtCore.Qt.RightDockWidgetArea, self.logoDock)
        self.dockWidget = QtWidgets.QDockWidget(self)
        self.dockWidget.setFeatures(QtWidgets.QDockWidget.NoDockWidgetFeatures)
        self.dockWidget.setAllowedAreas(QtCore.Qt.LeftDockWidgetArea) 
        self.dockWidget.setTitleBarWidget(QtWidgets.QWidget())
        self.dockWidget.setMinimumWidth(100)
        self.dockWidget.setMaximumWidth(100)
                                        # Only allow left docking
        self.dockWidget.setObjectName("dockWidget")
        self.dockWidget.setStyleSheet("""
        QDockWidget {
          background: #f8f8f2;
          border: 3px solid #b0b0b0;   /* Gray border */
          border-radius: 8px;
        }
        
        QWidget {
          background: #f8f8f2;
          border: 3px solid #b0b0b0;   /* Gray border */
          border-radius: 8px;
        }
        QPushButton {
          background-color: #e0e0e0;
          color: #282a36;
          border: 1px solid #cccccc;
          border-radius: 6px;
          font-size: 13px;
        }
        QPushButton:hover {
                background-color: #d6d6d6;
                border: 1.5px solid #6272a4;
                color: #22223b;
            }
        """)
        
        self.dockWidgetContents = QtWidgets.QWidget()
        self.dockWidgetContents.setObjectName("dockWidgetContents")
        def setAction(TYPE):

            self.editor.setAction(TYPE)
            self.statusbar.showMessage("Action %s selected" % TYPE)

        def setAtomType(TYPE):
            
            self.editor.setChemEntity(TYPE)
            # self.editor.setRingType(None)
            self.statusbar.showMessage("Atomtype %s selected" % TYPE)

        def setBondType(TYPE):
            self.editor.setChemEntity(TYPE)
            self.statusbar.showMessage("Bondtype %s selected" % TYPE)

        def setRingType(TYPE):
            
            self.editor.setChemEntity(TYPE)
            self.statusbar.showMessage("Bondtype %s selected" % TYPE)
        
        self.pushButton_11 = QtWidgets.QPushButton(self.dockWidgetContents)
        self.pushButton_11.setGeometry(QtCore.QRect(10, 30, 30, 30))
        self.pushButton_11.setText("")
        icon1 = QtGui.QIcon()
        icon1.addPixmap(QtGui.QPixmap(singlebondpng), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.pushButton_11.setIcon(icon1)
        self.pushButton_11.setIconSize(QtCore.QSize(20, 20))
        self.pushButton_11.setObjectName("pushButton_11")
        self.pushButton_11.clicked.connect(lambda: setBondType("SINGLE"))
        def common_x(checked, type):
             if checked:
                 setRingType(type)
             else:
                 setBondType('SINGLE')

        self.pushButton_11_r = QtWidgets.QPushButton(self.dockWidgetContents)
        self.pushButton_11_r.setGeometry(QtCore.QRect(50, 30, 30, 30))
        self.pushButton_11_r.setText("")
        self.pushButton_11_r.setStatusTip("Benzene")
        icon1 = QtGui.QIcon()
        icon1.addPixmap(QtGui.QPixmap(benzene), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.pushButton_11_r.setIcon(icon1)
        self.pushButton_11_r.setIconSize(QtCore.QSize(20, 20))
        self.pushButton_11_r.setObjectName("pushButton_11_r")
        self.pushButton_11_r.setCheckable(True)
        self.pushButton_11_r.setStyleSheet("""
                                      QPushButton:checked {
                                        background-color: #6272a4; /* Change to your preferred checked color */
                                        color: #fff;
                                        border: 2px solid #22223b;
                                    }
                                           """)
        self.pushButton_11_r.toggled.connect(lambda checked: common_x(checked, "benzene"))
       
        self.pushButton_12 = QtWidgets.QPushButton(self.dockWidgetContents)
        self.pushButton_12.setGeometry(QtCore.QRect(10, 70, 30, 30))
        self.pushButton_12.setText("")
        icon2 = QtGui.QIcon()
        icon2.addPixmap(QtGui.QPixmap(doublebondpng), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.pushButton_12.setIcon(icon2)
        self.pushButton_12.setIconSize(QtCore.QSize(20, 20))
        self.pushButton_12.setObjectName("pushButton_12")
        self.pushButton_12.clicked.connect(lambda: setBondType("DOUBLE"))
       
                
        self.pushButton_12_r = QtWidgets.QPushButton(self.dockWidgetContents)
        self.pushButton_12_r.setGeometry(QtCore.QRect(50, 70, 30, 30))
        self.pushButton_12_r.setText("")
        icon2 = QtGui.QIcon()
        icon2.addPixmap(QtGui.QPixmap(cyclohexane), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.pushButton_12_r.setIcon(icon2)
        self.pushButton_12_r.setIconSize(QtCore.QSize(20, 20))
        self.pushButton_12_r.setObjectName("pushButton_12_r")
        self.pushButton_12_r.setCheckable(True)
        self.pushButton_12_r.setStyleSheet("""
                                      QPushButton:checked {
                                        background-color: #6272a4; /* Change to your preferred checked color */
                                        color: #fff;
                                        border: 2px solid #22223b;
                                    }
                                           """)
        self.pushButton_12_r.setStatusTip("Cyclohexane")
        self.pushButton_12_r.toggled.connect(lambda checked: common_x(checked, "cyclohexane"))
        self.pushButton_13 = QtWidgets.QPushButton(self.dockWidgetContents)
        self.pushButton_13.setGeometry(QtCore.QRect(10, 110, 30, 30))
        self.pushButton_13.setText("")
        icon3 = QtGui.QIcon()
        icon3.addPixmap(QtGui.QPixmap(triplebondpng), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.pushButton_13.setIcon(icon3)
        self.pushButton_13.setIconSize(QtCore.QSize(20, 20))
        self.pushButton_13.setObjectName("pushButton_13")
        self.pushButton_13.clicked.connect(lambda: setBondType("TRIPLE"))
        self.pushButton_13_r = QtWidgets.QPushButton(self.dockWidgetContents)
        self.pushButton_13_r.setGeometry(QtCore.QRect(50, 110, 30, 30))
        self.pushButton_13_r.setText("")
        icon3 = QtGui.QIcon()
        icon3.addPixmap(QtGui.QPixmap(cyclopentane), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.pushButton_13_r.setIcon(icon3)
        self.pushButton_13_r.setCheckable(True)
        self.pushButton_13_r.setStyleSheet("""
                                          QPushButton:checked {
                                        background-color: #6272a4; /* Change to your preferred checked color */
                                        color: #fff;
                                        border: 2px solid #22223b;
                                    }
                                           """)
        self.pushButton_13_r.setIconSize(QtCore.QSize(20, 20))
        self.pushButton_13_r.setObjectName("pushButton_13")
        self.pushButton_13_r.setStatusTip("Cyclopentane")
        self.pushButton_13_r.toggled.connect(lambda checked: common_x(checked, "cyclopentane"))
        def setAtomTypeName(atomname):
            self.editor.setChemEntity(str(atomname))
            self.statusbar.showMessage("Atomtype %s selected" % atomname)
        def checkforatoms(checked, atomname, button):
            if checked:
                button.clicked.connect(lambda: setAtomTypeName(atomname))
                self.statusbar.showMessage(f"{atomname} selected")
            else:
                button.clicked.connect(lambda: setBondType("SINGLE"))

        self.pushButton_R = QtWidgets.QPushButton(self.dockWidgetContents)
        self.pushButton_R.setGeometry(QtCore.QRect(10, 150, 30, 30))
        font = QtGui.QFont()
        font.setPointSize(12)
        font.setBold(True)
        font.setWeight(QtGui.QFont.Weight.Bold)
        self.pushButton_R.setFont(font)
        self.pushButton_R.setIconSize(QtCore.QSize(20, 20))
        self.pushButton_R.setObjectName("pushButton_R")
        self.pushButton_R.setCheckable(True)
      
        self.pushButton_R.setStyleSheet("""
                                          QPushButton:checked {
                                        background-color: #6272a4; /* Change to your preferred checked color */
                                        color: #fff;
                                        border: 2px solid #22223b;
                                    }
                                        """)
        self.pushButton_R.toggled.connect(
            lambda checked: checkforatoms(checked, "R", self.pushButton_R)
        )
       
        self.pushButton_R_r = QtWidgets.QPushButton(self.dockWidgetContents)
        self.pushButton_R_r.setGeometry(QtCore.QRect(50, 150, 30, 30))
        icon3 = QtGui.QIcon()
        icon3.addPixmap(QtGui.QPixmap(cyclobutane), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.pushButton_R_r.setIcon(icon3)
        self.pushButton_R_r.setIconSize(QtCore.QSize(20, 20))
        self.pushButton_R_r.setObjectName("pushButton_R")
        self.pushButton_R_r.setStatusTip("Cyclobutane")
        self.pushButton_R_r.setCheckable(True)
        self.pushButton_R_r.setStyleSheet("""
                                      QPushButton:checked {
                                        background-color: #6272a4; /* Change to your preferred checked color */
                                        color: #fff;
                                        border: 2px solid #22223b;
                                    }
                                          """)

        self.pushButton_R_r.toggled.connect(lambda checked: common_x(checked, "cyclobutane"))
        self.pushButton_R_2 = QtWidgets.QPushButton(self.dockWidgetContents)
        self.pushButton_R_2.setGeometry(QtCore.QRect(10, 190, 30, 30))
        font = QtGui.QFont()
        font.setPointSize(12)
        font.setBold(True)
        font.setWeight(QtGui.QFont.Weight.Bold)
        self.pushButton_R_2.setFont(font)
        self.pushButton_R_2.setStyleSheet("QPushButton { color: gray; }")
        self.pushButton_R_2.setIconSize(QtCore.QSize(20, 20))
        self.pushButton_R_2.setObjectName("pushButton_R_2")
        self.pushButton_R_2.setCheckable(True)
        self.pushButton_R_2_r = QtWidgets.QPushButton(self.dockWidgetContents)
        self.pushButton_R_2_r.setGeometry(QtCore.QRect(50, 190, 30, 30))
        
        self.pushButton_R_2_r.setFont(font)
        self.pushButton_R_2_r.setStyleSheet("""
                                  QPushButton:checked {
                                        background-color: #6272a4; /* Change to your preferred checked color */
                                        color: #fff;
                                        border: 2px solid #22223b;
                                    }
                                            """)
        icon3 = QtGui.QIcon()
        icon3.addPixmap(QtGui.QPixmap(cyclopropane), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.pushButton_R_2_r.setIcon(icon3)
        self.pushButton_R_2_r.setIconSize(QtCore.QSize(20, 20))
        self.pushButton_R_2_r.setObjectName("pushButton_R_2")
        self.pushButton_R_2_r.setCheckable(True)
        self.pushButton_R_2_r.toggled.connect(lambda checked: common_x(checked, "cyclopropane"))
        self.pushButton_R_2_r.setStatusTip("Cyclopropane")   
        self.pushButton_R_2.setStyleSheet("""
                                              QPushButton:checked {
                                        background-color: #6272a4; /* Change to your preferred checked color */
                                        color: #fff;
                                        border: 2px solid #22223b;
                                    }
                                     
                                          """)
        self.pushButton_R_2.toggled.connect(
            lambda checked: checkforatoms(checked, "C", self.pushButton_R_2)
        )
        self.pushButton_R_3 = QtWidgets.QPushButton(self.dockWidgetContents)
        self.pushButton_R_3.setGeometry(QtCore.QRect(10, 230, 30, 30))
        font = QtGui.QFont()
        font.setPointSize(12)
        font.setBold(False)
        font.setWeight(QtGui.QFont.Weight.Normal)
        self.pushButton_R_3.setFont(font)
        self.pushButton_R_3.setStyleSheet("QPushButton rgb(130, 130, 130)")
        self.pushButton_R_3.setIconSize(QtCore.QSize(20, 20))
        self.pushButton_R_3.setObjectName("pushButton_R_3")
        self.pushButton_R_3.setCheckable(True)
       
        self.pushButton_R_3.setStyleSheet("""
                                            QPushButton:checked {
                                        background-color: #6272a4; /* Change to your preferred checked color */
                                        color: #fff;
                                        border: 2px solid #22223b;
                                    }
                                     
                
                                          """)
        self.pushButton_R_3.toggled.connect(
            lambda checked: checkforatoms(checked, "H", self.pushButton_R_3)
        )
        self.pushButton_R_3_r = QtWidgets.QPushButton(self.dockWidgetContents)
        self.pushButton_R_3_r.setGeometry(QtCore.QRect(50, 230, 30, 30))
        font = QtGui.QFont()
        font.setPointSize(12)
        font.setBold(True)
        font.setWeight(QtGui.QFont.Weight.Bold)
        self.pushButton_R_3_r.setFont(font)
        self.pushButton_R_3_r.setText("COOH")
        self.pushButton_R_3_r.setStyleSheet("""
                                            
                                              QPushButton { 
                                          color: red;
                                          font-size: 8px;
                                     
                                          }
                                            QPushButton:checked {
                                        background-color: #6272a4; /* Change to your preferred checked color */
                                        color: #fff;
                                        border: 2px solid #22223b;
                                    }
                
                                          """)

        self.pushButton_R_3_r.setObjectName("pushButton_R_3")
        self.pushButton_R_3_r.setCheckable(True)
        self.pushButton_R_3_r.toggled.connect(lambda checked: common_x(checked, "carboxylic acid"))
       
     
        self.pushButton_R_4 = QtWidgets.QPushButton(self.dockWidgetContents)
        self.pushButton_R_4.setGeometry(QtCore.QRect(10, 270, 30, 30))
        font = QtGui.QFont()
        font.setPointSize(12)
        font.setBold(True)
        font.setWeight(QtGui.QFont.Weight.Bold)
        self.pushButton_R_4.setFont(font)
        self.pushButton_R_4.setStyleSheet("QPushButton { color: red; }")
        self.pushButton_R_4.setCheckable(True)
       
        self.pushButton_R_4.setStyleSheet("""
                                          QPushButton { 
                                          color: red; 
                                          }
                                            QPushButton:checked {
                                        background-color: #6272a4; /* Change to your preferred checked color */
                                        color: #fff;
                                        border: 2px solid #22223b;
                                    }
                                     
                
                                          """)
        self.pushButton_R_4.toggled.connect(
            lambda checked: checkforatoms(checked, "O", self.pushButton_R_4)
        )
        self.pushButton_R_4.setIconSize(QtCore.QSize(20, 20))
        self.pushButton_R_4.setObjectName("pushButton_R_4")
        self.pushButton_R_5 = QtWidgets.QPushButton(self.dockWidgetContents)
        self.pushButton_R_5.setGeometry(QtCore.QRect(10, 310, 30, 30))
        font = QtGui.QFont()
        font.setPointSize(12)
        font.setBold(True)
        font.setWeight(QtGui.QFont.Weight.Bold)
        self.pushButton_R_5.setFont(font)
        self.pushButton_R_5.setIconSize(QtCore.QSize(20, 20))
        self.pushButton_R_5.setObjectName("pushButton_R_5")
        self.pushButton_R_5.setCheckable(True)
      
        self.pushButton_R_5.setStyleSheet("""
                                          QPushButton { 
                                          color: blue; 
                                          }
                                            QPushButton:checked {
                                        background-color: #6272a4; /* Change to your preferred checked color */
                                        color: #fff;
                                        border: 2px solid #22223b;
                                    }
                                     
                
                                          """)
        self.pushButton_R_5.toggled.connect(
            lambda checked: checkforatoms(checked, "N", self.pushButton_R_5)
        )
        self.pushButton_R_6 = QtWidgets.QPushButton(self.dockWidgetContents)
        self.pushButton_R_6.setGeometry(QtCore.QRect(10, 350, 30, 30))
        font = QtGui.QFont()
        font.setPointSize(12)
        font.setBold(True)
        font.setWeight(QtGui.QFont.Weight.Bold)
        self.pushButton_R_6.setFont(font)
        self.pushButton_R_6.setIconSize(QtCore.QSize(20, 20))
        self.pushButton_R_6.setObjectName("pushButton_R_6")
        self.pushButton_R_6.setCheckable(True)
    
        self.pushButton_R_6.setStyleSheet("""
                                          QPushButton { 
                                          color: yellow; 
                                          }
                                            QPushButton:checked {
                                        background-color: #6272a4; /* Change to your preferred checked color */
                                        color: #fff;
                                        border: 2px solid #22223b;
                                    }
                                     
                
                                          """)
        self.pushButton_R_6.toggled.connect(
            lambda checked: checkforatoms(checked, "S", self.pushButton_R_6)
        )

        self.pushButton_R_7 = QtWidgets.QPushButton(self.dockWidgetContents)
        self.pushButton_R_7.setGeometry(QtCore.QRect(10, 390, 30, 30))
        font = QtGui.QFont()
        font.setPointSize(12)
        font.setBold(False)
        font.setWeight(QtGui.QFont.Weight.Normal)
        self.pushButton_R_7.setFont(font)
        self.pushButton_R_7.setStyleSheet("QPushButton rgb(130, 130, 130)")
        self.pushButton_R_7.setIconSize(QtCore.QSize(20, 20))
        self.pushButton_R_7.setObjectName("pushButton_R_7")
        self.pushButton_R_7.setCheckable(True)
       
        self.pushButton_R_7.setStyleSheet("""
                                         
                                            QPushButton:checked {
                                        background-color: #6272a4; /* Change to your preferred checked color */
                                        color: #fff;
                                        border: 2px solid #22223b;
                                    }
                                     
                
                                          """)
        self.pushButton_R_7.toggled.connect(
            lambda checked: checkforatoms(checked, "Cl", self.pushButton_R_7)
        )
        
        self.pushButton_R_8 = QtWidgets.QPushButton(self.dockWidgetContents)
        self.pushButton_R_8.setGeometry(QtCore.QRect(10, 430, 30, 30))
        font = QtGui.QFont()
        font.setPointSize(12)
        font.setBold(False)
        font.setWeight(QtGui.QFont.Weight.Normal)
        self.pushButton_R_8.setFont(font)
        self.pushButton_R_8.setStyleSheet("QPushButton rgb(130, 130, 130)")
        self.pushButton_R_8.setIconSize(QtCore.QSize(20, 20))
        self.pushButton_R_8.setObjectName("pushButton_R_8")
        self.pushButton_R_8.setCheckable(True)
        
        self.pushButton_R_8.setStyleSheet("""
                                         
                                            QPushButton:checked {
                                        background-color: #6272a4; /* Change to your preferred checked color */
                                        color: #fff;
                                        border: 2px solid #22223b;
                                    }
                                     
                
                                          """)
        self.pushButton_R_8.toggled.connect(
            lambda checked: checkforatoms(checked, "F", self.pushButton_R_8)
        )
        self.pushButton_R_9 = QtWidgets.QPushButton(self.dockWidgetContents)
        self.pushButton_R_9.setGeometry(QtCore.QRect(10, 470, 30, 30))
        font = QtGui.QFont()
        font.setPointSize(12)
        font.setBold(False)
        font.setWeight(QtGui.QFont.Weight.Normal)
        self.pushButton_R_9.setFont(font)
        self.pushButton_R_9.setStyleSheet("QPushButton rgb(130, 130, 130)")
        self.pushButton_R_9.setIconSize(QtCore.QSize(20, 20))
        self.pushButton_R_9.setObjectName("pushButton_R_9")
        self.pushButton_R_9.setCheckable(True)
       
        self.pushButton_R_9.setStyleSheet("""
                                         
                                            QPushButton:checked {
                                        background-color: #6272a4; /* Change to your preferred checked color */
                                        color: #fff;
                                        border: 2px solid #22223b;
                                    }
                                     
                
                                          """)
        self.pushButton_R_9.toggled.connect(
            lambda checked: checkforatoms(checked, "Br", self.pushButton_R_9)
        )
        self.pushButton_R_10 = QtWidgets.QPushButton(self.dockWidgetContents)
        self.pushButton_R_10.setGeometry(QtCore.QRect(10, 510, 30, 30))
        font = QtGui.QFont()
        font.setPointSize(12)
        font.setBold(False)
        font.setWeight(QtGui.QFont.Weight.Normal)
        self.pushButton_R_10.setFont(font)
        self.pushButton_R_10.setStyleSheet("QPushButton rgb(130, 130, 130)")
        self.pushButton_R_10.setIconSize(QtCore.QSize(20, 20))
        self.pushButton_R_10.setObjectName("pushButton_R_10")
        self.pushButton_R_10.setCheckable(True)
       
        self.pushButton_R_10.setStyleSheet("""
                                         
                                            QPushButton:checked {
                                        background-color: #6272a4; /* Change to your preferred checked color */
                                        color: #fff;
                                        border: 2px solid #22223b;
                                    }
                                     
                
                                          """)
        self.pushButton_R_10.toggled.connect(
            lambda checked: checkforatoms(checked, "P", self.pushButton_R_10)
        )
        self.pushButton_R_11 = QtWidgets.QPushButton(self.dockWidgetContents)
        self.pushButton_R_11.setGeometry(QtCore.QRect(10, 550, 30, 30))
        font = QtGui.QFont()
        font.setPointSize(12)
        font.setBold(False)
        font.setWeight(QtGui.QFont.Weight.Normal)
        self.pushButton_R_11.setFont(font)
        self.pushButton_R_11.setStyleSheet("QPushButton rgb(130, 130, 130)")
        self.pushButton_R_11.setIconSize(QtCore.QSize(20, 20))
        self.pushButton_R_11.setObjectName("pushButton_R_11")
        self.pushButton_R_12 = QtWidgets.QPushButton(self.dockWidgetContents)
        self.pushButton_R_12.setGeometry(QtCore.QRect(10, 550, 30, 30))
        font = QtGui.QFont()
        font.setPointSize(12)
        font.setBold(False)
        font.setWeight(QtGui.QFont.Weight.Normal)
        def openPtable():
          self.ptable.show()
        self.pushButton_R_12.setFont(font)
        self.pushButton_R_12.setStyleSheet("QPushButton rgb(130, 130, 130)")
        self.pushButton_R_12.setText("")
        icon4 = QtGui.QIcon()
        icon4.addPixmap(QtGui.QPixmap(tablepng), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.pushButton_R_12.setIcon(icon4)
        self.pushButton_R_12.setIconSize(QtCore.QSize(20, 20))
        self.pushButton_R_12.setObjectName("pushButton_R_12")
        self.pushButton_R_12.clicked.connect(openPtable)
     
        self.dockWidget.setWidget(self.dockWidgetContents)
        self.addDockWidget(QtCore.Qt.LeftDockWidgetArea, self.dockWidget)

        # Bottom DockWidget
        self.dockWidget_2 = QtWidgets.QDockWidget(self)
        self.dockWidget_2.setAllowedAreas(QtCore.Qt.TopDockWidgetArea)
        self.dockWidget_2.setFeatures(QtWidgets.QDockWidget.NoDockWidgetFeatures)
        self.dockWidget_2.setTitleBarWidget(QtWidgets.QWidget())
        self.dockWidget_2.setMinimumHeight(80)
        self.dockWidget_2.setMinimumWidth(800)
        self.dockWidget_2.setStyleSheet("""
            QDockWidget {
                background: #f8f8f2;
                border: 2px solid #b0b0b0;
                border-radius: 10px;
                padding: 4px;
            }
            QWidget {
                background: #f8f8f2;
                border: 1px solid #b0b0b0;
                border-radius: 10px;
                padding: 1px;
            }
            
            QPushButton:hover {
                background-color: #d6d6d6;
                border: 1.5px solid #6272a4;
                color: #22223b;
            }
        """)
        self.dockWidgetContents_2 = QtWidgets.QWidget()
        self.dockWidgetContents_2.setObjectName("dockWidgetContents")

        def loadSmilesFile(filename):
            self.fileName = filename
            with open(self.fileName, "r") as file:
                lines = file.readlines()
                if len(lines) > 1:
                    self.editor.logger.warning("The SMILES file contains more than one line.")
                    self.statusBar().showMessage("The SMILES file contains more than one line.")
                    StyledMessageDialog("Chem-Parse", "The SMILES file contains more than one line."  , self ).exec()
                    return None
                smiles = lines[0].strip()
                mol = Chem.MolFromSmiles(smiles)
                self.editor.mol = mol
                self.statusBar().showMessage(f"SMILES file {filename} opened")
                StyledMessageDialog("Chem-Parse", f"SMILES file {filename} opened" , self ).exec()

        def loadMolFile(filename):
            self.fileName = filename
            mol = Chem.MolFromMolFile(str(self.fileName), sanitize=False, strictParsing=False)
            self.editor.mol = mol
            self.statusBar().showMessage(f"Mol file {filename} opened")
            StyledMessageDialog("Chem-Parse", f"Mol file {filename} opened", self ).exec()

        def loadSDFFile(filename):
            self.fileName = filename
            supplier = Chem.SDMolSupplier(str(self.fileName), sanitize=False)
            
            # Get the first valid molecule
            mol = None
            for m in supplier:
                if m is not None:
                    mol = m
                    break

            if mol:
                self.editor.mol = mol
                self.statusBar().showMessage(f"SDF file {filename} opened", 2000)
                StyledMessageDialog("Chem-Parse", f"SDF file {filename} opened", self).exec()
            else:
                self.statusBar().showMessage("No valid molecule found in SDF", 2000)
                StyledMessageDialog("Chem-Parse", "No valid molecule found in SDF", self).exec()

        def openFile():
            self.fileName, _ = QFileDialog.getOpenFileName(self, caption="Open file", filter=self.filters_open)
            return loadFile()

        def loadFile():
            if not self.fileName:
                self.editor.logger.warning("No file selected.")
                self.statusBar().showMessage("No file selected.")
                StyledMessageDialog("Chem-Parse", "No file selected" , self ).exec()
                return
            if self.fileName.lower().endswith(".mol"):
                loadMolFile(self.fileName)
            elif self.fileName.lower().endswith(".smi"):
                loadSmilesFile(self.fileName)
            elif self.fileName.lower().endswith(".sdf"):
                loadSDFFile(self.fileName)
            else:
                self.editor.logger.warning("Unknown file format. Assuming file as .mol format.")
                self.statusBar().showMessage("Unknown file format. Assuming file as .mol format.")
                StyledMessageDialog("Chem-Parse", "Unknown file format. Assuming file as .mol format." , self ).exec()
                loadMolFile(self.fileName)
                self.fileName += ".mol"
      
        self.pushButton = QtWidgets.QPushButton(self.dockWidgetContents_2)
        self.pushButton.setGeometry(QtCore.QRect(10, 10, 51, 51))
        self.pushButton.setText("")
        icon5 = QtGui.QIcon()
        icon5.addPixmap(QtGui.QPixmap(fileopen), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.pushButton.setIcon(icon5)
        self.pushButton.setIconSize(QtCore.QSize(40, 40))
        self.pushButton.setObjectName("pushButton")
        self.pushButton.clicked.connect(openFile)
        self.pushButton_2 = QtWidgets.QPushButton(self.dockWidgetContents_2)
        self.pushButton_2.setGeometry(QtCore.QRect(70, 10, 51, 51))
        self.pushButton_2.setText("")
        icon6 = QtGui.QIcon()
        icon6.addPixmap(QtGui.QPixmap(drawpng), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.pushButton_2.setIcon(icon6)
        self.pushButton_2.setIconSize(QtCore.QSize(36, 65))
        self.pushButton_2.setObjectName("pushButton_2")
        self.pushButton_2.clicked.connect(lambda: setBondType("SINGLE"))
        self.pushButton_3 = QtWidgets.QPushButton(self.dockWidgetContents_2)
        self.pushButton_3.setGeometry(QtCore.QRect(130, 10, 51, 51))
        self.pushButton_3.setText("")
        icon7 = QtGui.QIcon()
        icon7.addPixmap(QtGui.QPixmap(removeatom), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.pushButton_3.setIcon(icon7)
        self.pushButton_3.setIconSize(QtCore.QSize(36, 65))
        self.pushButton_3.setObjectName("pushButton_3")
        self.pushButton_3.setCheckable(True)
        self.pushButton_3.setStyleSheet("""
                                            QPushButton:checked {
                                        background-color: #6272a4; /* Change to your preferred checked color */
                                        color: #fff;
                                        border: 2px solid #22223b;
                                    }

                                        """)
        def checkforremove(checked):
            if checked:
                setAction("Remove")
                self.statusbar.showMessage("Remove Atom Mode")
            else:
                setAction("Add")
                self.statusbar.showMessage("Draw Mode")
        self.pushButton_3.toggled.connect(checkforremove)
        self.pushButton_3.toggled.connect(checkforremove)
        self.pushButton_4 = QtWidgets.QPushButton(self.dockWidgetContents_2)
        self.pushButton_4.setGeometry(QtCore.QRect(190, 10, 51, 51))
        self.pushButton_4.setText("")
        
        icon8 = QtGui.QIcon()
        icon8.addPixmap(QtGui.QPixmap(Stereopng), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        iconez = QtGui.QIcon()
        iconez.addPixmap(QtGui.QPixmap(Stereoezpng), QtGui.QIcon.Normal, QtGui.QIcon.Off)

        def checkforRS(checked):
            if checked:
                self.editor._chementitytype = "atom"
                setAction("RStoggle")
                self.statusbar.showMessage("R/S Toggle Mode")

            else:

                setAction("Add")
                self.statusbar.showMessage("Draw Mode")
            
        def checkforEZ(checked):
            if checked:
                self.editor._chementitytype = "atom"
                setAction("EZtoggle")
                self.statusbar.showMessage("E/Z Toggle Mode")

            else:

                setAction("Add")
                self.statusbar.showMessage("Draw Mode")
        self.pushButton_ez = QtWidgets.QPushButton(self.dockWidgetContents_2)
        self.pushButton_ez.setGeometry(QtCore.QRect(250, 10, 51, 51))
        self.pushButton_ez.setText("")
        self.pushButton_ez.setStatusTip("E/Z Toggle")
        self.pushButton_ez.toggled.connect(checkforEZ)
        self.pushButton_ez.setIcon(iconez)
        self.pushButton_ez.setIconSize(QtCore.QSize(36, 65))
        self.pushButton_ez.setObjectName("pushButton_4")
        self.pushButton_ez.setCheckable(True)
        self.pushButton_ez.setStyleSheet("""
                                          QPushButton:checked {
                                        background-color: #6272a4; /* Change to your preferred checked color */
                                        color: #fff;
                                        border: 2px solid #22223b;
                                    }

                                        """)
        

        self.pushButton_4.toggled.connect(checkforRS)
        self.pushButton_4.setIcon(icon8)
        self.pushButton_4.setIconSize(QtCore.QSize(36, 65))
        self.pushButton_4.setObjectName("pushButton_4")
        self.pushButton_4.setCheckable(True)
        self.pushButton_4.setStyleSheet("""
                                          QPushButton:checked {
                                        background-color: #6272a4; /* Change to your preferred checked color */
                                        color: #fff;
                                        border: 2px solid #22223b;
                                    }

                                        """)
        


        def cleanup_mol(self):
            mol = copy.deepcopy(self.mol)
            if self.sanitize_on_cleanup:
                Chem.SanitizeMol(mol)
            if self.kekulize_on_cleanup:
                Chem.Kekulize(mol)
            # if Chem.MolToCXSmiles(self.mol) != Chem.MolToCXSmiles(mol):
            self.mol = mol
        self.pushButton_5 = QtWidgets.QPushButton(self.dockWidgetContents_2)
        self.pushButton_5.setGeometry(QtCore.QRect(670, 10, 51, 51))
        self.pushButton_5.setText("")
        icon9 = QtGui.QIcon()
        icon9.addPixmap(QtGui.QPixmap(cleanstr), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.pushButton_5.setIcon(icon9)
        self.pushButton_5.setIconSize(QtCore.QSize(36, 65))
        self.pushButton_5.setObjectName("pushButton_5")
        self.pushButton_5.clicked.connect(self.editor.cleanup_mol)

        def canon_coords_and_draw(self):
            self.computeNewCoords(canonOrient=True, ignoreExisting=True)
            self._drawmol = copy.deepcopy(self._mol)  # Chem.Mol(self._mol.ToBinary())
            self.draw()


        self.pushButton_6 = QtWidgets.QPushButton(self.dockWidgetContents_2)
        self.pushButton_6.setGeometry(QtCore.QRect(310, 10, 51, 51))
        self.pushButton_6.setText("")
        icon10 = QtGui.QIcon()
        icon10.addPixmap(QtGui.QPixmap(xypng), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.pushButton_6.setIcon(icon10)
        self.pushButton_6.setIconSize(QtCore.QSize(36, 65))
        self.pushButton_6.setObjectName("pushButton_6")
        self.pushButton_6.clicked.connect(self.editor.canon_coords_and_draw)
        self.pushButton_7 = QtWidgets.QPushButton(self.dockWidgetContents_2)
        self.pushButton_7.setGeometry(QtCore.QRect(370, 10, 51, 51))
        self.pushButton_7.setText("")
        icon11 = QtGui.QIcon()
        icon11.addPixmap(QtGui.QPixmap(pluschargepng), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.pushButton_7.setIcon(icon11)
        self.pushButton_7.setIconSize(QtCore.QSize(43, 65))
        self.pushButton_7.setObjectName("pushButton_7")
        self.pushButton_7.setCheckable(True)
        self.pushButton_7.setStyleSheet("""
                                            QPushButton:checked {
                                        background-color: #6272a4; /* Change to your preferred checked color */
                                        color: #fff;
                                        border: 2px solid #22223b;
                                    }

                                        """)
        def checkforpos(checked):
            if checked:
                setAction("Increase Charge")
                self.statusbar.showMessage("Assign + Formal Charge")
            else:
                setAction("Add")
                self.statusbar.showMessage("Draw Mode")
        self.pushButton_7.toggled.connect(checkforpos)
        self.pushButton_8 = QtWidgets.QPushButton(self.dockWidgetContents_2)
        self.pushButton_8.setGeometry(QtCore.QRect(430, 10, 51, 51))
        self.pushButton_8.setText("")
        icon12 = QtGui.QIcon()
        icon12.addPixmap(QtGui.QPixmap(negativechargepng), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.pushButton_8.setIcon(icon12)
        self.pushButton_8.setIconSize(QtCore.QSize(43, 65))
        self.pushButton_8.setObjectName("pushButton_8")
        self.pushButton_8.setCheckable(True)
        self.pushButton_8.setStyleSheet("""
                                            QPushButton:checked {
                                        background-color: #6272a4; /* Change to your preferred checked color */
                                        color: #fff;
                                        border: 2px solid #22223b;
                                    }

                                        """)
        def checkforneg(checked):
            if checked:
                setAction("Decrease Charge")
                self.statusbar.showMessage("Assign - Formal Charge")
            else:
                setAction("Add")
                self.statusbar.showMessage("Draw Mode")
        self.pushButton_8.toggled.connect(checkforneg)
       

        self.pushButton_9 = QtWidgets.QPushButton(self.dockWidgetContents_2)
        self.pushButton_9.setGeometry(QtCore.QRect(490, 10, 51, 51))
        self.pushButton_9.setText("")
        icon13 = QtGui.QIcon()
        icon13.addPixmap(QtGui.QPixmap(hashmappng), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.pushButton_9.setIcon(icon13)
        self.pushButton_9.setIconSize(QtCore.QSize(36, 65))
        self.pushButton_9.setObjectName("pushButton_9")
        self.pushButton_9.setCheckable(True)
        self.pushButton_9.setStyleSheet("""
                                         QPushButton:checked {
                                        background-color: #6272a4; /* Change to your preferred checked color */
                                        color: #fff;
                                        border: 2px solid #22223b;
                                    }
                                        """)
        def setActionChecked(checked):
            if checked:
                self.pushButton_9.clicked.connect(lambda: setAction("Number Atom"))
                self.statusbar.showMessage("Atom Mapping Mode")
            else:
                self.pushButton_9.clicked.connect(lambda: setAction("Add"))
                self.statusbar.showMessage("Draw Mode")

        self.pushButton_9.toggled.connect(setActionChecked)
       
        def clearCanvas():
            self.editor.clearAtomSelection()
            self.editor.mol = None
            self.fileName = None
            self.statusBar().showMessage("Canvas Cleared")
            StyledMessageDialog("Chem-Parse", "Canvas Cleared" , self).exec()
        self.pushButton_10 = QtWidgets.QPushButton(self.dockWidgetContents_2)
        self.pushButton_10.setGeometry(QtCore.QRect(550, 10, 51, 51))
        self.pushButton_10.setText("")
        self.pushButton_10.clicked.connect(clearCanvas)
        icon14 = QtGui.QIcon()
        icon14.addPixmap(QtGui.QPixmap(trashpng), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.pushButton_10.setIcon(icon14)
        self.pushButton_10.setIconSize(QtCore.QSize(36, 65))
        self.pushButton_10.setObjectName("pushButton_10")
        self.pushButton_14 = QtWidgets.QPushButton(self.dockWidgetContents_2)
        self.pushButton_14.setGeometry(QtCore.QRect(610, 10, 51, 51))
        self.pushButton_14.setText("")
        icon15 = QtGui.QIcon()
        icon15.addPixmap(QtGui.QPixmap(cursormode), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.pushButton_14.setIcon(icon15)
        self.pushButton_14.setIconSize(QtCore.QSize(36, 65))
        self.pushButton_14.setObjectName("pushButton_14")
        self.pushButton_14.setCheckable(True)
        self.pushButton_14.setStyleSheet("""
                                        QPushButton:checked {
                                        background-color: #6272a4; /* Change to your preferred checked color */
                                        color: #fff;
                                        border: 2px solid #22223b;
                                    }

                                         """)
        self.pushButton_100 = QtWidgets.QPushButton(self.dockWidgetContents_2)
        self.pushButton_100.setGeometry(QtCore.QRect(730, 10, 51, 51))
        self.pushButton_100.setText("")
        icon15 = QtGui.QIcon()
        icon15.addPixmap(QtGui.QPixmap(copypng), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.pushButton_100.setIcon(icon15)
        self.pushButton_100.setIconSize(QtCore.QSize(36, 65))
        self.pushButton_100.setObjectName("pushButton_14")
        self.pushButton_100.setStyleSheet("""
                                        QPushButton:checked {
                                        background-color: #6272a4; /* Change to your preferred checked color */
                                        color: #fff;
                                        border: 2px solid #22223b;
                                    }

                                         """)
       
        self.pushButton_100.setToolTip("Copy SMILES to Table")
        
        def COPYtotable():
            selected_text = Chem.MolToSmiles(self.editor.mol, isomericSmiles=True)
            clipboard = QApplication.clipboard()
            clipboard.setText(selected_text)
            text = clipboard.text()
            mol = Chem.MolFromSmiles(text, sanitize=False)
            if mol:

                table = self.tableWidget  # Adjust if your table widget has a different name
                row_count = table.rowCount()
                inserted = False
                for row in range(row_count):
                    item = table.item(row, 1)  # 1 is the SMILES column
                    if item is None or not item.text().strip():
                        table.setItem(row, 0, QtWidgets.QTableWidgetItem(f"SMILES{row + 1}"))  # Name column
                        table.setItem(row, 1, QtWidgets.QTableWidgetItem(selected_text))        # SMILES column
                        inserted = True
                        break
                if not inserted:
                    # If no blank row, add a new row
                    table.insertRow(row_count)
                    table.setItem(row_count, 0, QtWidgets.QTableWidgetItem(f"SMILES{row_count + 1}"))  # Name column
                    table.setItem(row_count, 1, QtWidgets.QTableWidgetItem(selected_text))  
            else:
                self.editor.logger.warning(f"Failed to parse the content of the clipboard as a SMILES: {repr(text)}")

        self.pushButton_100.clicked.connect(COPYtotable)
        def setforselect(checked):
            if checked:
                self.pushButton_14.clicked.connect(lambda: setAction("Select"))
                self.statusbar.showMessage("Highlight Mode")
            else:
                self.pushButton_14.clicked.connect(lambda: setAction("Add"))
                self.statusbar.showMessage("Draw Mode")
        self.pushButton_14.toggled.connect(setforselect)
        self.dockWidget_2.setWidget(self.dockWidgetContents_2)
        self.addDockWidget(QtCore.Qt.TopDockWidgetArea, self.dockWidget_2)
        self.actionOpen = QtGui.QAction(self, triggered=openFile, shortcut=QKeySequence.Open)
        self.actionOpen.setObjectName("actionOpen")
        def undo(self):
          self.mol = self._prevmol
        def saveAsFile():
            self.fileName, filterName = QFileDialog.getSaveFileName(self, filter=self.filters_save, caption="Save file as")
            if self.fileName != "":
                if filterName == "MOL Files (*.mol *.mol)":
                    if not self.fileName.lower().endswith(".mol"):
                        self.fileName = self.fileName + ".mol"
                    Chem.MolToMolFile(self.editor.mol, str(self.fileName))
                    self.statusBar().showMessage("File saved as MolFile", 2000)
                    StyledMessageDialog("Chem-Parse", "File saved as MolFile", self ).exec()
                elif filterName == "SMILES Files (*.smi *.smi)":
                    if  self.fileName.lower().endswith(".smi"):
                        self.fileName = self.fileName + ".smi"
                    smiles = Chem.MolToSmiles(self.editor.mol)
                    with open(self.fileName, "w") as file:
                        file.write(smiles + "\n")
                    self.statusBar().showMessage("File saved as SMILES", 2000)
                    StyledMessageDialog("Chem-Parse", "File saved as SMILES", self ).exec()
                elif filterName == "PNG File (*.png *.png)":
                    if not self.fileName.lower().endswith(".png"):
                        self.fileName = self.fileName + ".png"
                    mol = self.editor.mol
                    if mol:
                        # Generate coordinates if needed
                        Chem.rdDepictor.Compute2DCoords(mol)
                         # Generate 2D coordinates if not present
                        rdDepictor.Compute2DCoords(mol)

                        # Set up drawer with transparent background
                        drawer = rdMolDraw2D.MolDraw2DCairo(500, 500)  # Adjust size as needed
                        drawer.drawOptions().clearBackground = False   # Important for transparency

                        # Draw the molecule
                        rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
                        drawer.FinishDrawing()

                        # Save PNG file
                        with open(self.fileName, "wb") as f:
                            f.write(drawer.GetDrawingText())

                        self.statusBar().showMessage("File saved as PNG", 2000)
                        StyledMessageDialog("Chem-Parse", "File saved as PNG", self).exec()
                    else:
                        self.statusBar().showMessage("No molecule to save", 2000)
                        StyledMessageDialog("Chem-Parse", "No molecule loaded to export", self).exec()
                elif filterName == "SDF File (*.sdf *.sdf)":
                    mol = self.editor.mol
                    if mol:
                        # Add hydrogens for 3D optimization
                        molH = Chem.AddHs(mol)

                        # Generate 3D coordinates
                        AllChem.EmbedMolecule(molH, AllChem.ETKDG())


                        # Write to SDF
                        writer = Chem.SDWriter(self.fileName)
                        writer.write(molH)
                        writer.close()

                        self.statusBar().showMessage("File saved as SDF", 2000)
                        StyledMessageDialog("Chem-Parse", "File saved as SDF", self).exec()
                    else:
                        self.statusBar().showMessage("No molecule to save", 2000)
                        StyledMessageDialog("Chem-Parse", "No molecule loaded to export", self).exec()
                else:
                    self.statusBar().showMessage("Invalid file format", 2000)
                    StyledMessageDialog("Chem-Parse", "Invalid file format", self ).exec()
        def saveFile():
            if self.fileName is not None:
                Chem.MolToMolFile(self.editor.mol, str(self.fileName))
            else:
                saveAsFile()

        
        def copy():
            selected_text = Chem.MolToSmiles(self.editor.mol, isomericSmiles=True)
            clipboard = QApplication.clipboard()
            clipboard.setText(selected_text)

        def paste():
            clipboard = QApplication.clipboard()
            text = clipboard.text()
            mol = Chem.MolFromSmiles(text, sanitize=False)
            if mol:
                try:
                    Chem.SanitizeMol(copy.deepcopy(mol))  # ).ToBinary()))
                except Exception as e:
                    self.editor.logger.warning(f"Pasted SMILES is not sanitizable: {e}")

                self.editor.assign_stereo_atoms(mol)
                Chem.rdmolops.SetBondStereoFromDirections(mol)

                self.editor.mol = mol
            else:
                self.editor.logger.warning(f"Failed to parse the content of the clipboard as a SMILES: {repr(text)}")

        def parse_smiles():

            selected_text = Chem.MolToSmiles(self.editor.mol, isomericSmiles=True)
            clipboard = QApplication.clipboard()
            clipboard.setText(selected_text)
            text = clipboard.text()
            mol = Chem.MolFromSmiles(text, sanitize=False)
            if mol:

                table = self.tableWidget  # Adjust if your table widget has a different name
                row_count = table.rowCount()
                inserted = False
                for row in range(row_count):
                    item = table.item(row, 1)  # 1 is the SMILES column
                    if item is None or not item.text().strip():
                        table.setItem(row, 0, QtWidgets.QTableWidgetItem(f"SMILES{row + 1}"))  # Name column
                        table.setItem(row, 1, QtWidgets.QTableWidgetItem(selected_text))        # SMILES column
                        inserted = True
                        break
                if not inserted:
                    # If no blank row, add a new row
                    table.insertRow(row_count)
                    table.setItem(row_count, 0, QtWidgets.QTableWidgetItem(f"SMILES{row_count + 1}"))  # Name column
                    table.setItem(row_count, 1, QtWidgets.QTableWidgetItem(selected_text))  
            else:
                self.editor.logger.warning(f"Failed to parse the content of the clipboard as a SMILES: {repr(text)}")


        self.actionSave = QtGui.QAction(self, triggered=saveFile, shortcut=QKeySequence.Save)
        self.actionSave.setObjectName("actionSave")
        self.actionSave_As = QtGui.QAction(self, triggered=saveAsFile, shortcut=QKeySequence.SaveAs)
        self.actionSave_As.setObjectName("actionSave_As")
        self.ActionUndo = QtGui.QAction(self, triggered=self.editor.undo, shortcut=QKeySequence.Undo)
        self.actionClean_Structure = QtGui.QAction(self, triggered=self.editor.cleanup_mol, shortcut=QKeySequence("Ctrl+Shift+C"))
        self.actionClean_Structure.setObjectName("actionClean_Structure")
        # ...existing code...
        def calc_prop():
          
          try:
            mol = self.editor.mol
            if mol:
                from rdkit.Chem import Descriptors, rdMolDescriptors

                mw = Descriptors.MolWt(mol)
                exact_mass = Descriptors.ExactMolWt(mol)
                formula = rdMolDescriptors.CalcMolFormula(mol)
                logp = Descriptors.MolLogP(mol)
                h_donors = Descriptors.NumHDonors(mol)
                h_acceptors = Descriptors.NumHAcceptors(mol)

                msg = (
                    f"<b>Chemical Formula:</b> {formula}<br>"
                    f"<b>Molecular Weight:</b> {mw:.2f} g/mol<br>"
                    f"<b>Exact Mass:</b> {exact_mass:.4f} Da<br>"
                    f"<b>LogP:</b> {logp:.2f}<br>"
                    f"<b>H-bond Donors:</b> {h_donors}<br>"
                    f"<b>H-bond Acceptors:</b> {h_acceptors}<br>"
    
                )

                dlg = StyledMessageDialog("Molecular Properties", msg, self)
                dlg.exec()
            else:
                StyledMessageDialog("Molecular Properties", "No molecule loaded.", self).exec()

          except Exception as e:
              StyledMessageDialog("Molecular Properties", "ERROR in calculation", self).exec()
       

        self.actionCalculate_Molecular_Properties = QtGui.QAction(self, triggered=calc_prop, shortcut=QKeySequence("Ctrl+M"))
        self.actionCalculate_Molecular_Properties.setObjectName("actionCalculate_Molecular_Properties")
      
        self.copysmiles = QtGui.QAction(self, triggered=copy, shortcut=QKeySequence.Copy)
        self.pastesmiles = QtGui.QAction(self, triggered=paste, shortcut=QKeySequence.Paste)
        self.parsesmiles = QtGui.QAction(self, triggered=parse_smiles, shortcut=QKeySequence("Ctrl+P"))
        self.clearcanva = QtGui.QAction(self, triggered=clearCanvas, shortcut=QKeySequence("Ctrl+L"))
        self.menuFile.addAction(self.actionOpen)
        self.menuFile.addAction(self.actionSave)
        self.menuFile.addAction(self.actionSave_As)
        self.menuFile.addAction(self.ActionUndo)
        self.menuHelp.addAction(self.actionClean_Structure)
        self.menuHelp.addAction(self.actionCalculate_Molecular_Properties)
       
        self.menuHelp.addAction(self.copysmiles)
        self.menuHelp.addAction(self.pastesmiles)
        self.menuHelp.addAction(self.parsesmiles)
        self.menuHelp.addAction(self.clearcanva)
       
        self.menubar.addAction(self.menuFile.menuAction())
        self.menubar.addAction(self.menuHelp.menuAction())
        #self.menubar.addAction(self.menuHelp_2.menuAction())
        
       
        self.retranslateUi()
        QtCore.QMetaObject.connectSlotsByName(self)
        self.showMaximized()
        self.show()

   
   

   

    def retranslateUi(self):
        _translate = QtCore.QCoreApplication.translate
        self.setWindowTitle(_translate("MainWindow", "Chem-Parse"))
        self.menuFile.setTitle(_translate("MainWindow", "File"))
        self.menuHelp.setTitle(_translate("MainWindow", "Structure"))
        self.menuHelp_2.setTitle(_translate("MainWindow", "Help"))
        self.pushButton_11.setStatusTip(_translate("MainWindow", "Single Bond"))
        self.pushButton_11.setWhatsThis(_translate("MainWindow", "Apply Single Bond (C-C)"))
        self.pushButton_12.setStatusTip(_translate("MainWindow", "Double Bond"))
        self.pushButton_12.setWhatsThis(_translate("MainWindow", "Click to Apply Double Bond"))
        self.pushButton_13.setStatusTip(_translate("MainWindow", "Triple Bond"))
        self.pushButton_13.setWhatsThis(_translate("MainWindow", "Click to Apply Triple Bond"))
        self.pushButton_R.setStatusTip(_translate("MainWindow", "R Group"))
        self.pushButton_R.setWhatsThis(_translate("MainWindow", "Add R Group"))
        self.pushButton_R.setText(_translate("MainWindow", "R"))
        self.pushButton_R_2.setStatusTip(_translate("MainWindow", "Carbon Atom"))
        self.pushButton_R_2.setWhatsThis(_translate("MainWindow", "Add C Atom"))
        self.pushButton_R_2.setText(_translate("MainWindow", "C"))
        self.pushButton_R_3.setStatusTip(_translate("MainWindow", "Hydrogen Atom"))
        self.pushButton_R_3.setWhatsThis(_translate("MainWindow", "Add H Atom"))
        self.pushButton_R_3.setText(_translate("MainWindow", "H"))
        self.pushButton_R_4.setStatusTip(_translate("MainWindow", "Oxygen Atom"))
        self.pushButton_R_4.setWhatsThis(_translate("MainWindow", "Add Oxygen Atom"))
        self.pushButton_R_4.setText(_translate("MainWindow", "O"))
        self.pushButton_R_5.setStatusTip(_translate("MainWindow", "Nitrogen Atom"))
        self.pushButton_R_5.setWhatsThis(_translate("MainWindow", "Add Nitrogen Atom"))
        self.pushButton_R_5.setText(_translate("MainWindow", "N"))
        self.pushButton_R_6.setStatusTip(_translate("MainWindow", "Sulphur Atom"))
        self.pushButton_R_6.setWhatsThis(_translate("MainWindow", "Add Sulphur Atom"))
        self.pushButton_R_6.setText(_translate("MainWindow", "S"))
        self.pushButton_R_7.setStatusTip(_translate("MainWindow", "Chlorine Atom"))
        self.pushButton_R_7.setWhatsThis(_translate("MainWindow", "Add Chlorine Atom"))
        self.pushButton_R_7.setText(_translate("MainWindow", "Cl"))
        self.pushButton_R_8.setStatusTip(_translate("MainWindow", "Fluorine Atom"))
        self.pushButton_R_8.setWhatsThis(_translate("MainWindow", "Add Fluorine Atom"))
        self.pushButton_R_8.setText(_translate("MainWindow", "F"))
        self.pushButton_R_9.setStatusTip(_translate("MainWindow", "Bromine Atom"))
        self.pushButton_R_9.setWhatsThis(_translate("MainWindow", "Add Bromine Atom"))
        self.pushButton_R_9.setText(_translate("MainWindow", "Br"))
        self.pushButton_R_10.setStatusTip(_translate("MainWindow", "Phosphorous Atom"))
        self.pushButton_R_10.setWhatsThis(_translate("MainWindow", "Add Phosphorous Atom"))
        self.pushButton_R_10.setText(_translate("MainWindow", "P"))
        self.pushButton_R_11.setStatusTip(_translate("MainWindow", "Add Iodine Atom"))
        self.pushButton_R_11.setWhatsThis(_translate("MainWindow", "Add Iodine Atom"))
        self.pushButton_R_11.setText(_translate("MainWindow", "I"))
        self.pushButton_R_12.setStatusTip(_translate("MainWindow", "Periodic Table"))
        self.pushButton_R_12.setWhatsThis(_translate("MainWindow", "Open Periodic Table"))
        self.pushButton.setStatusTip(_translate("MainWindow", "Open File ( .mol, .smi )"))
        self.pushButton_2.setStatusTip(_translate("MainWindow", "Draw Canvas Mode"))
        self.pushButton_2.setWhatsThis(_translate("MainWindow", "Draw Canvas Mode"))
        self.pushButton_3.setStatusTip(_translate("MainWindow", "Remove Atoms Mode"))
        self.pushButton_3.setWhatsThis(_translate("MainWindow", "Click Atom to Remove"))
        self.pushButton_4.setStatusTip(_translate("MainWindow", "Stereo R\\S"))
        self.pushButton_4.setWhatsThis(_translate("MainWindow", "Toggle Dashed / Wedged Bonds"))
        self.pushButton_5.setStatusTip(_translate("MainWindow", "Clean Structure"))
        self.pushButton_5.setWhatsThis(_translate("MainWindow", "Click to Clean Chemical Struture"))
        self.pushButton_6.setStatusTip(_translate("MainWindow", "Recalculate Coordinates"))
        self.pushButton_6.setWhatsThis(_translate("MainWindow", "Recalculate Coordinates"))
        self.pushButton_7.setStatusTip(_translate("MainWindow", " Positive Formal Charges"))
        self.pushButton_7.setWhatsThis(_translate("MainWindow", "Add Positive Formal Charges"))
        self.pushButton_8.setStatusTip(_translate("MainWindow", "Negative Formal Charges"))
        self.pushButton_8.setWhatsThis(_translate("MainWindow", "Add Negative Formal Charges"))
        self.pushButton_9.setStatusTip(_translate("MainWindow", "Atom Mapping"))
        self.pushButton_9.setWhatsThis(_translate("MainWindow", "Click Atom to Map Atom with Index"))
        self.pushButton_10.setStatusTip(_translate("MainWindow", "Clear Canvas"))
        self.pushButton_10.setWhatsThis(_translate("MainWindow", "Click to Clear The Canvas"))
        self.pushButton_14.setStatusTip(_translate("MainWindow", "Highlight Atoms"))
        self.pushButton_14.setWhatsThis(_translate("MainWindow", "Highlight Atoms"))
        self.actionOpen.setText(_translate("MainWindow", "Open"))
        self.actionSave.setText(_translate("MainWindow", "Save"))
        self.actionSave_As.setText(_translate("MainWindow", "Save As"))
        self.ActionUndo.setText(_translate("MainWindow", "Undo"))
        self.actionClean_Structure.setText(_translate("MainWindow", "Clean Structure"))
        self.actionCalculate_Molecular_Properties.setText(_translate("MainWindow", "Calculate Molecular Properties"))
        self.copysmiles.setText(_translate("MainWindow", "Copy SMILES"))
        self.pastesmiles.setText(_translate("MainWindow", "Paste SMILES"))
        self.parsesmiles.setText(_translate("MainWindow", "Parse SMILES to Table"))
        self.clearcanva.setText(_translate("MainWindow", "Clear Canvas"))


def show_exit_confirmation(self, event):
        msg_box = QtWidgets.QMessageBox()
        icon1 = QtGui.QIcon()
        icon1.addPixmap(QtGui.QPixmap(appicon), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        msg_box.setWindowTitle("Chem-Parse")
        msg_box.setWindowIcon(icon1)
        pixmap = QtGui.QPixmap(appicon).scaled(80, 80, QtCore.Qt.KeepAspectRatio, QtCore.Qt.SmoothTransformation)
        msg_box.setIconPixmap(pixmap)
        msg_box.setText("Are you sure you want to exit Chem-Parse?")
        msg_box.setStandardButtons(QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
        msg_box.setStyleSheet("""
            QMessageBox {
                background-color: #232629;
                color: #f8f8f2;
                font-size: 15px;
                border-radius: 10px;
            }
            QLabel {
                color: #f8f8f2;
                font-size: 15px;
            }
            QPushButton {
                background-color: #44475a;
                color: #f8f8f2;
                border-radius: 6px;
                padding: 6px 18px;
                font-size: 14px;
                min-width: 80px;
            }
            QPushButton:hover {
                background-color: #6272a4;
            }
        """)

        reply = msg_box.exec()
        if reply == QtWidgets.QMessageBox.Yes:
            QtWidgets.QApplication.quit()
        else:
            event.ignore()

   


def launch(loglevel="WARNING"):
    "Function that launches the mainWindow Application"
    try:
        myApp = QApplication(sys.argv)
        if len(sys.argv) > 1:
            mainWindow = MainWindow(fileName=sys.argv[1], loglevel=loglevel)
        else:
            mainWindow = MainWindow(loglevel=loglevel)

        # Properly override closeEvent
        def closeEvent(event):
            show_exit_confirmation(mainWindow, event)
        mainWindow.closeEvent = closeEvent

        if getattr(sys, 'frozen', False):
            pyi_splash.close()
        mainWindow.show()
        myApp.exec()
        sys.exit(0)
    except NameError:
        print("Name Error:", sys.exc_info()[1])
    except SystemExit:
        print("Closing Window...")
    except Exception:
        print(sys.exc_info()[1])


if __name__ == "__main__":
    launch(loglevel="DEBUG")