import sys
import json
from pathlib import Path

from PySide6.QtCore import QProcess
from PySide6.QtGui import QTextCursor
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QFormLayout, QLineEdit, QPushButton, QFileDialog, QTextEdit,
    QComboBox, QCheckBox, QLabel, QGroupBox, QMessageBox, QListWidget
)


def application_directory():
    if getattr(sys, "frozen", False):
        return Path(sys.executable).resolve().parent
    return Path(__file__).resolve().parent


ENERGY_UNITS = {
    "kJ/mol": 1,
    "kcal/mol": 2,
    "Ry": 3,
    "eV": 4,
    "Hartree": 5,
}


class EnGOATGUI(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("EnGOAT GUI")
        self.resize(600, 700)

        self.gui_directory = application_directory()
        self.cube_path = None

        self.process = QProcess(self)
        self.process.readyReadStandardOutput.connect(self.read_process_output)
        self.process.readyReadStandardError.connect(self.read_process_error)
        self.process.finished.connect(self.process_finished)
        self.process.errorOccurred.connect(self.process_error)

        self.create_ui()

    def create_ui(self):
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QVBoxLayout(central_widget)

        # Cube file
        cube_group = QGroupBox("Cube File")
        cube_layout = QVBoxLayout()
        cube_row = QHBoxLayout()

        self.cube_edit = QLineEdit()
        self.cube_edit.setPlaceholderText("Select a .cube file...")
        self.cube_edit.setReadOnly(True)

        browse_button = QPushButton("Browse")
        browse_button.clicked.connect(self.select_cube_file)

        cube_row.addWidget(self.cube_edit)
        cube_row.addWidget(browse_button)
        cube_layout.addLayout(cube_row)

        self.directory_label = QLabel("Working directory: Not selected")
        self.directory_label.setWordWrap(True)
        cube_layout.addWidget(self.directory_label)

        cube_group.setLayout(cube_layout)
        main_layout.addWidget(cube_group)

        # TuTraSt parameters
        tutrast_group = QGroupBox("TuTraSt Parameters")
        tutrast_form = QFormLayout()

        self.e_unit = QComboBox()
        self.e_unit.addItems(ENERGY_UNITS.keys())
        self.e_unit.setCurrentText("kJ/mol")
        tutrast_form.addRow("Energy unit:", self.e_unit)

        self.e_step = QLineEdit("2")
        tutrast_form.addRow("Energy step [kJ/mol]:", self.e_step)

        self.e_cutoff = QLineEdit("50.0")
        tutrast_form.addRow("Energy cutoff [kJ/mol]:", self.e_cutoff)

        self.dynamic_e_step = QComboBox()
        self.dynamic_e_step.addItems(["off", "on"])
        self.dynamic_e_step.setCurrentText("off")
        tutrast_form.addRow("Dynamic E step:", self.dynamic_e_step)

        tutrast_group.setLayout(tutrast_form)
        main_layout.addWidget(tutrast_group)

        # KMC parameters
        kmc_group = QGroupBox("kMC Parameters")
        kmc_form = QFormLayout()

        self.kmc_run = QCheckBox("Enable kMC")
        self.kmc_run.setChecked(True)
        kmc_form.addRow("kMC:", self.kmc_run)

        temperature_widget = QWidget()
        temperature_layout = QVBoxLayout(temperature_widget)
        temperature_layout.setContentsMargins(0, 0, 0, 0)

        self.temperature_list = QListWidget()
        self.temperature_list.addItem("300")
        temperature_layout.addWidget(self.temperature_list)

        temperature_controls = QHBoxLayout()

        self.temperature_input = QLineEdit()
        self.temperature_input.setPlaceholderText("Temperature [K]")
        self.temperature_input.returnPressed.connect(self.add_temperature)

        add_temperature_button = QPushButton("Add")
        add_temperature_button.clicked.connect(self.add_temperature)

        remove_temperature_button = QPushButton("Remove")
        remove_temperature_button.clicked.connect(self.remove_temperature)

        temperature_controls.addWidget(self.temperature_input)
        temperature_controls.addWidget(add_temperature_button)
        temperature_controls.addWidget(remove_temperature_button)
        temperature_layout.addLayout(temperature_controls)

        kmc_form.addRow("Temperatures (K):", temperature_widget)

        self.steps = QLineEdit("1000000")
        kmc_form.addRow("No. Steps:", self.steps)

        self.runs = QLineEdit("5")
        kmc_form.addRow("No. Runs:", self.runs)

        self.particle_mass = QLineEdit("6.941e-3")
        kmc_form.addRow("Particle mass [kg/mol]:", self.particle_mass)

        kmc_group.setLayout(kmc_form)
        main_layout.addWidget(kmc_group)

        # Status
        status_layout = QHBoxLayout()

        status_layout.addWidget(QLabel("Status:"))

        self.status_label = QLabel("Ready")
        status_layout.addWidget(self.status_label)

        status_layout.addStretch()
        main_layout.addLayout(status_layout)

        # Run / Stop buttons
        button_layout = QHBoxLayout()

        self.run_button = QPushButton("Run EnGOAT")
        self.run_button.setMinimumHeight(40)
        self.run_button.clicked.connect(self.run_engoat)

        self.stop_button = QPushButton("Stop")
        self.stop_button.setMinimumHeight(40)
        self.stop_button.setEnabled(False)
        self.stop_button.clicked.connect(self.stop_engoat)

        button_layout.addWidget(self.run_button)
        button_layout.addWidget(self.stop_button)
        main_layout.addLayout(button_layout)

        # Output
        output_group = QGroupBox("Output")
        output_layout = QVBoxLayout()

        self.output = QTextEdit()
        self.output.setReadOnly(True)
        self.output.setMinimumHeight(200)

        output_layout.addWidget(self.output)
        output_group.setLayout(output_layout)
        main_layout.addWidget(output_group)

    def select_cube_file(self):
        if self.process.state() != QProcess.NotRunning:
            QMessageBox.warning(
                self,
                "EnGOAT is running",
                "Please stop EnGOAT before selecting another cube file."
            )
            return

        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Select Cube File",
            "",
            "Cube Files (*.cube);;All Files (*)"
        )

        if not file_path:
            return

        self.cube_path = Path(file_path)
        self.cube_edit.setText(str(self.cube_path))

        working_directory = self.cube_path.parent
        self.directory_label.setText(f"Working directory: {working_directory}")

        self.output.append("Selected cube file:")
        self.output.append(str(self.cube_path))
        self.output.append("")
        self.output.append("Working directory:")
        self.output.append(str(working_directory))
        self.output.append("")

    def add_temperature(self):
        text = self.temperature_input.text().strip()

        if not text:
            return

        try:
            temperature = int(text)
            if temperature <= 0:
                raise ValueError
        except ValueError:
            QMessageBox.warning(
                self,
                "Invalid temperature",
                "Please enter a positive integer temperature."
            )
            return

        existing_temperatures = [
            int(self.temperature_list.item(i).text())
            for i in range(self.temperature_list.count())
        ]

        if temperature in existing_temperatures:
            QMessageBox.information(
                self,
                "Temperature already exists",
                f"{temperature} K is already in the list."
            )
            return

        self.temperature_list.addItem(str(temperature))
        self.temperature_input.clear()
        self.temperature_input.setFocus()

    def remove_temperature(self):
        selected_items = self.temperature_list.selectedItems()

        if not selected_items:
            return

        for item in selected_items:
            self.temperature_list.takeItem(self.temperature_list.row(item))

    def get_temperatures(self):
        return [
            int(self.temperature_list.item(i).text())
            for i in range(self.temperature_list.count())
        ]

    def create_input_data(self):
        temperatures = sorted(self.get_temperatures())

        if not temperatures:
            raise ValueError("At least one temperature must be specified.")

        return {
            "TuTraSt": {
                "E_unit": self.e_unit.currentText(),
                "E_step": float(self.e_step.text()),
                "E_cutoff": float(self.e_cutoff.text()),
                "Dynamic_E_step": self.dynamic_e_step.currentText(),
            },
            "kmc": {
                "run": self.kmc_run.isChecked(),
                "temperatures": temperatures,
                "steps": int(self.steps.text()),
                "runs": int(self.runs.text()),
                "particle_mass": float(self.particle_mass.text()),
            },
            "cube_file": self.cube_path.name,
        }

    def run_engoat(self):
        if self.cube_path is None:
            QMessageBox.warning(
                self,
                "No cube file",
                "Please select a cube file first."
            )
            return

        if self.process.state() != QProcess.NotRunning:
            QMessageBox.warning(
                self,
                "EnGOAT is already running",
                "EnGOAT is already running."
            )
            return

        working_directory = self.cube_path.parent
        engoat_script = self.gui_directory / "EnGOAT.py"

        if not engoat_script.exists():
            QMessageBox.critical(
                self,
                "EnGOAT.py not found",
                f"Could not find EnGOAT.py at:\n\n{engoat_script}"
            )
            return

        input_file = working_directory / "input.json"

        try:
            input_data = self.create_input_data()

            with open(input_file, "w") as f:
                json.dump(input_data, f, indent=4)

        except ValueError as e:
            QMessageBox.critical(
                self,
                "Invalid input",
                f"Please check the input values.\n\n{e}"
            )
            return

        except OSError as e:
            QMessageBox.critical(
                self,
                "Could not create input file",
                str(e)
            )
            return

        self.output.clear()
        self.output.append("Starting EnGOAT...\n")
        self.output.append(f"Cube file:\n{self.cube_path}\n")
        self.output.append(f"Input file:\n{input_file}\n")
        self.output.append(f"EnGOAT.py:\n{engoat_script}\n")
        self.output.append(f"Working directory:\n{working_directory}\n")
        self.output.append("----------------------------------------\n")

        self.process.setWorkingDirectory(str(working_directory))

        if getattr(sys, "frozen", False):
            worker_executable = self.gui_directory / "EnGOAT_worker.exe"

            if not worker_executable.exists():
                QMessageBox.critical(
                    self,
                    "EnGOAT worker not found",
                    f"Could not find EnGOAT_worker.exe at:\n\n{worker_executable}"
                )
                return

            self.process.setProgram(str(worker_executable))
            self.process.setArguments([])

        else:
            self.process.setProgram(sys.executable)
            self.process.setArguments(["-u", str(engoat_script)])

        self.run_button.setEnabled(False)
        self.stop_button.setEnabled(True)
        self.status_label.setText("Running...")

        self.process.start()

    def read_process_output(self):
        data = self.process.readAllStandardOutput()
        text = bytes(data).decode("utf-8", errors="replace")
        self.append_output(text)

    def read_process_error(self):
        data = self.process.readAllStandardError()
        text = bytes(data).decode("utf-8", errors="replace")
        self.append_output(text)

    def append_output(self, text):
        self.output.moveCursor(QTextCursor.End)
        self.output.insertPlainText(text)
        self.output.moveCursor(QTextCursor.End)
        self.output.ensureCursorVisible()

    def process_finished(self, exit_code, exit_status):
        self.output.append("\n----------------------------------------")

        if exit_code == 0:
            self.output.append("EnGOAT finished successfully.")
            self.status_label.setText("Finished")
        else:
            self.output.append(f"EnGOAT exited with code {exit_code}.")
            self.status_label.setText(f"Failed (exit code {exit_code})")

        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)

    def process_error(self, error):
        self.status_label.setText("Error")

        self.output.append("\n----------------------------------------")
        self.output.append("Error while running EnGOAT:")
        self.output.append(self.process.errorString())

        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)

    def stop_engoat(self):
        if self.process.state() == QProcess.NotRunning:
            return

        reply = QMessageBox.question(
            self,
            "Stop EnGOAT?",
            "Are you sure you want to stop the EnGOAT calculation?"
        )

        if reply != QMessageBox.Yes:
            return

        self.output.append("\nStopping EnGOAT...")

        self.process.terminate()

        if not self.process.waitForFinished(2000):
            self.process.kill()

        self.status_label.setText("Stopped")
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)

    def closeEvent(self, event):
        if self.process.state() != QProcess.NotRunning:
            reply = QMessageBox.question(
                self,
                "EnGOAT is running",
                "EnGOAT is still running.\n\n"
                "Do you want to stop it and close the GUI?"
            )

            if reply != QMessageBox.Yes:
                event.ignore()
                return

            self.process.terminate()

            if not self.process.waitForFinished(2000):
                self.process.kill()

        event.accept()


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = EnGOATGUI()
    window.show()
    sys.exit(app.exec())