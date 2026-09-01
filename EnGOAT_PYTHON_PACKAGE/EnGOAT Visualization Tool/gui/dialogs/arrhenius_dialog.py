from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QGridLayout,
    QGroupBox, QLabel, QPushButton, QCheckBox,
    QSizePolicy, QFrame
)
from PySide6.QtCore import Qt
import numpy as np
from plotting.matplotlib.plot_arrhenius import plot_arrhenius


class ArrheniusDialog(QDialog):
    def __init__(self, project, parent=None):
        super().__init__(parent)
        self.project = project
        self.direction_checkboxes = {}
        self.activation_labels = {}
        self.setWindowTitle("Arrhenius Analysis")
        self.resize(620, 560)
        self.setMinimumSize(560, 480)
        self.setStyleSheet("""
            QDialog { background-color: #f7f7f7; }
            QGroupBox { font-weight: bold; border: 1px solid #c8c8c8; border-radius: 7px; margin-top: 10px; padding-top: 10px; background-color: white; }
            QGroupBox::title { subcontrol-origin: margin; left: 10px; padding: 0 5px; background-color: #f7f7f7; }
            QLabel { color: #303030; }
            QLabel#description { color: #666666; font-weight: normal; }
            QLabel#direction { font-size: 13px; font-weight: bold; }
            QLabel#axis_type { color: #777777; font-size: 11px; font-weight: normal; }
            QLabel#barrier { color: #303030; font-size: 12px; font-weight: normal; }
            QLabel#barrier_value { color: #222222; font-size: 12px; font-weight: bold; }
            QCheckBox { spacing: 6px; }
            QPushButton { min-height: 30px; padding-left: 14px; padding-right: 14px; border: 1px solid #b8b8b8; border-radius: 5px; background-color: white; }
            QPushButton:hover { background-color: #eeeeee; }
            QPushButton:pressed { background-color: #e2e2e2; }
            QPushButton#plot_button { font-weight: bold; min-height: 34px; }
            QFrame#direction_card { border: 1px solid #dddddd; border-radius: 5px; background-color: #fafafa; }
        """)
        self.create_ui()

    def create_ui(self):
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(14, 14, 14, 14)
        main_layout.setSpacing(10)

        title = QLabel("Arrhenius Analysis")
        title.setStyleSheet("QLabel { font-size: 18px; font-weight: bold; color: #222222; }")
        main_layout.addWidget(title)

        description = QLabel("Select the diffusion directions to include in the Arrhenius plot.")
        description.setObjectName("description")
        main_layout.addWidget(description)

        direction_group = QGroupBox("Diffusion directions")
        direction_layout = QVBoxLayout(direction_group)
        direction_layout.setContentsMargins(10, 10, 10, 10)
        direction_layout.setSpacing(6)

        header = QGridLayout()
        header.setColumnMinimumWidth(0, 150)
        header.setColumnMinimumWidth(1, 170)
        header.setColumnMinimumWidth(2, 150)
        header.addWidget(self.create_header_label("Direction"), 0, 0)
        header.addWidget(self.create_header_label("Activation barrier"), 0, 1)
        header.addWidget(self.create_header_label("Select"), 0, 2, alignment=Qt.AlignCenter)
        direction_layout.addLayout(header)

        direction_layout.addWidget(self.create_section_label("Crystal axes"))
        for direction in ["a", "b", "c"]:
            direction_layout.addWidget(self.create_direction_row(direction, "Crystal axis"))

        direction_layout.addSpacing(5)
        direction_layout.addWidget(self.create_section_label("Cartesian axes"))
        for direction in ["x", "y", "z"]:
            direction_layout.addWidget(self.create_direction_row(direction, "Cartesian axis"))

        direction_layout.addSpacing(5)
        direction_layout.addWidget(self.create_section_label("Isotropic"))
        direction_layout.addWidget(self.create_direction_row("3D", "Isotropic diffusion"))
        direction_layout.addStretch()
        main_layout.addWidget(direction_group, stretch=1)

        activation_energies = self.calculate_activation_energies(self.project.kMC_data)
        self.update_activation_barriers(activation_energies)

        button_layout = QHBoxLayout()
        button_layout.setSpacing(8)
        button_layout.addStretch()

        plot_button = QPushButton("Arrhenius plot")
        plot_button.setObjectName("plot_button")
        plot_button.setSizePolicy(QSizePolicy.Minimum, QSizePolicy.Fixed)
        plot_button.clicked.connect(self.plot_arrhenius)
        button_layout.addWidget(plot_button)

        main_layout.addLayout(button_layout)

    def create_direction_row(self, direction, axis_type):
        frame = QFrame()
        frame.setObjectName("direction_card")

        layout = QGridLayout(frame)
        layout.setContentsMargins(8, 5, 8, 5)
        layout.setColumnMinimumWidth(0, 150)
        layout.setColumnMinimumWidth(1, 170)
        layout.setColumnMinimumWidth(2, 150)

        direction_layout = QHBoxLayout()
        direction_layout.setSpacing(8)

        checkbox = QCheckBox()
        checkbox.setChecked(False)

        direction_label = QLabel(direction)
        direction_label.setObjectName("direction")

        axis_label = QLabel(axis_type)
        axis_label.setObjectName("axis_type")

        direction_layout.addWidget(checkbox)
        direction_layout.addWidget(direction_label)
        direction_layout.addWidget(axis_label)
        direction_layout.addStretch()
        layout.addLayout(direction_layout, 0, 0)

        barrier_layout = QHBoxLayout()
        barrier_layout.setSpacing(4)

        barrier_label = QLabel("—")
        barrier_label.setObjectName("barrier_value")

        barrier_unit = QLabel("kJ/mol")
        barrier_unit.setObjectName("barrier")

        barrier_layout.addWidget(barrier_label)
        barrier_layout.addWidget(barrier_unit)
        barrier_layout.addStretch()
        layout.addLayout(barrier_layout, 0, 1)

        selection_label = QLabel("Not selected")
        selection_label.setStyleSheet("QLabel { color: #888888; font-size: 11px; }")

        selection_layout = QHBoxLayout()
        selection_layout.addStretch()
        selection_layout.addWidget(selection_label)
        layout.addLayout(selection_layout, 0, 2, alignment=Qt.AlignRight)

        checkbox.toggled.connect(lambda checked, label=selection_label: label.setText("Selected" if checked else "Not selected"))

        self.direction_checkboxes[direction] = checkbox
        self.activation_labels[direction] = barrier_label

        return frame

    @staticmethod
    def create_header_label(text):
        label = QLabel(f"<b>{text}</b>")
        label.setStyleSheet("QLabel { color: #666666; font-size: 11px; }")
        return label

    @staticmethod
    def create_section_label(text):
        label = QLabel(text)
        label.setStyleSheet("QLabel { color: #555555; font-size: 11px; font-weight: bold; padding-top: 3px; padding-bottom: 1px; }")
        return label

    def get_selected_directions(self):
        return [direction for direction, checkbox in self.direction_checkboxes.items() if checkbox.isChecked()]

    def update_activation_barriers(self, activation_energies):
        for direction, label in self.activation_labels.items():
            if direction in activation_energies:
                label.setText(f"{activation_energies[direction]:.2f}")
                label.setStyleSheet("QLabel { color: #222222; font-size: 12px; font-weight: bold; }")
            else:
                label.setText("—")

    def plot_arrhenius(self):
        selected = self.get_selected_directions()
        if selected:
            plot_arrhenius(self.project.kMC_data, selected)

    @staticmethod
    def calculate_activation_energies(kMC_data, directions=None):
        R = 8.314462618
        if directions is None:
            directions = ["a", "b", "c", "x", "y", "z", "3D"]

        temperatures = []
        for temperature in kMC_data:
            try:
                temperatures.append(float(temperature))
            except (TypeError, ValueError):
                continue

        temperatures = np.array(sorted(temperatures))
        activation_energies = {}

        for direction in directions:
            D_values = []
            valid_temperatures = []

            for temperature in temperatures:
                temperature_key = str(int(temperature))
                if temperature_key in kMC_data:
                    data = kMC_data[temperature_key]
                elif temperature in kMC_data:
                    data = kMC_data[temperature]
                else:
                    continue

                try:
                    D = float(data["D_tot"][direction]["D"])
                except (KeyError, TypeError, ValueError):
                    continue

                if D <= 0:
                    continue

                valid_temperatures.append(temperature)
                D_values.append(D)

            if len(D_values) < 2:
                continue

            x = 1000.0 / np.asarray(valid_temperatures)
            y = np.log10(np.asarray(D_values))
            slope, _ = np.polyfit(x, y, 1)
            activation_energies[direction] = -slope * 2.303 * R

        return activation_energies