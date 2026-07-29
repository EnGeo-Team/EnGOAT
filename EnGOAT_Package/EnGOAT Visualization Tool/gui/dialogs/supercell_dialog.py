import numpy as np

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QDialog,
    QVBoxLayout,
    QHBoxLayout,
    QGridLayout,
    QGroupBox,
    QLabel,
    QPushButton,
    QSpinBox,
    QCheckBox
)


class SupercellDialog(QDialog):
    def __init__(self, project, parent=None):
        super().__init__(parent)

        self.project = project

        self.setWindowTitle("Create Supercell")
        self.resize(450, 300)
        self.setMinimumSize(400, 250)

        self.setStyleSheet("""
            QGroupBox {
                font-weight: bold;
                border: 1px solid #c8c8c8;
                border-radius: 6px;
                margin-top: 8px;
                padding-top: 8px;
            }

            QGroupBox::title {
                subcontrol-origin: margin;
                left: 10px;
                padding: 0 4px;
            }
        """)

        self.create_ui()

    def create_ui(self):

        layout = QVBoxLayout(self)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(8)

        layout.addWidget(self.create_unit_cell_group())
        layout.addWidget(self.create_supercell_group())

        layout.addStretch()

        button_layout = QHBoxLayout()
        button_layout.addStretch()

        create_button = QPushButton("Create Supercell")
        create_button.clicked.connect(self.create_supercell)

        button_layout.addWidget(create_button)

        layout.addLayout(button_layout)

        
    def create_unit_cell_group(self):

        group = QGroupBox("Unit Cell")

        layout = QGridLayout(group)
        layout.setHorizontalSpacing(24)
        layout.setVerticalSpacing(8)

        metadata = self.project.metadata

        a_vec = np.asarray(metadata["grid_vectors"][0]) * metadata["grid_shape"][0]
        b_vec = np.asarray(metadata["grid_vectors"][1]) * metadata["grid_shape"][1]
        c_vec = np.asarray(metadata["grid_vectors"][2]) * metadata["grid_shape"][2]

        a = np.linalg.norm(a_vec)
        b = np.linalg.norm(b_vec)
        c = np.linalg.norm(c_vec)

        alpha = self.angle(b_vec, c_vec)
        beta = self.angle(a_vec, c_vec)
        gamma = self.angle(a_vec, b_vec)

        entries = [
            ("a", f"{a:.4f} Å", "α", f"{alpha:.2f}°"),
            ("b", f"{b:.4f} Å", "β", f"{beta:.2f}°"),
            ("c", f"{c:.4f} Å", "γ", f"{gamma:.2f}°"),
        ]

        for row, (l1, v1, l2, v2) in enumerate(entries):

            label1 = QLabel(f"<b>{l1}</b>")
            value1 = QLabel(v1)
            label2 = QLabel(f"<b>{l2}</b>")
            value2 = QLabel(v2)

            for widget in (label1, value1, label2, value2):
                widget.setMinimumHeight(22)

            layout.addWidget(label1, row, 0)
            layout.addWidget(value1, row, 1)
            layout.addWidget(label2, row, 2)
            layout.addWidget(value2, row, 3)

        layout.setColumnStretch(1, 1)
        layout.setColumnStretch(3, 1)

        return group


    def create_supercell_group(self):

        group = QGroupBox("Supercell")

        main_layout = QVBoxLayout(group)
        main_layout.setContentsMargins(10, 10, 10, 10)
        main_layout.setSpacing(8)

        # a b c row
        cell_layout = QHBoxLayout()
        cell_layout.setSpacing(18)

        Na, Nb, Nc = self.project.supercell

        self.a_spin = self.create_spinbox(Na)
        self.b_spin = self.create_spinbox(Nb)
        self.c_spin = self.create_spinbox(Nc)

        cell_layout.addWidget(QLabel("<b>a:</b>"))
        cell_layout.addWidget(self.a_spin)

        cell_layout.addWidget(QLabel("<b>b:</b>"))
        cell_layout.addWidget(self.b_spin)

        cell_layout.addWidget(QLabel("<b>c:</b>"))
        cell_layout.addWidget(self.c_spin)

        cell_layout.addStretch()

        main_layout.addLayout(cell_layout)

        # UC grid row
        grid_layout = QHBoxLayout()

        self.uc_grid_checkbox = QCheckBox("Display UC grid")

        current_state = self.project.visibility.get("UC_grid", False)
        self.uc_grid_checkbox.setChecked(current_state)

        self.uc_grid_checkbox.stateChanged.connect(
            lambda state: self.project.set_visibility(
                "UC_grid",
                state
            )
        )

        grid_layout.addWidget(self.uc_grid_checkbox)
        grid_layout.addStretch()

        main_layout.addLayout(grid_layout)

        return group
    
    def create_spinbox(self, N):

        box = QSpinBox()

        box.setMinimum(1)
        box.setSingleStep(1)
        box.setValue(N)

        return box

    def create_supercell(self):

        supercell = np.array([
            int(self.a_spin.value()),
            int(self.b_spin.value()),
            int(self.c_spin.value()),
        ])

        self.project.create_supercell(supercell)

        self.accept()

    @staticmethod
    def angle(v1, v2):

        v1 = np.asarray(v1, dtype=float)
        v2 = np.asarray(v2, dtype=float)

        cosang = np.dot(v1, v2) / (
            np.linalg.norm(v1) * np.linalg.norm(v2)
        )

        cosang = np.clip(cosang, -1.0, 1.0)

        return np.degrees(np.arccos(cosang))

