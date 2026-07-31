from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QGridLayout, QGroupBox, QLabel,
    QPushButton, QCheckBox, QScrollArea, QWidget, QFrame, QSizePolicy
)
from PySide6.QtCore import Qt
import numpy as np

from plotting.matplotlib.plot_msd import plot_tunnel_msd
import matplotlib.pyplot as plt


class DiffusionDialog(QDialog):
    def __init__(self, project, temperature, parent=None):
        super().__init__(parent)
        self.project = project
        self.temperature = temperature
        self.data = project.kMC_data[temperature]
        self.selected_directions = {}

        self.setWindowTitle(f"Diffusion — T = {float(temperature):g}")
        self.resize(720, 800)
        self.setMinimumSize(600, 500)

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
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(10, 10, 10, 10)
        main_layout.setSpacing(8)

        main_layout.addWidget(self.create_diffusion_group("Total Diffusion Coefficient", self.data["D_tot"]))

        tunnel_group = QGroupBox("Tunnel Systems")
        tunnel_layout = QVBoxLayout(tunnel_group)
        tunnel_layout.setContentsMargins(6, 8, 6, 6)

        tunnel_scroll = QScrollArea()
        tunnel_scroll.setWidgetResizable(True)
        tunnel_scroll.setFrameShape(QFrame.NoFrame)

        tunnel_container = QWidget()
        tunnel_container_layout = QVBoxLayout(tunnel_container)
        tunnel_container_layout.setContentsMargins(2, 2, 2, 2)
        tunnel_container_layout.setSpacing(6)

        for tunnel_ID, tunnel_data in self.data["tunnel_systems"].items():
            tunnel_container_layout.addWidget(self.create_tunnel_system_widget(tunnel_ID, tunnel_data))

        tunnel_container_layout.addStretch()
        tunnel_scroll.setWidget(tunnel_container)
        tunnel_layout.addWidget(tunnel_scroll)
        main_layout.addWidget(tunnel_group, stretch=1)

    def create_diffusion_group(self, title, data, selectable=False):
        group = QGroupBox(title)
        layout = QVBoxLayout(group)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(4)
        layout.addLayout(self.create_direction_table(data, selectable))
        return group

    def create_direction_table(self, data, selectable=False, tunnel_ID=None):
        table = QGridLayout()
        table.setHorizontalSpacing(14)
        table.setVerticalSpacing(1)
        table.setColumnMinimumWidth(0, 80)
        table.setColumnMinimumWidth(1, 125)
        table.setColumnMinimumWidth(2, 125)

        if selectable:
            table.setColumnMinimumWidth(3, 34)

        table.addWidget(self.create_header_label("Direction"), 0, 0)
        table.addWidget(self.create_header_label("D [cm²/s]"), 0, 1)
        table.addWidget(self.create_header_label("SD [cm²/s]"), 0, 2)

        if selectable:
            plot_selected_button = QPushButton("Plot selected MSD")
            plot_selected_button.setSizePolicy(QSizePolicy.Maximum, QSizePolicy.Fixed)
            plot_selected_button.clicked.connect(lambda checked=False, tunnel_ID=tunnel_ID: self.plot_selected_msd(tunnel_ID))
            table.addWidget(plot_selected_button, 0, 3, alignment=Qt.AlignCenter)

        direction_groups = [
            ("Crystal axes", ["a", "b", "c"]),
            ("Cartesian axes", ["x", "y", "z"]),
            ("Isotropic", ["3D"]),
        ]

        row = 1

        for group_name, group_directions in direction_groups:
            available_directions = [direction for direction in group_directions if direction in data]

            if not available_directions:
                continue

            group_label = QLabel(f"<b>{group_name}</b>")
            group_label.setProperty("class", "section-label")
            table.addWidget(group_label, row, 0, 1, 4 if selectable else 3)
            row += 1

            for direction in available_directions:
                direction_label = "D" if direction == "3D" else direction

                table.addWidget(QLabel(direction_label), row, 0)
                table.addWidget(self.create_value_label(data[direction]["D"]), row, 1)
                table.addWidget(self.create_value_label(data[direction]["sd"], prefix="± "), row, 2)

                if selectable:
                    checkbox = QCheckBox()
                    table.addWidget(checkbox, row, 3, alignment=Qt.AlignCenter)
                    self.selected_directions[tunnel_ID][direction] = checkbox

                row += 1

        return table

    def create_tunnel_system_widget(self, tunnel_ID, tunnel_data):
        group = QGroupBox(f"Tunnel system {tunnel_ID}")

        layout = QVBoxLayout(group)
        layout.setContentsMargins(8, 7, 8, 7)
        layout.setSpacing(4)

        header_layout = QHBoxLayout()
        header_layout.setSpacing(12)

        weight_label = QLabel(f"<b>Weight:</b> {self.format_weight(tunnel_data['weight'])}")
        header_layout.addWidget(weight_label)
        header_layout.addStretch()
        layout.addLayout(header_layout)

        self.selected_directions[tunnel_ID] = {}

        layout.addLayout(self.create_direction_table(tunnel_data["directions"], selectable=True, tunnel_ID=tunnel_ID))

        return group

    def get_selected_directions(self, tunnel_ID):
        return [
            direction
            for direction, checkbox in self.selected_directions[tunnel_ID].items()
            if checkbox.isChecked()
        ]

    def plot_selected_msd(self, tunnel_ID):
        selected = self.get_selected_directions(tunnel_ID)

        if not selected:
            return

        directions = ["a", "b", "c", "x", "y", "z", "3D"]
        conditions = self.get_msd_fit_conditions(self.project.metadata["grid_vectors"], self.project.metadata["grid_shape"])
        tab10 = plt.get_cmap("tab10")
        direction_colors = {
            direction: tab10(i)
            for i, direction in enumerate(directions)
        }

        plot_tunnel_msd(self.data, tunnel_ID, selected, direction_colors, conditions)
        return None


    def create_header_label(self, text):
        label = QLabel(f"<b>{text}</b>")
        label.setAlignment(Qt.AlignLeft)
        return label

    def create_value_label(self, value, prefix=""):
        label = QLabel(f"{prefix}{self.format_diffusion(value)}")
        label.setProperty("class", "value-label")
        label.setAlignment(Qt.AlignLeft)
        return label

    @staticmethod
    def format_diffusion(value):
        return f"{value:.4e}"

    @staticmethod
    def format_weight(value):
        return f"{value:.6g}"


    def get_msd_fit_conditions(self, grid_vectors, grid_shape):
        """
        Calculate MSD fitting start/end thresholds.

        Parameters
        ----------
        grid_vectors : array-like
            Three grid vectors [a_vec, b_vec, c_vec].
        grid_shape : array-like
            Number of grid points along each direction [Na, Nb, Nc].

        Returns
        -------
        conditions : dict
            Start/end MSD thresholds for a,b,c, x,y,z and isotropic diffusion.
        """

        a_vec = np.asarray(grid_vectors[0], dtype=float)
        b_vec = np.asarray(grid_vectors[1], dtype=float)
        c_vec = np.asarray(grid_vectors[2], dtype=float)

        Na, Nb, Nc = grid_shape

        # Real-space unit cell vectors
        a = a_vec * Na
        b = b_vec * Nb
        c = c_vec * Nc

        # Cell lengths
        cell_lengths = np.array([
            np.linalg.norm(a),
            np.linalg.norm(b),
            np.linalg.norm(c)
        ])

        # Cartesian projections
        xyz_lengths = np.array([
            np.linalg.norm([a[0], b[0], c[0]]),
            np.linalg.norm([a[1], b[1], c[1]]),
            np.linalg.norm([a[2], b[2], c[2]])
        ])

        conditions = {}

        # a,b,c directions
        abc_names = ["a", "b", "c"]

        for name, length in zip(abc_names, cell_lengths):
            conditions[name] = {
                "start": length**2,
                "end": 4 * length**2
            }

        # x,y,z directions
        xyz_names = ["x", "y", "z"]

        for name, length in zip(xyz_names, xyz_lengths):
            conditions[name] = {
                "start": length**2,
                "end": 4 * length**2
            }

        # Isotropic 3D
        diagonal_squared = (
            cell_lengths[0]**2 +
            cell_lengths[1]**2 +
            cell_lengths[2]**2
        )

        conditions["iso"] = {
            "start": diagonal_squared,
            "end": 4 * diagonal_squared
        }

        return conditions