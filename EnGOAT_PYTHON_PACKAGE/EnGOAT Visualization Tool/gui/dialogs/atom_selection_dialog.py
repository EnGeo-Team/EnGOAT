from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QCheckBox, QDialog, QDoubleSpinBox, QGridLayout, QGroupBox, QLabel,
    QPushButton, QLayout, QScrollArea, QSizePolicy, QVBoxLayout, QWidget
)
import numpy as np


class AtomsDialog(QDialog):
    PANEL_MAX_HEIGHT = 300
    COLUMN_WIDTHS = {0: 80, 1: 45, 2: 180}

    def __init__(self, project):
        super().__init__()

        self.project = project
        self.project.selected_basin_changed.connect(self.update_selection)
        self.setWindowTitle("Atoms")
        self.resize(500, 600)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(4)

        self.atom_widgets = {}

        self.build_atom_panels(layout)
        self.build_basin_neighbourhood_panel(layout)

        layout.addStretch(1)

    # ------------------------------------------------------------------
    # Panel construction
    # ------------------------------------------------------------------

    def build_atom_panels(self, parent_layout):
        for atom_type in self.get_atom_types():
            self.atom_widgets[atom_type] = self.build_atom_panel(parent_layout, atom_type)

    def build_atom_panel(self, parent_layout, atom_type):
        atom_ids = self.get_atom_ids(atom_type)
        widgets = {"master_row": {}, "individual_atoms": {}}

        group_box = QGroupBox(atom_type)
        group_layout = QVBoxLayout(group_box)
        group_layout.setContentsMargins(6, 6, 6, 6)
        group_layout.setSpacing(4)
        group_layout.setSizeConstraint(QLayout.SetMinimumSize)
        group_box.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Maximum)

        group_layout.addLayout(self.build_master_row(atom_type, atom_ids, widgets))

        rows_layout = QGridLayout()
        rows_layout.setContentsMargins(8, 0, 8, 4)
        rows_layout.setHorizontalSpacing(2)
        rows_layout.setVerticalSpacing(3)

        for row, atom_id in enumerate(atom_ids):
            row_widgets = self.build_atom_row(atom_id)
            self.add_atom_row(rows_layout, row, row_widgets)
            widgets["individual_atoms"][atom_id] = row_widgets

        self.set_column_widths(rows_layout, self.COLUMN_WIDTHS)

        content_widget = QWidget()
        content_widget.setLayout(rows_layout)
        group_layout.addWidget(self.create_scroll_area(content_widget, self.PANEL_MAX_HEIGHT))

        parent_layout.addWidget(group_box)
        return widgets

    # ------------------------------------------------------------------
    # Select Basin Neighbourhood
    # ------------------------------------------------------------------

    def build_basin_neighbourhood_panel(self, parent_layout):
        group_box = QGroupBox("Select Basin Neighbourhood")
        layout = QGridLayout(group_box)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setHorizontalSpacing(8)
        layout.setVerticalSpacing(6)

        self.selected_basin_label = QLabel()
        layout.addWidget(self.selected_basin_label, 0, 0, 1, 3)

        radius_label = QLabel("R:")

        self.radius_spinbox = QDoubleSpinBox()
        self.radius_spinbox.setRange(0.0, 1e6)
        self.radius_spinbox.setDecimals(3)
        self.radius_spinbox.setSingleStep(0.1)
        self.radius_spinbox.setValue(self.project.select_radius)
        self.radius_spinbox.setSuffix(" Å")
        self.radius_spinbox.valueChanged.connect(self.on_select_radius_changed)

        layout.addWidget(radius_label, 1, 0, Qt.AlignRight)
        layout.addWidget(self.radius_spinbox, 1, 1)

        self.select_neighbourhood_button = QPushButton("Select neighbourhood")
        self.select_neighbourhood_button.clicked.connect(self.select_basin_neighbourhood)
        layout.addWidget(self.select_neighbourhood_button, 1, 2)

        self.update_selection()
        parent_layout.addWidget(group_box)

    def on_select_radius_changed(self, value):
        self.project.select_radius = float(value)

    def update_selection(self):
        selected_basin = self.project.selected_basin

        self.selected_basin_label.setText(f"Selected basin: {selected_basin}")

        self.select_neighbourhood_button.setEnabled(
            selected_basin is not None and selected_basin in self.project.basin_data
        )

    def select_basin_neighbourhood(self):
        selected_basin = self.project.selected_basin

        if selected_basin is None or selected_basin not in self.project.basin_data:
            return

        radius = float(self.project.select_radius)
        basin_center = self.get_basin_center_cartesian(selected_basin)
        visibility = self.project.visibility["individual_atoms"].copy()

        for atom_id, atom in self.project.atoms.items():
            atom_center = np.asarray(atom["center"], dtype=float)
            distance = self.get_pbc_distance(atom_center, basin_center)
            visibility[atom_id] = bool(distance <= radius)

        self.project.set_visibility("individual_atoms", visibility)
        self.update_atom_checkboxes(visibility)

    def update_atom_checkboxes(self, visibility):
        for atom_type, widgets in self.atom_widgets.items():
            atom_ids = self.get_atom_ids(atom_type)

            for atom_id in atom_ids:
                checkbox = widgets["individual_atoms"][atom_id]["visible"]
                checkbox.blockSignals(True)
                checkbox.setChecked(bool(visibility[atom_id]))
                checkbox.blockSignals(False)

            master_checkbox = widgets["master_row"]["visible"]
            states = [bool(visibility[atom_id]) for atom_id in atom_ids]

            master_checkbox.blockSignals(True)
            master_checkbox.setChecked(bool(states and any(states)))
            master_checkbox.blockSignals(False)

    # ------------------------------------------------------------------
    # Basin / PBC helpers
    # ------------------------------------------------------------------

    def get_basin_center_cartesian(self, basin_id):
        """
        Convert a basin center from grid indices to Cartesian coordinates.

        basin_data[basin_id]["center"] contains grid indices [ix, iy, iz].
        """
        grid_center = np.asarray(self.project.basin_data[basin_id]["center"], dtype=float)
        origin = np.asarray(self.project.metadata["origin"], dtype=float)

        grid_vectors = [
            np.asarray(vector, dtype=float)
            for vector in self.project.metadata["grid_vectors"]
        ]

        return origin + sum(
            index * vector
            for index, vector in zip(grid_center, grid_vectors)
        )

    def get_lattice_matrix(self):
        grid_shape = self.project.metadata["grid_shape"]
        grid_vectors = self.project.metadata["grid_vectors"]

        lattice_vectors = [
            grid_size * np.asarray(grid_vector, dtype=float)
            for grid_size, grid_vector in zip(grid_shape, grid_vectors)
        ]

        return np.column_stack(lattice_vectors)

    def get_pbc_distance(self, coord1, coord2):
        lattice_matrix = self.get_lattice_matrix()

        coord1 = np.asarray(coord1, dtype=float)
        coord2 = np.asarray(coord2, dtype=float)

        fractional_delta = np.linalg.solve(lattice_matrix, coord1 - coord2)
        fractional_delta -= np.round(fractional_delta)

        cartesian_delta = lattice_matrix @ fractional_delta

        return float(np.linalg.norm(cartesian_delta))

    # ------------------------------------------------------------------
    # Master row
    # ------------------------------------------------------------------

    def build_master_row(self, atom_type, atom_ids, widgets):
        layout = QGridLayout()
        layout.setContentsMargins(8, 2, 8, 2)
        layout.setHorizontalSpacing(4)
        layout.setVerticalSpacing(0)

        show_all = QCheckBox("Show:")
        show_all_label = QCheckBox("Label")

        self.make_bold(show_all)
        self.make_bold(show_all_label)

        visibility = self.project.visibility["individual_atoms"]

        # Original behaviour: checked if ANY atom of this type is visible.
        show_all.setChecked(any(visibility[atom_id] for atom_id in atom_ids))

        label_visibility = self.project.visibility["atom_labels"]
        show_all_label.setChecked(
            any(label_visibility[atom_id] for atom_id in atom_ids)
        )

        fractional_coordinates = QLabel("Fractional Coordinates")
        self.make_bold(fractional_coordinates)

        layout.addWidget(show_all, 0, 0, Qt.AlignLeft)
        layout.addWidget(show_all_label, 0, 1, Qt.AlignCenter)
        layout.addWidget(fractional_coordinates, 0, 2, Qt.AlignCenter)

        self.set_column_widths(layout, self.COLUMN_WIDTHS)

        show_all.toggled.connect(
            lambda state: self.on_master_visibility_changed(atom_type, state)
        )

        show_all_label.toggled.connect(
            lambda state: self.on_master_label_changed(atom_type, state)
        )

        widgets["master_row"] = {
            "visible": show_all,
            "label": show_all_label
        }

        return layout

    # ------------------------------------------------------------------
    # Atom rows
    # ------------------------------------------------------------------

    def build_atom_row(self, atom_id):
        atom = self.project.atoms[atom_id]
        visible = self.project.visibility["individual_atoms"][atom_id]
        label_visible = self.project.visibility["atom_labels"][atom_id]

        visible_checkbox = QCheckBox(atom_id)
        visible_checkbox.setChecked(bool(visible))
        self.make_bold(visible_checkbox)

        label_checkbox = QCheckBox()
        label_checkbox.setChecked(bool(label_visible))

        fractional_coordinates = self.get_fractional_coordinates(atom["center"])

        coordinates_label = QLabel(
            f"({fractional_coordinates[0]:.3f}, "
            f"{fractional_coordinates[1]:.3f}, "
            f"{fractional_coordinates[2]:.3f})"
        )
        coordinates_label.setAlignment(Qt.AlignCenter)

        visible_checkbox.toggled.connect(
            lambda state: self.on_visibility_changed(atom_id, state)
        )

        label_checkbox.toggled.connect(
            lambda state: self.on_label_changed(atom_id, state)
        )

        return {
            "visible": visible_checkbox,
            "label": label_checkbox,
            "coordinates": coordinates_label,
        }

    def add_atom_row(self, layout, row, widgets):
        layout.addWidget(widgets["visible"], row, 0, Qt.AlignLeft)
        layout.addWidget(widgets["label"], row, 1, Qt.AlignCenter)
        layout.addWidget(widgets["coordinates"], row, 2, Qt.AlignCenter)

    # ------------------------------------------------------------------
    # Visibility and labels
    # ------------------------------------------------------------------

    def on_master_visibility_changed(self, atom_type, state):
        atom_ids = self.get_atom_ids(atom_type)
        visibility = self.project.visibility["individual_atoms"].copy()

        for atom_id in atom_ids:
            visibility[atom_id] = bool(state)

            checkbox = (
                self.atom_widgets[atom_type]
                ["individual_atoms"][atom_id]
                ["visible"]
            )

            checkbox.blockSignals(True)
            checkbox.setChecked(bool(state))
            checkbox.blockSignals(False)

        self.project.set_visibility("individual_atoms", visibility)

    def on_visibility_changed(self, atom_id, state):
        visibility = self.project.visibility["individual_atoms"].copy()
        visibility[atom_id] = bool(state)

        self.project.set_visibility("individual_atoms", visibility)
        self.update_master_visibility_checkbox(atom_id)

    def update_master_visibility_checkbox(self, atom_id):
        atom_type = self.project.atoms[atom_id]["type"]
        atom_ids = self.get_atom_ids(atom_type)
        visibility = self.project.visibility["individual_atoms"]

        master_checkbox = self.atom_widgets[atom_type]["master_row"]["visible"]

        master_checkbox.blockSignals(True)
        master_checkbox.setChecked(
            any(visibility[atom_id] for atom_id in atom_ids)
        )
        master_checkbox.blockSignals(False)

    def on_master_label_changed(self, atom_type, state):
        atom_ids = self.get_atom_ids(atom_type)
        visibility = self.project.visibility["atom_labels"].copy()

        for atom_id in atom_ids:
            visibility[atom_id] = bool(state)

            checkbox = (
                self.atom_widgets[atom_type]
                ["individual_atoms"][atom_id]
                ["label"]
            )

            checkbox.blockSignals(True)
            checkbox.setChecked(bool(state))
            checkbox.blockSignals(False)

        self.project.set_visibility("atom_labels", visibility)

    def on_label_changed(self, atom_id, state):
        visibility = self.project.visibility["atom_labels"].copy()
        visibility[atom_id] = bool(state)

        self.project.set_visibility("atom_labels", visibility)

    # ------------------------------------------------------------------
    # Data helpers
    # ------------------------------------------------------------------

    def get_atom_types(self):
        return sorted({atom["type"] for atom in self.project.atoms.values()})

    def get_atom_ids(self, atom_type):
        return [
            atom_id
            for atom_id, atom in self.project.atoms.items()
            if atom["type"] == atom_type
        ]

    def get_fractional_coordinates(self, cart_coord):
        lattice_matrix = self.get_lattice_matrix()
        return np.linalg.solve(
            lattice_matrix,
            np.asarray(cart_coord, dtype=float)
        )

    # ------------------------------------------------------------------
    # Widget helpers
    # ------------------------------------------------------------------

    @staticmethod
    def create_scroll_area(widget, max_height=300):
        scroll_area = QScrollArea()
        scroll_area.setWidget(widget)
        scroll_area.setWidgetResizable(True)
        scroll_area.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        scroll_area.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        scroll_area.setMaximumHeight(max_height)
        return scroll_area

    @staticmethod
    def set_column_widths(layout, column_widths):
        for column, width in column_widths.items():
            layout.setColumnMinimumWidth(column, width)
            layout.setColumnStretch(column, 0)

    @staticmethod
    def make_bold(widget):
        font = widget.font()
        font.setBold(True)
        widget.setFont(font)