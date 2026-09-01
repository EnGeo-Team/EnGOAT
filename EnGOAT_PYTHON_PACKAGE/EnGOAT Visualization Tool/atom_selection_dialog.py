from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QCheckBox, QDialog, QGridLayout, QGroupBox, QLabel,
    QLayout, QScrollArea, QSizePolicy, QVBoxLayout, QWidget
)
import numpy as np


class AtomsDialog(QDialog):
    PANEL_MAX_HEIGHT = 300
    COLUMN_WIDTHS = {0: 80, 1: 45, 2: 180}

    def __init__(self, project):
        super().__init__()
        self.project = project
        self.setWindowTitle("Atoms")
        self.resize(500, 600)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(4)

        self.atom_widgets = {}
        self.build_atom_panels(layout)
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

        show_all.setChecked(any(visibility[atom_id] for atom_id in atom_ids))
        show_all_label.setChecked(any(self.project.visibility["atom_labels"][atom_id] for atom_id in atom_ids))

        fractional_coordinates = QLabel("Fractional Coordinates")
        self.make_bold(fractional_coordinates)

        layout.addWidget(show_all, 0, 0, Qt.AlignLeft)
        layout.addWidget(show_all_label, 0, 1, Qt.AlignCenter)
        layout.addWidget(fractional_coordinates, 0, 2, Qt.AlignCenter)

        self.set_column_widths(layout, self.COLUMN_WIDTHS)

        show_all.toggled.connect(lambda state: self.on_master_visibility_changed(atom_type, state))
        show_all_label.toggled.connect(lambda state: self.on_master_label_changed(atom_type, state))

        widgets["master_row"] = {"visible": show_all, "label": show_all_label}
        return layout

    # ------------------------------------------------------------------
    # Atom rows
    # ------------------------------------------------------------------

    def build_atom_row(self, atom_id):
        atom = self.project.atoms[atom_id]
        visible = self.project.visibility["individual_atoms"][atom_id]
        label_visible = self.project.visibility["atom_labels"][atom_id]

        visible_checkbox = QCheckBox(atom_id)
        visible_checkbox.setChecked(visible)
        self.make_bold(visible_checkbox)

        label_checkbox = QCheckBox()
        label_checkbox.setChecked(label_visible)

        fractional_coordinates = self.get_fractional_coordinates(atom["center"])
        coordinates_label = QLabel(
            f"({fractional_coordinates[0]:.3f}, "
            f"{fractional_coordinates[1]:.3f}, "
            f"{fractional_coordinates[2]:.3f})"
        )
        coordinates_label.setAlignment(Qt.AlignCenter)

        visible_checkbox.toggled.connect(lambda state: self.on_visibility_changed(atom_id, state))
        label_checkbox.toggled.connect(lambda state: self.on_label_changed(atom_id, state))

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
            visibility[atom_id] = state
            self.atom_widgets[atom_type]["individual_atoms"][atom_id]["visible"].blockSignals(True)
            self.atom_widgets[atom_type]["individual_atoms"][atom_id]["visible"].setChecked(state)
            self.atom_widgets[atom_type]["individual_atoms"][atom_id]["visible"].blockSignals(False)

        self.project.set_visibility("individual_atoms", visibility)

    def on_visibility_changed(self, atom_id, state):
        visibility = self.project.visibility["individual_atoms"].copy()
        visibility[atom_id] = state
        self.project.set_visibility("individual_atoms", visibility)

    def on_master_label_changed(self, atom_type, state):
        atom_ids = self.get_atom_ids(atom_type)
        visibility = self.project.visibility["atom_labels"].copy()

        for atom_id in atom_ids:
            visibility[atom_id] = state
            self.atom_widgets[atom_type]["individual_atoms"][atom_id]["label"].blockSignals(True)
            self.atom_widgets[atom_type]["individual_atoms"][atom_id]["label"].setChecked(state)
            self.atom_widgets[atom_type]["individual_atoms"][atom_id]["label"].blockSignals(False)

        self.project.set_visibility("atom_labels", visibility)

    def on_label_changed(self, atom_id, state):
        visibility = self.project.visibility["atom_labels"].copy()
        visibility[atom_id] = state
        self.project.set_visibility("atom_labels", visibility)

    # ------------------------------------------------------------------
    # Data helpers
    # ------------------------------------------------------------------

    def get_atom_types(self):
        return sorted({atom["type"] for atom in self.project.atoms.values()})

    def get_atom_ids(self, atom_type):
        return [atom_id for atom_id, atom in self.project.atoms.items() if atom["type"] == atom_type]
    
    def get_fractional_coordinates(self, cart_coord):
        grid_shape = self.project.metadata["grid_shape"]
        grid_vectors = self.project.metadata["grid_vectors"]
    
        lattice_vectors = [
            grid_size * np.asarray(grid_vector)
            for grid_size, grid_vector in zip(grid_shape, grid_vectors)
        ]
    
        lattice_matrix = np.column_stack(lattice_vectors)
    
        return np.linalg.solve(lattice_matrix, np.asarray(cart_coord))

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
