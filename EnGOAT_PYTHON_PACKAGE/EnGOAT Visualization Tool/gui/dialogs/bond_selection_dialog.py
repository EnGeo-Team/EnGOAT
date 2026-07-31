from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QCheckBox, QDialog, QGridLayout, QGroupBox, QLabel,
    QLayout, QScrollArea, QSizePolicy, QVBoxLayout, QWidget
)


class BondsDialog(QDialog):
    PANEL_MAX_HEIGHT = 300
    COLUMN_WIDTHS = {0: 100, 1: 45, 2: 180, 3: 80}

    def __init__(self, project):
        super().__init__()
        self.project = project
        self.setWindowTitle("Bonds")
        self.resize(550, 600)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(4)

        self.bond_widgets = {}
        self.build_bond_panels(layout)

        layout.addStretch(1)

    # ------------------------------------------------------------------
    # Panel construction
    # ------------------------------------------------------------------

    def build_bond_panels(self, parent_layout):
        for bond_type in self.get_bond_types():
            self.bond_widgets[bond_type] = self.build_bond_panel(parent_layout, bond_type)

    def build_bond_panel(self, parent_layout, bond_type):
        bond_ids = self.get_bond_ids(bond_type)
        widgets = {"master_row": {}, "individual_bonds": {}}

        group_box = QGroupBox(str(bond_type))
        group_layout = QVBoxLayout(group_box)
        group_layout.setContentsMargins(6, 6, 6, 6)
        group_layout.setSpacing(4)
        group_layout.setSizeConstraint(QLayout.SetMinimumSize)
        group_box.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Maximum)

        group_layout.addLayout(self.build_master_row(bond_type, bond_ids, widgets))

        rows_layout = QGridLayout()
        rows_layout.setContentsMargins(8, 0, 8, 4)
        rows_layout.setHorizontalSpacing(2)
        rows_layout.setVerticalSpacing(3)

        for row, bond_id in enumerate(bond_ids):
            row_widgets = self.build_bond_row(bond_id)
            self.add_bond_row(rows_layout, row, row_widgets)
            widgets["individual_bonds"][bond_id] = row_widgets

        self.set_column_widths(rows_layout, self.COLUMN_WIDTHS)

        content_widget = QWidget()
        content_widget.setLayout(rows_layout)
        group_layout.addWidget(self.create_scroll_area(content_widget, self.PANEL_MAX_HEIGHT))

        parent_layout.addWidget(group_box)
        return widgets

    # ------------------------------------------------------------------
    # Master row
    # ------------------------------------------------------------------

    def build_master_row(self, bond_type, bond_ids, widgets):
        layout = QGridLayout()
        layout.setContentsMargins(8, 2, 8, 2)
        layout.setHorizontalSpacing(4)
        layout.setVerticalSpacing(0)

        show_all = QCheckBox("Show:")
        show_all_label = QCheckBox("Label")

        self.make_bold(show_all)
        self.make_bold(show_all_label)

        visibility = self.project.visibility["individual_bonds"]
        labels = self.project.visibility["bond_labels"]

        visible_ids = [
            bond_id for bond_id in bond_ids
            if bond_id not in self.project.invisible_bonds
        ]

        show_all.setChecked(any(visibility[bond_id] for bond_id in visible_ids))
        show_all_label.setChecked(any(labels[bond_id] for bond_id in visible_ids))

        connecting_atoms = QLabel("Connecting atoms")
        distance = QLabel("Distance [Å]")

        self.make_bold(connecting_atoms)
        self.make_bold(distance)

        layout.addWidget(show_all, 0, 0, Qt.AlignLeft)
        layout.addWidget(show_all_label, 0, 1, Qt.AlignCenter)
        layout.addWidget(connecting_atoms, 0, 2, Qt.AlignCenter)
        layout.addWidget(distance, 0, 3, Qt.AlignCenter)

        self.set_column_widths(layout, self.COLUMN_WIDTHS)

        show_all.toggled.connect(
            lambda state: self.on_master_visibility_changed(bond_type, state)
        )

        show_all_label.toggled.connect(
            lambda state: self.on_master_label_changed(bond_type, state)
        )

        widgets["master_row"] = {
            "visible": show_all,
            "label": show_all_label,
        }

        return layout

    # ------------------------------------------------------------------
    # Bond rows
    # ------------------------------------------------------------------

    def build_bond_row(self, bond_id):
        bond = self.project.bonds[bond_id]
        invisible = bond_id in self.project.invisible_bonds

        visible = (
            self.project.visibility["individual_bonds"][bond_id]
            if not invisible else False
        )

        label_visible = (
            self.project.visibility["bond_labels"][bond_id]
            if not invisible else False
        )

        visible_checkbox = QCheckBox(bond_id)
        visible_checkbox.setChecked(visible)
        self.make_bold(visible_checkbox)

        label_checkbox = QCheckBox()
        label_checkbox.setChecked(label_visible)

        atom1, atom2 = bond["atoms"]

        atoms_label = QLabel(f"{atom1}, {atom2}")
        atoms_label.setAlignment(Qt.AlignCenter)

        distance_label = QLabel(f"{bond['distance']:.3f}")
        distance_label.setAlignment(Qt.AlignCenter)

        if invisible:
            visible_checkbox.setEnabled(False)
            label_checkbox.setEnabled(False)
            visible_checkbox.setStyleSheet("color: gray; font-weight: bold;")
            atoms_label.setStyleSheet("color: gray;")
            distance_label.setStyleSheet("color: gray;")
        else:
            visible_checkbox.toggled.connect(
                lambda state: self.on_visibility_changed(bond_id, state)
            )

            label_checkbox.toggled.connect(
                lambda state: self.on_label_changed(bond_id, state)
            )

        return {
            "visible": visible_checkbox,
            "label": label_checkbox,
            "atoms": atoms_label,
            "distance": distance_label,
        }

    def add_bond_row(self, layout, row, widgets):
        layout.addWidget(widgets["visible"], row, 0, Qt.AlignLeft)
        layout.addWidget(widgets["label"], row, 1, Qt.AlignCenter)
        layout.addWidget(widgets["atoms"], row, 2, Qt.AlignCenter)
        layout.addWidget(widgets["distance"], row, 3, Qt.AlignCenter)

    # ------------------------------------------------------------------
    # Visibility and labels
    # ------------------------------------------------------------------

    def on_master_visibility_changed(self, bond_type, state):
        bond_ids = self.get_bond_ids(bond_type)
        visibility = self.project.visibility["individual_bonds"].copy()

        for bond_id in bond_ids:
            if bond_id in self.project.invisible_bonds:
                visibility[bond_id] = False
                continue

            visibility[bond_id] = state
            self.bond_widgets[bond_type]["individual_bonds"][bond_id]["visible"].blockSignals(True)
            self.bond_widgets[bond_type]["individual_bonds"][bond_id]["visible"].setChecked(state)
            self.bond_widgets[bond_type]["individual_bonds"][bond_id]["visible"].blockSignals(False)

        self.project.set_visibility("individual_bonds", visibility)

    def on_visibility_changed(self, bond_id, state):
        visibility = self.project.visibility["individual_bonds"].copy()
        visibility[bond_id] = state
        self.project.set_visibility("individual_bonds", visibility)

    def on_master_label_changed(self, bond_type, state):
        bond_ids = self.get_bond_ids(bond_type)
        visibility = self.project.visibility["bond_labels"].copy()

        for bond_id in bond_ids:
            if bond_id in self.project.invisible_bonds:
                visibility[bond_id] = False
                continue

            visibility[bond_id] = state
            self.bond_widgets[bond_type]["individual_bonds"][bond_id]["label"].blockSignals(True)
            self.bond_widgets[bond_type]["individual_bonds"][bond_id]["label"].setChecked(state)
            self.bond_widgets[bond_type]["individual_bonds"][bond_id]["label"].blockSignals(False)

        self.project.set_visibility("bond_labels", visibility)

    def on_label_changed(self, bond_id, state):
        visibility = self.project.visibility["bond_labels"].copy()
        visibility[bond_id] = state
        self.project.set_visibility("bond_labels", visibility)

    # ------------------------------------------------------------------
    # Data helpers
    # ------------------------------------------------------------------

    def get_bond_types(self):
        return sorted(
            {bond["type"] for bond in self.project.bonds.values()},
            key=str,
        )

    def get_bond_ids(self, bond_type):
        return [
            bond_id
            for bond_id, bond in self.project.bonds.items()
            if bond["type"] == bond_type
        ]

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