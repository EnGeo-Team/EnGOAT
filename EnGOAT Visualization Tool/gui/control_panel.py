from PySide6.QtWidgets import (
    QWidget,
    QPushButton,
    QVBoxLayout,
    QCheckBox
)
from gui.dialogs.atom_selection_dialog import AtomSelectionDialog
from gui.dialogs.bond_selection_dialog import BondSelectionDialog

from PySide6.QtWidgets import QLabel
from PySide6.QtWidgets import QGroupBox, QHBoxLayout
from PySide6.QtWidgets import QSlider, QLineEdit
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QComboBox, QRadioButton

from PySide6.QtWidgets import QScrollArea

class ControlPanel(QWidget):


    def __init__(self):
        super().__init__()

        self.viewer = None
        
        main_layout = QVBoxLayout(self)

        # ✅ Scroll area
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.setMinimumWidth(300)   # try 280–350 depending on your UI


        # ✅ Inner container (this holds all your groups)
        content_widget = QWidget()
        content_layout = QVBoxLayout(content_widget)
        content_layout.setAlignment(Qt.AlignTop)

        scroll.setWidget(content_widget)

        # ✅ Put scroll into main layout
        main_layout.addWidget(scroll)

        self.setStyleSheet("""
            QGroupBox {
                font-weight: bold;
                border: 1px solid #ccc;
                margin-top: 6px;
            }

            QGroupBox::title {
                subcontrol-origin: margin;
                left: 8px;
                padding: 0 3px 0 3px;
            }
        """)

        # ======================
        # 🔹 UNIT CELL GROUP
        # ======================
        uc_group = QGroupBox("Unit cell")
        uc_layout = QVBoxLayout()

        self.unit_cell_toggle = QCheckBox("Show unit cell")
        self.unit_cell_toggle.setChecked(True)

        self.uc_params_toggle = QCheckBox("Show unit cell parameters")
        self.uc_params_toggle.setChecked(True)

        uc_layout.addWidget(self.unit_cell_toggle)
        uc_layout.addWidget(self.uc_params_toggle)

        # --- ATOMS ROW ---
        atoms_row = QHBoxLayout()

        self.atoms_toggle = QCheckBox("Show atoms")
        self.atoms_toggle.setChecked(True)

        self.choose_atoms_button = QPushButton("Choose")
        self.choose_atoms_button.setFixedWidth(70)  # ✅ smaller width
        self.choose_atoms_button.setEnabled(True)
        self.choose_atoms_button.clicked.connect(self.open_atom_dialog)

        atoms_row.addWidget(self.atoms_toggle)
        atoms_row.addStretch()
        atoms_row.addWidget(self.choose_atoms_button)

        uc_layout.addLayout(atoms_row)

        # --- BONDS ROW ---
        bonds_row = QHBoxLayout()

        self.bonds_toggle = QCheckBox("Show bonds")
        self.bonds_toggle.setChecked(True)

        self.choose_bonds_button = QPushButton("Choose")
        self.choose_bonds_button.setFixedWidth(70)
        self.choose_bonds_button.setEnabled(True)
        self.choose_bonds_button.clicked.connect(self.open_bond_dialog)

        bonds_row.addWidget(self.bonds_toggle)
        bonds_row.addStretch()
        bonds_row.addWidget(self.choose_bonds_button)

        uc_layout.addLayout(bonds_row)

        #Main thing
        uc_group.setLayout(uc_layout)
        content_layout.addWidget(uc_group)

        # ======================
        # 🔹 ISOSURFACE GROUP
        # ======================
        iso_group = QGroupBox("Isosurface")
        iso_layout = QVBoxLayout()

        # ✅ toggle
        self.iso_toggle = QCheckBox("Show isosurface")
        self.iso_toggle.setChecked(False)

        iso_layout.addWidget(self.iso_toggle)

        # ✅ slider label
        self.iso_label = QLabel("Energy:   /   kJ/mol")

        # ✅ slider
        self.iso_slider = QSlider(Qt.Horizontal)
        self.iso_slider.setMinimum(1)
        self.iso_slider.setMaximum(10)   # default, update later
        self.iso_slider.setValue(1)

        self.energy_step = 1.0
        self.iso_slider.valueChanged.connect(self.update_iso_label)

        iso_layout.addWidget(self.iso_label)
        iso_layout.addWidget(self.iso_slider)


        #color and opacity
        row = QHBoxLayout()

        # --- COLOR ---
        row.addWidget(QLabel("Color:"))

        self.color_button = QPushButton()
        self.color_button.setFixedSize(30, 20)  # ✅ small square
        self.color_button.setStyleSheet("""
            border: 1px solid #888;
            border-radius: 3px;
            background-color: rgb(135, 206, 235);
            """)
        self.iso_color = (0.529, 0.808, 0.922)

        row.addWidget(self.color_button)

        # spacing
        row.addSpacing(15)

        # --- OPACITY ---
        row.addWidget(QLabel("Opacity:"))

        self.opacity_input = QLineEdit("0.5")
        self.opacity_input.setFixedWidth(50)  # ✅ narrow box

        row.addWidget(self.opacity_input)

        row.addStretch()

        iso_layout.addLayout(row)


        # ✅ plot button
        self.plot_iso_button = QPushButton("Plot isosurface")

        iso_layout.addWidget(self.plot_iso_button)

        iso_group.setLayout(iso_layout)
        content_layout.addWidget(iso_group)

        # initialize disabled state
        self.iso_slider.setEnabled(False)
        self.opacity_input.setEnabled(False)
        self.color_button.setEnabled(False)
        self.plot_iso_button.setEnabled(False)


        # ======================
        # 🔹 BASINS & TUNNELS
        # ======================
        bt_group = QGroupBox("Basins and Tunnels")
        bt_layout = QVBoxLayout()

        # ✅ dropdown
        self.view_mode_combo = QComboBox()
        self.view_mode_combo.addItems(["Volumetric view", "Graph view"])

        bt_layout.addWidget(self.view_mode_combo)

        
        # ======================
        # 🔹 LABELS ROW
        # ======================
        labels_row = QHBoxLayout()

        labels_row.addWidget(QLabel("Show labels:"))
        labels_row.addStretch()

        self.labels_basins_cb = QCheckBox("Basins")
        self.labels_transitions_cb = QCheckBox("Transitions")

        labels_row.addWidget(self.labels_basins_cb)
        labels_row.addWidget(self.labels_transitions_cb)

        bt_layout.addLayout(labels_row)


        # ======================
        # 🔹 ENERGIES ROW
        # ======================
        energies_row = QHBoxLayout()

        energies_row.addWidget(QLabel("Show energies:"))
        energies_row.addStretch()

        self.energies_basins_cb = QCheckBox("Basins")
        self.energies_transitions_cb = QCheckBox("Transitions")

        energies_row.addWidget(self.energies_basins_cb)
        energies_row.addWidget(self.energies_transitions_cb)

        bt_layout.addLayout(energies_row)



        bt_group.setLayout(bt_layout)
        content_layout.addWidget(bt_group)

        self.tunnel_container_layout = QVBoxLayout()
        bt_layout.addLayout(self.tunnel_container_layout)


        # ======================
        # 🔹 SPACER
        # ======================
        content_layout.addStretch()

        # ======================
        # 🔹 LOAD BUTTON (BOTTOM CENTER)
        # ======================
        load_layout = QHBoxLayout()

        self.load_button = QPushButton("Load structure")

        self.load_button.setFixedHeight(40)   # taller button

        self.load_button.setStyleSheet("""
            QPushButton {
                font-size: 11px;
            }
        """)

        load_layout.addStretch()
        load_layout.addWidget(self.load_button)
        load_layout.addStretch()

        main_layout.addLayout(load_layout)

#
#ATOMS AND BONDS
#
    def open_atom_dialog(self):

        if self.viewer is None:
            return

        if not hasattr(self.viewer, "atom_actors"):
            return

        dialog = AtomSelectionDialog(
            elements=self.viewer.atom_actors.keys(),
            selected_elements=self.viewer.visible_elements
        )

        if dialog.exec():
            selected = dialog.get_selection()
            self.viewer.update_visible_atoms(selected)

    def open_bond_dialog(self):

        if self.viewer is None:
            return

        if not hasattr(self.viewer, "bond_actors"):
            return

        dialog = BondSelectionDialog(
            bonds=self.viewer.bond_actors.keys(),
            selected_bonds=self.viewer.visible_bonds
        )

        if dialog.exec():
            selected = dialog.get_selection()
            self.viewer.update_visible_bonds(selected)

#
#ISOSURFACE
#
    def get_opacity(self):
        try:
            value = float(self.opacity_input.text())
            return max(0.0, min(1.0, value))
        except:
            return 0.5

    def update_iso_label(self, val):
        max_val = self.iso_slider.maximum()
        self.iso_label.setText(
            f"Energy: {val * self.energy_step:.2f} / {max_val * self.energy_step:.2f} kJ/mol"
        )

    def update_iso_range(self, max_level, energy_step):
        self.energy_step = energy_step
        self.iso_slider.setMaximum(max_level)
        self.update_iso_label(self.iso_slider.value())

    def update_tunnel_systems(self, tunnel_layout_dict, isolated_clusters, project):

        # ✅ clear old UI
        while self.tunnel_container_layout.count():
            item = self.tunnel_container_layout.takeAt(0)
            widget = item.widget()
            if widget:
                widget.deleteLater()

        self.tunnel_widgets = {}

        for name, meps in tunnel_layout_dict.items():

            group = QGroupBox(name)
            layout = QVBoxLayout()

            # --- Show toggle ---
            show_cb = QCheckBox(f"Show {name}")
            layout.addWidget(show_cb)

            show_cb.blockSignals(True)
            show_cb.setChecked(
                project.tunnel_systems_plotting[name]["basins"]["visible"]
            )
            show_cb.blockSignals(False)

            # --- Color + opacity ---
            row = QHBoxLayout()

            color_btn = QPushButton()

            color_btn.setFixedSize(25, 15)

            color = project.tunnel_systems_plotting[name]["basins"]["color"]
            opacity = project.tunnel_systems_plotting[name]["basins"]["opacity"]

            r, g, b = color

            color_btn = QPushButton()
            color_btn.setFixedSize(25, 15)

            color_btn.setStyleSheet(f"""
                border: 1px solid #888;
                border-radius: 3px;
                background-color: rgb(
                    {int(r*255)},
                    {int(g*255)},
                    {int(b*255)}
                );
            """)

            # Keep color available for plot button
            color_btn.selected_color = color

            opacity_input = QLineEdit(str(opacity))
            opacity_input.setFixedWidth(50)

            row.addWidget(QLabel("Color:"))
            row.addWidget(color_btn)
            row.addSpacing(10)
            row.addWidget(QLabel("Opacity:"))
            row.addWidget(opacity_input)
            row.addStretch()

            layout.addLayout(row)

            # --- Plot button ---
            plot_btn = QPushButton(f"Plot {name}")
            layout.addWidget(plot_btn)
            
            # --- MEP toggle + E diagram button ---
            mep_row_top = QHBoxLayout()

            mep_cb = QCheckBox("Show MEPs")

            mep_cb.blockSignals(True)
            mep_cb.setChecked(
                project.tunnel_systems_plotting[name]["show_MEP"]
            )
            mep_cb.blockSignals(False)

            plot_E_btn = QPushButton("Plot E diagram")
            plot_E_btn.setFixedWidth(110)  # optional, keeps it compact

            mep_row_top.addWidget(mep_cb)
            mep_row_top.addStretch()
            mep_row_top.addWidget(plot_E_btn)

            layout.addLayout(mep_row_top)

            # --- MEP radios ---
            mep_row = QHBoxLayout()
            mep_buttons = {}

            saved_direction = project.tunnel_systems_plotting[name]["MEP_direction"]
            first_available = None

            for key in ["a", "b", "c"]:
            
                rb = QRadioButton(key)

                if meps[key] == False:
                    rb.setEnabled(False)
                else:
                
                    if first_available is None:
                        first_available = key

                mep_row.addWidget(rb)

                mep_buttons[key] = rb

            if (
                saved_direction in mep_buttons
                and meps[saved_direction] is not None
            ):
                mep_buttons[saved_direction].setChecked(True)

            elif first_available is not None:
                mep_buttons[first_available].setChecked(True)

            mep_row.addStretch()
            layout.addLayout(mep_row)

            # --- Finalize group ---
            group.setLayout(layout)
            self.tunnel_container_layout.addWidget(group)

            # store references (important later)
            self.tunnel_widgets[name] = {
                "show": show_cb,
                "color": color_btn,
                "opacity": opacity_input,
                "plot": plot_btn,
                "plot_E": plot_E_btn,
                "mep_toggle": mep_cb,
                "mep_buttons": mep_buttons,
            }

        # ==========================
        # ✅ ISOLATED CLUSTERS ENTRY
        # ==========================
        if isolated_clusters:
        
            group = QGroupBox("Isolated clusters")
            layout = QVBoxLayout()

            # --- Show toggle ---
            show_cb = QCheckBox("Show isolated clusters")
            show_cb.blockSignals(True)
            show_cb.setChecked(
                project.isolated_clusters_plotting["basins"]["visible"]
            )
            show_cb.blockSignals(False)
            layout.addWidget(show_cb)

            # --- Color + opacity ---
            row = QHBoxLayout()

            color_btn = QPushButton()
            color_btn.setFixedSize(25, 15)
            color = project.isolated_clusters_plotting["basins"]["color"]
            opacity = project.isolated_clusters_plotting["basins"]["opacity"]

            r, g, b = color

            color_btn = QPushButton()
            color_btn.setFixedSize(25, 15)

            color_btn.setStyleSheet(f"""
                border: 1px solid #888;
                border-radius: 3px;
                background-color: rgb(
                    {int(r*255)},
                    {int(g*255)},
                    {int(b*255)}
                );
            """)

            color_btn.selected_color = color

            opacity_input = QLineEdit(str(opacity))
            opacity_input.setFixedWidth(50)

            row.addWidget(QLabel("Color:"))
            row.addWidget(color_btn)
            row.addSpacing(10)
            row.addWidget(QLabel("Opacity:"))
            row.addWidget(opacity_input)
            row.addStretch()

            layout.addLayout(row)

            # --- Plot button ---
            plot_btn = QPushButton("Plot isolated clusters")
            layout.addWidget(plot_btn)

            group.setLayout(layout)
            self.tunnel_container_layout.addWidget(group)

            # ✅ store separately (recommended)
            self.isolated_cluster_widgets = {
                "show": show_cb,
                "color": color_btn,
                "opacity": opacity_input,
                "plot": plot_btn
            }








