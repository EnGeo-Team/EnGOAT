from PySide6.QtWidgets import (
    QGroupBox, 
    QHBoxLayout,
    QSlider, 
    QLabel, 
    QComboBox,
    QRadioButton, 
    QDoubleSpinBox, 
    QScrollArea,
    QWidget,
    QPushButton,
    QVBoxLayout,
    QCheckBox, 
    QColorDialog
)

from PySide6.QtCore import Qt


class ControlPanel(QWidget):

    def __init__(self):
        super().__init__()

        #-----------------------------------------------------
        # Main_window will asign the project to a control panel when build_control_panel(project) is called
        #-----------------------------------------------------
        self.project = None
        
        #-----------------------------------------------------
        # Set style of boxes
        #-----------------------------------------------------
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
        
        #-----------------------------------------------------
        # Create layout
        #-----------------------------------------------------
        main_layout = QVBoxLayout(self)

        # Scroll bar for the inner content
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.setMinimumWidth(330)   # try 280–350 depending on your UI

        # The inner content where everything is added
        content_widget = QWidget()
        self.content_layout = QVBoxLayout(content_widget)
        self.content_layout.setAlignment(Qt.AlignTop)

        # Add the inner content to scroll
        scroll.setWidget(content_widget)

        # Add the scroll to main layout
        main_layout.addWidget(scroll)

        #-----------------------------------------------------
        # Create the LOAD button
        #-----------------------------------------------------
        load_layout = QHBoxLayout()

        self.load_button = QPushButton("Load structure")
        self.load_button.setFixedHeight(40)   # taller button
        self.load_button.setStyleSheet("""QPushButton {font-size: 11px;}""")

        load_layout.addStretch()
        load_layout.addWidget(self.load_button)
        load_layout.addStretch()

        main_layout.addLayout(load_layout)


    def build_control_panel(self, project):

        # Clear old widgets
        while self.content_layout.count():
            item = self.content_layout.takeAt(0)
            if item.widget():
                item.widget().deleteLater()

        # set the project
        self.project = project

        # A dictionary where all buttons are stored
        widgets = {}

        if project is not None:



            # Once the project is set, build all the control buttons and set their states
            widgets["UC_panel"] = self.build_UC_panel()
            widgets["isosurface_panel"] = self.build_isosurface_panel()
            widgets["TuTraSt_panel"] = self.build_TuTraSt_panel()

            self.content_layout.addStretch()

        return widgets

    #-----------------------------------------------------
    # Build the UC section of the control panel
    #-----------------------------------------------------

    def build_UC_panel(self):

        UC_group = QGroupBox("Atoms and Bonds")
        UC_layout = QVBoxLayout()

        UC_widgets = {}

        # UC outline toggle button
        #UC_outline_toggle = QCheckBox("Show unit cell outline")
        #UC_outline_toggle.setChecked(self.project.visibility["UC_outline"])
        #UC_widgets["UC_outline_toggle"] = UC_outline_toggle
        #UC_layout.addWidget(UC_outline_toggle)

        ## UC parameters toggle button
        #UC_parameters_toggle = QCheckBox("Show unit cell parameters")
        #UC_parameters_toggle.setChecked(self.project.visibility["UC_parameters"])
        #UC_widgets["UC_parameters_toggle"] = UC_parameters_toggle
        #UC_layout.addWidget(UC_parameters_toggle)

        # --- ATOMS ROW ---
        atoms_row = QHBoxLayout()

        atoms_toggle = QCheckBox("Display atoms")
        atoms_toggle.setChecked(self.project.visibility["atoms"])
        UC_widgets["atoms_toggle"] = atoms_toggle
        atoms_row.addWidget(atoms_toggle)

        atoms_row.addStretch()

        choose_atoms_button = QPushButton("...")
        choose_atoms_button.setFixedWidth(30)  # ✅ smaller width
        UC_widgets["choose_atoms_button"] = choose_atoms_button
        atoms_row.addWidget(choose_atoms_button)

        UC_layout.addLayout(atoms_row)

        # --- BONDS ROW ---
        bonds_row = QHBoxLayout()

        bonds_toggle = QCheckBox("Display bonds")
        bonds_toggle.setChecked(self.project.visibility["bonds"])
        UC_widgets["bonds_toggle"] = bonds_toggle
        bonds_row.addWidget(bonds_toggle)

        bonds_row.addStretch()

        choose_bonds_button = QPushButton("...")
        choose_bonds_button.setFixedWidth(30)
        UC_widgets["choose_bonds_button"] = choose_bonds_button      
        bonds_row.addWidget(choose_bonds_button)

        UC_layout.addLayout(bonds_row)

        # Add UC layout to UC group box, then add the box to the inner content layout
        UC_group.setLayout(UC_layout)
        self.content_layout.addWidget(UC_group)

        return UC_widgets

    #-----------------------------------------------------
    # Build the isosurface section of the control panel
    #-----------------------------------------------------

    def build_isosurface_panel(self):
        isosurface_group = QGroupBox("Isosurface")
        isosurface_layout = QVBoxLayout()

        isosurface_widgets = {}

        # Isosurface toggle button
        isosurface_toggle = QCheckBox("Show")
        isosurface_toggle.setChecked(self.project.visibility["isosurface"])
        isosurface_widgets["isosurface_toggle"] = isosurface_toggle
        isosurface_layout.addWidget(isosurface_toggle)

        # Isosurface level label + slider
        N_levels = self.project.metadata["N_levels"]
        E_step = self.project.metadata["E_step"]

        isosurface_level_label = QLabel(f"Energy: {self.project.plotting_data["isosurface"]["level"]*E_step:.2f} / {N_levels * E_step:.2f} kJ/mol")
        isosurface_layout.addWidget(isosurface_level_label)
        isosurface_widgets["isosurface_level_label"] = isosurface_level_label

        isosurface_level_slider = QSlider(Qt.Horizontal)
        isosurface_level_slider.setMinimum(1)
        isosurface_level_slider.setMaximum(N_levels)
        isosurface_level_slider.setValue(self.project.plotting_data["isosurface"]["level"])
        isosurface_widgets["isosurface_level_slider"] = isosurface_level_slider
        isosurface_level_slider.setMaximumWidth(200)
        isosurface_layout.addWidget(isosurface_level_slider)

        # --- COLOR AND OPACITY ROW ---
        row, color_button, opacity_slider = self.create_color_opacity_row(self.project.plotting_data["isosurface"]["color"], self.project.plotting_data["isosurface"]["opacity"])

        isosurface_widgets["color_button"] = color_button
        isosurface_widgets["opacity_slider"] = opacity_slider
        isosurface_layout.addLayout(row)

        # Add isosurface layout to isosurface group box, then add the box to the inner content layout
        isosurface_group.setLayout(isosurface_layout)
        self.content_layout.addWidget(isosurface_group)

        return isosurface_widgets

    #-----------------------------------------------------
    # Build the TuTraSt section of the control panel
    #-----------------------------------------------------

    def build_TuTraSt_panel(self):

        TuTraSt_group = QGroupBox("Tunnels and Isolated Groups")
        TuTraSt_layout = QVBoxLayout()

        TuTraSt_widgets = {}

        # View mode / merge tree row
        TuTraSt_row = QHBoxLayout()

        view_mode_dropdown = QComboBox()
        view_mode_dropdown.addItems(["Volumetric view", "Graph view"])
        TuTraSt_widgets["view_mode_dropdown"] = view_mode_dropdown
        TuTraSt_row.addWidget(view_mode_dropdown, stretch=1)

        merge_tree_button = QPushButton("Merge tree")
        merge_tree_button.setFixedWidth(90)
        TuTraSt_widgets["merge_tree_button"] = merge_tree_button
        TuTraSt_row.addWidget(merge_tree_button)

        TuTraSt_layout.addLayout(TuTraSt_row)

        # --- INDIVIDUAL TUNNEL SYSTEMS ---
        tunnel_systems_widgets = {}
        for tunnel_system in self.project.tunnel_systems:
            tunnel_layout, tunnel_widgets = self.build_group_subpanel("tunnel_systems", tunnel_system)
            tunnel_systems_widgets[str(tunnel_system)] = tunnel_widgets
            TuTraSt_layout.addWidget(tunnel_layout)
        
        TuTraSt_widgets["tunnel_systems"] = tunnel_systems_widgets

        # --- INDIVIDUAL ISOLATED GROUPS ---
        isolated_groups_widgets = {}
        for isolated_group in self.project.isolated_groups:
            iso_layout, iso_widgets = self.build_group_subpanel("isolated_groups", isolated_group)
            isolated_groups_widgets[str(isolated_group)] = iso_widgets
            TuTraSt_layout.addWidget(iso_layout)

        TuTraSt_widgets["isolated_groups"] = isolated_groups_widgets

        # Add TuTraSt layout to isosurface group box, then add the box to the inner content layout
        TuTraSt_group.setLayout(TuTraSt_layout)
        self.content_layout.addWidget(TuTraSt_group)

        return TuTraSt_widgets

    def build_group_subpanel(self, group_type, ID):

        if group_type == "tunnel_systems":
            label = f"Tunnel system {ID}"
        elif group_type == "isolated_groups":
            label = f"Isolated group {ID}"

        group_group = QGroupBox(f"{label}")
        group_layout = QVBoxLayout()

        group_widgets = {}

        group_row = QHBoxLayout()

        # Show group row
        group_toggle = QCheckBox(f"Show")
        group_toggle.setChecked(self.project.visibility[group_type][str(ID)])
        group_widgets["toggle"] = group_toggle
        group_row.addWidget(group_toggle)

        group_row.addStretch()

        choose_objects_button = QPushButton("...")
        choose_objects_button.setFixedWidth(30)
        group_widgets["choose_objects_button"] = choose_objects_button
        group_row.addWidget(choose_objects_button)

        group_layout.addLayout(group_row)

        # Color and opacity row
        row, color_button, opacity_slider = self.create_color_opacity_row(self.project.plotting_data[group_type][str(ID)]["color"], self.project.plotting_data[group_type][str(ID)]["opacity"])
        group_widgets["color_button"] = color_button
        group_widgets["opacity_slider"] = opacity_slider
        group_layout.addLayout(row)

        # --- MEP section ---
        if group_type == "tunnel_systems":
            MEP_row = QHBoxLayout()

            # MEP toggle button
            MEP_toggle = QCheckBox("MEP")
            MEP_toggle.setChecked(self.project.visibility["MEPs"][str(ID)])
            group_widgets["MEP_toggle"] = MEP_toggle
            MEP_row.addWidget(MEP_toggle)

            #MEP_row.addStretch()

            # MEP selection radio buttons
            radio_row = QHBoxLayout()
            radios = {}

            for direction, path in self.project.tunnel_systems[str(ID)]["MEPs"].items():
                direction_radio = QRadioButton(direction)

                if path is None:    #Disable non existing directions
                    direction_radio.setEnabled(False)
                if direction == self.project.plotting_data["tunnel_systems"][str(ID)]["MEP"]:    #Set checked the first direction
                    direction_radio.setChecked(True)
                
                radios[direction] = direction_radio
                radio_row.addWidget(direction_radio)
            
            MEP_row.addLayout(radio_row)

            MEP_row.addStretch()

            # Plot E diagram button
            plot_E_diagram_button = QPushButton("Plot")
            plot_E_diagram_button.setFixedWidth(50)
            group_widgets["plot_E_diagram_button"] = plot_E_diagram_button
            MEP_row.addWidget(plot_E_diagram_button)

            group_layout.addLayout(MEP_row)

            group_widgets["MEP_radio_button"] = radios
            #group_layout.addLayout(radio_row)

        group_group.setLayout(group_layout)

        return group_group, group_widgets
    

    @staticmethod
    def create_color_opacity_row(color, opacity):

        row = QHBoxLayout()

        row.addWidget(QLabel("Color:"))
        color_button = QPushButton()
        color_button.setFixedSize(30, 20)
        ControlPanel.set_button_color(color_button, color)

        row.addWidget(color_button)

        row.addSpacing(15)

        row.addWidget(QLabel("Opacity:"))

        opacity_slider = QDoubleSpinBox()
        opacity_slider.setRange(0.0, 1.0)
        opacity_slider.setSingleStep(0.1)
        opacity_slider.setValue(opacity)
        
        row.addWidget(opacity_slider)

        row.addStretch()
        return row, color_button, opacity_slider

    @staticmethod
    def set_button_color(button, color):#maybe move somewhere else idk

        r, g, b = map(
            lambda x: int(x * 255),
            color
        )

        button.setStyleSheet(
            f"border:1px solid #888;"
            f"border-radius:3px;"
            f"background-color:rgb({r},{g},{b});"
        )
    
    @staticmethod
    def choose_color(current_color):

        color = QColorDialog.getColor()

        if not color.isValid():
            return current_color

        return color.getRgbF()[:3]
