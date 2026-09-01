from PySide6.QtWidgets import (
    QMainWindow,
    QWidget,
    QHBoxLayout,
    QFileDialog,
    QTabWidget
)
from PySide6.QtGui import QAction
from PySide6.QtCore import Qt

import os
from pathlib import Path

from gui.control_panel import ControlPanel
from gui.project import Project

from gui.dialogs.atom_selection_dialog import AtomsDialog
from gui.dialogs.bond_selection_dialog import BondsDialog
from gui.dialogs.group_dialog import GroupDialog
from gui.dialogs.diffusion_dialog import DiffusionDialog
from gui.dialogs.arrhenius_dialog import ArrheniusDialog
from gui.dialogs.supercell_dialog import SupercellDialog

from plotting.matplotlib.plot_merge_tree import create_merge_trees
from plotting.matplotlib.plot_MEP import plot_MEP_energy_diagram

class MainWindow(QMainWindow):

    def __init__(self):
        super().__init__()

        # --- MAIN WINDOW LAYOUT ---
        self.setWindowTitle("EnGOAT_Visualization")
        self.resize(1400, 900)
        central = QWidget()
        self.setCentralWidget(central)
        layout = QHBoxLayout(central)

        # Open dialogue windows - Enables rendering when multiple windows are open
        self.open_dialogs = {}

        # --- SET UP CONTROL PANEL ---
        self.controls = ControlPanel()
        layout.addWidget(self.controls)
        # Connect the LOAD button
        self.controls.load_button.clicked.connect(self.load_file)

        # --- SET UP DROPDOWN MENUS ---
        self.create_menus()

        # --- SET UP TABS ---
        self.tabs = QTabWidget()
        layout.addWidget(self.tabs, stretch=1)
        # Sync UI when changing tabs
        self.tabs.currentChanged.connect(self.on_tab_changed)
        # Closing tabs
        self.tabs.setTabsClosable(True)
        self.tabs.tabCloseRequested.connect(self.close_tab)

    #-----------------------------------------------------
    # FILE LOADING FUNCTION
    #-----------------------------------------------------

    def load_file(self):

        filename, _ = QFileDialog.getOpenFileName(
            self,
            "Open cube file",
            "",
            "Cube files (*.cube);;All files (*)"
        )

        print("Loading:", filename)
        directory = Path(filename).parent

        # Create a project, build visualization, connect it to a tab
        project = Project()
        project.load_structure(directory)
        index = self.tabs.addTab(project, os.path.basename(filename))
        self.tabs.setCurrentIndex(index)
    
        # Create and connect the control panel
        self.create_control_panel()
        self.update_menus()

        return None

    #-----------------------------------------------------
    # SCREENSHOT FUNCTION
    #-----------------------------------------------------

    def save_screenshot(self):
        project = self.current_project()
        filename, _ = QFileDialog.getSaveFileName(
            self,
            "Save screenshot",
            "screenshot.png",
            "PNG Images (*.png)"
        )
        if not filename:
            return
        project.viewer.plotter.screenshot(filename)
        return None

    #-----------------------------------------------------
    # TABS FUNCTIONS
    #-----------------------------------------------------

    # Simple tab closer
    def close_tab(self, index):
        widget = self.tabs.widget(index)
        widget.deleteLater()
        self.tabs.removeTab(index)

        return None

    # Tells you which project you are in 
    def current_project(self):
        return self.tabs.currentWidget()
    
    def on_tab_changed(self):
        self.create_control_panel()
        self.update_menus()

        return None

    #-----------------------------------------------------
    # CONTROL PANEL FUNCTIONS
    #-----------------------------------------------------

    def create_control_panel(self):
        project = self.current_project()

        # Create all widgets
        widgets = self.controls.build_control_panel(project)

        if widgets == {}:
            return None

        widgets["UC_panel"]["atoms_toggle"].stateChanged.connect(lambda state: project.set_visibility("atoms", bool(state)))
        widgets["UC_panel"]["choose_atoms_button"].clicked.connect(lambda _: self.open_atom_dialog())

        widgets["UC_panel"]["bonds_toggle"].stateChanged.connect(lambda state: project.set_visibility("bonds", bool(state)))
        widgets["UC_panel"]["choose_bonds_button"].clicked.connect(lambda: self.open_bond_dialog())

        # Connect widgets on the isosurface panel
        widgets["isosurface_panel"]["isosurface_toggle"].stateChanged.connect(lambda state: project.set_visibility("isosurface", bool(state)))
        widgets["isosurface_panel"]["isosurface_level_slider"].valueChanged.connect(lambda value: self.update_isosurface("level", value, widgets["isosurface_panel"]))
        widgets["isosurface_panel"]["color_button"].clicked.connect(lambda: self.update_isosurface("color", ControlPanel.choose_color(project.plotting_data["isosurface"]["color"]), widgets["isosurface_panel"]))
        widgets["isosurface_panel"]["opacity_slider"].valueChanged.connect(lambda value: self.update_isosurface("opacity", value, widgets["isosurface_panel"]))

        # Connect the widgets in the TuTraSt panel
        widgets["TuTraSt_panel"]["view_mode_dropdown"].currentTextChanged.connect(lambda text: project.set_view_mode("volume" if text == "Volumetric view" else "graph"))
        widgets["TuTraSt_panel"]["merge_tree_button"].clicked.connect(lambda _: self.plot_merge_tree())

        # Connect widgets in tunnel system subpanels
        for tunnel_system, tunnel_system_widgets in widgets["TuTraSt_panel"]["tunnel_systems"].items():
            tunnel_system_widgets["toggle"].stateChanged.connect(lambda state, ID = tunnel_system: project.set_group_visibility("tunnel_systems", ID, bool(state)))
            tunnel_system_widgets["choose_objects_button"].clicked.connect(lambda _, ID = tunnel_system: self.open_group_dialog(f"tunnel_systems", ID))
            tunnel_system_widgets["color_button"].clicked.connect(lambda _, ID = tunnel_system, ts_w = tunnel_system_widgets: self.update_group("tunnel_systems", ID, "color", ControlPanel.choose_color(project.plotting_data["tunnel_systems"][str(ID)]["color"]), ts_w))
            tunnel_system_widgets["opacity_slider"].valueChanged.connect(lambda value, ID = tunnel_system, ts_w = tunnel_system_widgets: self.update_group("tunnel_systems", ID, "opacity", value, ts_w))

            # MEPs section; the radio buttons use the "checked and f()" to only parse the checked button data
            tunnel_system_widgets["MEP_toggle"].stateChanged.connect(lambda state, ID = tunnel_system: project.set_group_visibility("MEPs", ID, bool(state)))
            for direction, direction_radio in tunnel_system_widgets["MEP_radio_button"].items():
                direction_radio.toggled.connect(lambda checked, d=direction, ID = tunnel_system, ts_w = tunnel_system_widgets: checked and self.update_group("tunnel_systems", ID, "MEP", d, ts_w))
            tunnel_system_widgets["plot_E_diagram_button"].clicked.connect(lambda _, ID = tunnel_system: self.plot_E_diagram(ID))

        # Connect widgets in isolated groups subpanels
        for isolated_group, isolated_group_widgets in widgets["TuTraSt_panel"]["isolated_groups"].items():
            isolated_group_widgets["toggle"].stateChanged.connect(lambda state, ID = isolated_group: project.set_group_visibility("isolated_groups", ID, bool(state)))
            isolated_group_widgets["choose_objects_button"].clicked.connect(lambda _, ID = isolated_group: self.open_group_dialog("isolated_groups", ID))
            isolated_group_widgets["color_button"].clicked.connect(lambda _, ID = isolated_group, iso_w = isolated_group_widgets: self.update_group("isolated_groups", ID, "color", ControlPanel.choose_color(project.plotting_data["isolated_groups"][str(ID)]["color"]), iso_w))
            isolated_group_widgets["opacity_slider"].valueChanged.connect(lambda value, ID = isolated_group, iso_w = isolated_group_widgets: self.update_group("isolated_groups", ID, "opacity", value, iso_w))


    def open_dialog(self, key, dialog_factory):
        dialog = self.open_dialogs.get(key)

        if dialog is not None:
            dialog.raise_()
            dialog.activateWindow()
            return None

        dialog = dialog_factory()
        self.open_dialogs[key] = dialog

        dialog.finished.connect(
            lambda: self.open_dialogs.pop(key, None)
        )

        dialog.show()
        dialog.raise_()
        dialog.activateWindow()

        return None
    
    def open_atom_dialog(self):
        project = self.current_project()
        self.open_dialog(
            "atoms",
            lambda: AtomsDialog(project),
        )
    
        return None
    
    def open_bond_dialog(self):
        project = self.current_project()
        self.open_dialog(
            "bonds",
            lambda: BondsDialog(project))

        return None
    
    def open_group_dialog(self, group_type, ID):
        project = self.current_project()
        key = f"group_{group_type}_{ID}"
    
        self.open_dialog(
            key,
            lambda: GroupDialog(project, group_type, ID),
        )

        return None

    def update_isosurface(self, attribute, value, iso_widgets):
        project = self.current_project()

        # Change the plotting data
        isosurface_plotting_data = project.plotting_data["isosurface"]
        isosurface_plotting_data[attribute] = value
        project.set_plotting_data("isosurface", isosurface_plotting_data)

        # Update the label
        iso_widgets["isosurface_level_label"].setText(f"Energy: {isosurface_plotting_data['level']*project.metadata['E_step']:.2f} / {project.metadata['N_levels']*project.metadata['E_step']:.2f} kJ/mol")

        # Update the color button
        ControlPanel.set_button_color(iso_widgets["color_button"], isosurface_plotting_data["color"])

        return None

    def update_group(self, group_type, ID, attribute, value, group_widgets):
        project = self.current_project()

        # Change the plotting data
        group_plotting_data = project.plotting_data[group_type][str(ID)]
        group_plotting_data[attribute] = value
        project.set_group_plotting_data(group_type, ID, group_plotting_data)

        # Update the color button
        ControlPanel.set_button_color(group_widgets["color_button"], group_plotting_data["color"])

        return None
    
    def plot_merge_tree(self):
        project = self.current_project()
        create_merge_trees(project)

        return None

    def plot_E_diagram(self, ID):
        project = self.current_project()
        
        direction = project.plotting_data["tunnel_systems"][ID]["MEP"]
        MEP = project.tunnel_systems[ID]["MEPs"][direction]

        plot_MEP_energy_diagram(MEP, direction, project.plotting_data, project.basin_data, project.TS_data, project.metadata["grid_shape"])

        return None
        
    #-----------------------------------------------------
    # DROPDOWN MENU FUNCTIONS
    #-----------------------------------------------------

    def create_menus(self):

        menubar = self.menuBar()
        menubar.setNativeMenuBar(False)

        menus = {}

        # --- File menu ---

        file_menu = menubar.addMenu("File")
        
        file_actions = {}

        load_structure_action = QAction("Load structure", self)
        save_screenshot_action = QAction("Save screenshot", self)
        export_data_action = QAction("Export data", self)
        export_data_action.setEnabled(False)

        file_actions["load_structure"] = load_structure_action
        file_actions["save_screenshot"] = save_screenshot_action
        file_actions["export_data"] = export_data_action

        file_menu.addAction(load_structure_action)
        file_menu.addAction(save_screenshot_action)
        file_menu.addAction(export_data_action)

        load_structure_action.triggered.connect(self.load_file)
        save_screenshot_action.triggered.connect(self.save_screenshot)

        menus["file"] = {"menu": file_menu, "actions": file_actions}

        # --- Unit cell menu ---

        unit_cell_menu = menubar.addMenu("Unit cell")

        unit_cell_menu.setEnabled(False)

        menus["unit_cell"] = {"menu": unit_cell_menu, "actions": {}}

        # --- Diffusion menu ---

        diffusion_menu = menubar.addMenu("Diffusion")

        diffusion_menu.setEnabled(False)

        menus["diffusion"] = {"menu": diffusion_menu, "actions": {}}

        # --- save menus ---

        self.menus = menus

        return None
    
    def update_menus(self):
        project = self.current_project()

        if project is None:
            self.menus["unit_cell"]["menu"].clear()
            self.menus["unit_cell"]["menu"].setEnabled(False)
            self.menus["diffusion"]["menu"].clear()
            self.menus["diffusion"]["menu"].setEnabled(False)
            
            return None
 
        # --- Unit cell menu ---
 
        self.menus["unit_cell"]["menu"].clear()
        self.menus["unit_cell"]["menu"].setEnabled(True)
 
        unit_cell_actions = {}
 
        show_uc_outline_action = QAction("Show UC outline", self)
        show_uc_outline_action.setCheckable(True)
        show_uc_outline_action.setChecked(project.visibility["UC_outline"])
        show_uc_outline_action.toggled.connect(lambda state: project.set_visibility("UC_outline", state))
 
        show_uc_parameters_action = QAction("Show UC parameters", self)
        show_uc_parameters_action.setCheckable(True)
        show_uc_parameters_action.setChecked(project.visibility["UC_parameters"])
        show_uc_parameters_action.toggled.connect(lambda state: project.set_visibility("UC_parameters", state))
 
        create_supercell_action = QAction("Create a supercell", self)
        create_supercell_action.triggered.connect(lambda _: self.open_supercell_dialog())
        
 
        self.menus["unit_cell"]["menu"].addAction(show_uc_outline_action)
        self.menus["unit_cell"]["menu"].addAction(show_uc_parameters_action)
        self.menus["unit_cell"]["menu"].addSeparator()
        self.menus["unit_cell"]["menu"].addAction(create_supercell_action)
 
        unit_cell_actions["show_uc_outline"] = show_uc_outline_action
        unit_cell_actions["show_uc_parameters"] = show_uc_parameters_action
        unit_cell_actions["create_supercell"] = create_supercell_action
 
        self.menus["unit_cell"]["actions"] = unit_cell_actions
 
        # --- Diffusion menu ---
 
        self.menus["diffusion"]["menu"].clear()
 
        if project.kMC_data is None:
            self.menus["diffusion"]["menu"].setEnabled(False)
        else:
            self.menus["diffusion"]["menu"].setEnabled(True)

            diffusion_actions = {}

            for temperature in project.kMC_data:
                
                action = QAction(f"T = {float(temperature):g}", self)
                action.triggered.connect(lambda _, temperature=temperature: self.open_diffusion_dialog(temperature))
                self.menus["diffusion"]["menu"].addAction(action)
                diffusion_actions[temperature] = action

            action = QAction(f"Arrhenius plot", self)
            action.triggered.connect(lambda _: self.open_arrhenius_dialog())
            if len(project.kMC_data) < 2:
                action.setEnabled(False)
            self.menus["diffusion"]["menu"].addAction(action)

            self.menus["diffusion"]["actions"] = diffusion_actions


        return None

    def open_diffusion_dialog(self, temperature):
        project = self.current_project()
    
        self.open_dialog(
            f"diffusion_{temperature}",
            lambda: DiffusionDialog(project, temperature, parent=self)
        )
    
        return None

    def open_arrhenius_dialog(self):
        project = self.current_project()

        self.open_dialog(
            f"Arrhenius Plot",
            lambda: ArrheniusDialog(project, parent=self)
        )
        return None

    def open_supercell_dialog(self):
        project = self.current_project()

        self.open_dialog(
            "supercell",
            lambda: SupercellDialog(project, parent=self)
        )

        return None
    

    