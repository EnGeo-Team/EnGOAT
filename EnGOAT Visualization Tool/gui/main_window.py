from PySide6.QtWidgets import (
    QMainWindow,
    QWidget,
    QHBoxLayout,
    QFileDialog,
    QTabWidget, QColorDialog
)
from PySide6.QtGui import QAction
from PySide6.QtCore import Qt

import os

from gui.control_panel import ControlPanel
from gui.pyvista_view import PyVistaView
from gui.project import Project
from data.parsers import get_them_files
from pathlib import Path

from gui.dialogs.tunnel_info_dialog import TunnelInfoDialog
from gui.dialogs.isolated_clusters_dialog import IsolatedClustersDialog
from gui.dialogs.basin_selection_dialog import BasinSelectionDialog
from gui.dialogs.merge_tree_selection_dialog import (MergeTreeSelectionDialog)

from gui.dialogs.transition_selection_dialog import TransitionSelectionDialog

from plotting.matplotlib.plot_MEP import plot_MEP_energy_diagram

class MainWindow(QMainWindow):

    def __init__(self):
        super().__init__()

        self.setWindowTitle("EnGOAT_Visualization")

        self.resize(1400, 900)

        self.create_menus()

        central = QWidget()
        self.setCentralWidget(central)

        layout = QHBoxLayout(central)

        self.controls = ControlPanel()
        self.tabs = QTabWidget()

        layout.addWidget(self.controls)
        layout.addWidget(self.tabs, stretch=1)

        #Connect toggles for UC
        self.controls.unit_cell_toggle.stateChanged.connect(self.toggle_unit_cell)
        self.controls.uc_params_toggle.stateChanged.connect(self.toggle_uc_params)

        self.controls.atoms_toggle.stateChanged.connect(self.toggle_atoms)
        self.controls.bonds_toggle.stateChanged.connect(self.toggle_bonds)

        self.controls.iso_toggle.stateChanged.connect(self.toggle_isosurface)
        self.controls.color_button.clicked.connect(self.pick_color)

        self.controls.plot_iso_button.clicked.connect(self.update_iso_params)

        self.controls.view_mode_combo.currentIndexChanged.connect(self.on_view_mode_changed)


        self.controls.labels_basins_cb.toggled.connect(
            lambda val: self._toggle_labels_basins(val)
        )

        self.controls.labels_transitions_cb.toggled.connect(
            lambda val: self._toggle_labels_TS(val)
        )

        self.controls.energies_basins_cb.toggled.connect(
            lambda val: self._toggle_energies_basins(val)
        )

        self.controls.energies_transitions_cb.toggled.connect(
            lambda val: self._toggle_energies_TS(val)
)

        #Load structure button connected action
        self.controls.load_button.clicked.connect(self.load_file)

        #Sync UI when changing tabs
        self.tabs.currentChanged.connect(self.on_tab_changed)

        #Closing tabs
        self.tabs.setTabsClosable(True)
        self.tabs.tabCloseRequested.connect(self.close_tab)

        #Open dialogue windows
        self.open_dialogs = []

        
    def connect_tunnel_widgets(self):

        for tunnel_name, widget_set in self.controls.tunnel_widgets.items():

            widget_set["show"].toggled.connect(
                lambda checked, tunnel=tunnel_name:
                    self.on_tunnel_visibility_changed(
                        tunnel,
                        checked
                    )
            )

            widget_set["color"].clicked.connect(
                 lambda checked=False, btn=widget_set["color"]:
                     self.choose_tunnel_color(btn)
             )
            
            widget_set["plot"].clicked.connect(
                lambda checked=False, tunnel=tunnel_name:
                    self.on_plot_tunnel_system(
                        tunnel
                    )
            )
            
            widget_set["mep_toggle"].toggled.connect(
                lambda checked, tunnel=tunnel_name:
                    self.on_mep_toggle_changed(
                        tunnel,
                        checked
                    )
            )
            for direction, rb in widget_set["mep_buttons"].items():

                rb.toggled.connect(
                    lambda checked,
                           tunnel=tunnel_name,
                           d=direction:
                    self.on_mep_direction_changed(
                        tunnel,
                        d,
                        checked
                    )
                )

            
            widget_set["plot_E"].clicked.connect(
                lambda checked=False, tunnel=tunnel_name:
                    self.on_plot_MEP_energy_diagram(tunnel)
            )


        if hasattr(self.controls, "isolated_cluster_widgets"):

            self.controls.isolated_cluster_widgets["show"].toggled.connect(
                self.on_isolated_cluster_visibility_changed
            )

            self.controls.isolated_cluster_widgets["color"].clicked.connect(
                    lambda checked=False,
                    btn=self.controls.isolated_cluster_widgets["color"]:
                        self.choose_tunnel_color(btn)
                )

            self.controls.isolated_cluster_widgets["plot"].clicked.connect(
                self.on_plot_isolated_clusters
            )


    #Tells you which project you are in 
    def current_project(self):
        return self.tabs.currentWidget()

    def on_tab_changed(self):
        project = self.current_project()
        if project is None:
            return
        
        self.controls.viewer = project.viewer   #Set the viewer for the control panel (upon changing tabs, change the viewer)
        
        #Called upon  tab change. For all toggles, freeze them, change them to the correct state, and then unfreeze them
        self.controls.unit_cell_toggle.blockSignals(True)
        self.controls.unit_cell_toggle.setChecked(project.show_unit_cell)
        self.controls.unit_cell_toggle.blockSignals(False)

        self.controls.uc_params_toggle.blockSignals(True)
        self.controls.uc_params_toggle.setChecked(project.show_uc_params)
        self.controls.uc_params_toggle.blockSignals(False)

        self.controls.atoms_toggle.blockSignals(True)
        self.controls.atoms_toggle.setChecked(project.show_atoms)
        self.controls.atoms_toggle.blockSignals(False)
        self.controls.choose_atoms_button.setEnabled(project.show_atoms)

        self.controls.bonds_toggle.blockSignals(True)
        self.controls.bonds_toggle.setChecked(project.show_bonds)
        self.controls.bonds_toggle.blockSignals(False)
        self.controls.choose_bonds_button.setEnabled(project.show_bonds)

        
        self.controls.iso_toggle.blockSignals(True)
        self.controls.iso_toggle.setChecked(project.show_isosurface)
        self.controls.iso_toggle.blockSignals(False)

        self.controls.iso_slider.setEnabled(project.show_isosurface)
        self.controls.opacity_input.setEnabled(project.show_isosurface)
        self.controls.color_button.setEnabled(project.show_isosurface)
        self.controls.plot_iso_button.setEnabled(project.show_isosurface)

        # ✅ restore isosurface UI
        self.controls.iso_slider.blockSignals(True)
        lv_Max = project.E_levels["N_steps"]
        E_step = project.E_levels["E_step"] 
        
        self.controls.iso_slider.blockSignals(True)
        self.controls.iso_slider.setValue(project.isosurface_params["level"])
        self.controls.iso_slider.blockSignals(False)

        self.controls.update_iso_range(lv_Max, E_step)

        self.controls.iso_slider.blockSignals(False)

        self.controls.opacity_input.blockSignals(True)
        self.controls.opacity_input.setText(str(project.isosurface_params["opacity"]))
        self.controls.opacity_input.blockSignals(False)

        self.controls.iso_color = project.isosurface_params["color"]

        r, g, b = project.isosurface_params["color"]
        self.controls.color_button.setStyleSheet(
            f"background-color: rgb({int(r*255)}, {int(g*255)}, {int(b*255)});"
        )

        #Restore graph/volumetric view

        self.controls.view_mode_combo.blockSignals(True)
        if project.view_mode == "graph":
            self.controls.view_mode_combo.setCurrentText("Graph view")
        else:
            self.controls.view_mode_combo.setCurrentText("Volumetric view")
        self.controls.view_mode_combo.blockSignals(False)

        #Tunnel system sections restore view
        self.controls.update_tunnel_systems(project.tunnel_system_layout, project.isolated_clusters, self.current_project())
        #Connect the tunnel system buttons
        self.connect_tunnel_widgets()

        #Update info menu
        self.update_info_menu()

    #Toggle functions for UC parameters
    def toggle_unit_cell(self, state):
        project = self.current_project()
        if project is None:
            return
        visible = bool(state)
        project.show_unit_cell = visible            #Set the state
        project.set_unit_cell_visibility(visible)   #Apply the state

    def toggle_uc_params(self, state):
        project = self.current_project()
        if project is None:
            return
        visible = bool(state)
        project.show_uc_params = visible
        project.set_unit_cell_params_visibility(visible)

    def toggle_atoms(self, state):
        project = self.current_project()
        if project is None:
            return

        visible = bool(state)
        project.set_atoms_visibility(visible)
        self.controls.choose_atoms_button.setEnabled(visible)   #Grey out the atoms button if toggle is switched off

    def toggle_bonds(self, state):
        project = self.current_project()
        if project is None:
            return

        visible = bool(state)
        project.set_bonds_visibility(visible)
        self.controls.choose_bonds_button.setEnabled(visible)

    def toggle_isosurface(self, state):
        project = self.current_project()
        if project is None:
            return
        
        visible = bool(state)
        project.set_isosurface_visibility(visible)

        self.controls.iso_slider.setEnabled(visible)
        self.controls.opacity_input.setEnabled(visible)
        self.controls.color_button.setEnabled(visible)
        self.controls.plot_iso_button.setEnabled(visible)


    def _toggle_labels_basins(self, val):
        project = self.current_project()
        project.show_labels_basins = val
        project.viewer.update_TuTraSt_plots(project)

    def _toggle_labels_TS(self, val):
        project = self.current_project()
        project.show_labels_TS = val
        project.viewer.update_TuTraSt_plots(project)

    def _toggle_energies_basins(self, val):
        project = self.current_project()
        project.show_energies_basins = val
        project.viewer.update_TuTraSt_plots(project)

    def _toggle_energies_TS(self, val):
        project = self.current_project()
        project.show_energies_TS = val
        project.viewer.update_TuTraSt_plots(project)


    def update_iso_params(self):
        project = self.current_project()
        if project is None:
            return

        project.update_isosurface_params(
            self.controls.iso_slider.value(),
            self.controls.get_opacity(),
            self.controls.iso_color
        )

        project.update_isosurface()

    def pick_color(self):
        project = self.current_project()
        if project is None:
            return
        color = QColorDialog.getColor()
        if color.isValid():
            self.controls.iso_color = (
                color.redF(),
                color.greenF(),
                color.blueF()
            )
            # ✅ show selected color
            self.controls.color_button.setStyleSheet(
                f"background-color: rgb({color.red()}, {color.green()}, {color.blue()});"
            )

    def on_view_mode_changed(self):

        text = self.controls.view_mode_combo.currentText()

        if text == "Graph view":
            mode = "graph"
        else:
            mode = "volume"

        project = self.current_project()
        project.set_view_mode(mode)
        project.viewer.update_basins(project)
        project.viewer.update_transitions(project)
        #project.viewer.update_basins(self.project)


    #Create an overhead menu
    def create_menus(self):
        menu_bar = self.menuBar()

        # --- FILE MENU ---

        file_menu = menu_bar.addMenu("File")

        load_action = QAction("Load", self)
        save_action = QAction("Save", self)

        load_action.triggered.connect(self.load_file)
        save_action.triggered.connect(self.save_file)

        file_menu.addAction(load_action)
        file_menu.addAction(save_action)

        # ======================
        # ✅ ADVANCED MENU
        # ======================
        self.advanced_menu = menu_bar.addMenu("Advanced")

        plot_basins_action = QAction("Plot individual basins", self)
        plot_transitions_action = QAction("Plot individual transitions", self)
        merge_tree_action = QAction("Merge tree diagram", self)

        # connect later
        plot_basins_action.triggered.connect(self.plot_individual_basins)
        plot_transitions_action.triggered.connect(self.plot_individual_transitions)
        merge_tree_action.triggered.connect(self.plot_merge_tree)

        self.advanced_menu.addAction(plot_basins_action)
        self.advanced_menu.addAction(plot_transitions_action)
        self.advanced_menu.addAction(merge_tree_action)

        # ======================
        # ✅ INFO MENU (dynamic)
        # ======================
        self.info_menu = menu_bar.addMenu("Info")




    def plot_individual_basins(self):

        project = self.current_project()
        if project is None:
            return

        dialog = BasinSelectionDialog(project)
        dialog.setAttribute(Qt.WA_DeleteOnClose)
        dialog.show()

        self.open_dialogs.append(dialog)

        dialog.destroyed.connect(
            lambda: self.open_dialogs.remove(dialog) if dialog in self.open_dialogs else None
        )

    
    def plot_individual_transitions(self):
        project = self.current_project()
        if not project:
            return

        dialog = TransitionSelectionDialog(project)
        dialog.setAttribute(Qt.WA_DeleteOnClose)
        dialog.show()

    def plot_merge_tree(self):
    
        project = self.current_project()
    
        if project is None:
            return
    
        dialog = MergeTreeSelectionDialog(
            project,
            self
        )
    
        dialog.setAttribute(Qt.WA_DeleteOnClose)
    
        dialog.show()
    
        self.open_dialogs.append(dialog)
    
        dialog.destroyed.connect(
            lambda:
            self.open_dialogs.remove(dialog)
            if dialog in self.open_dialogs
            else None
        )


    def update_info_menu(self):

        self.info_menu.clear()

        project = self.current_project()
        if project is None:
            return

        if not project.tunnel_system_layout:
            return

        for name in project.tunnel_system_layout.keys():

            action = QAction(name, self)

            # ✅ connect action
            action.triggered.connect(
                lambda checked, n=name: self.on_tunnel_info_selected(n)
            )

            self.info_menu.addAction(action)

        # ✅ Isolated clusters
        if project.isolated_clusters:

            self.info_menu.addSeparator()

            iso_action = QAction("Isolated clusters", self)
            iso_action.triggered.connect(self.on_isolated_clusters_selected)

            self.info_menu.addAction(iso_action)


    def on_tunnel_info_selected(self, name):

        project = self.current_project()
        if project is None:
            return

        # ✅ get tunnel info
        tunnel_info = project.tunnel_systems.get(name)
        if tunnel_info is None:
            return

        basin_ids = set(tunnel_info["basins"])
        ts_ids = set(tunnel_info["transitions"])

        # ✅ filter Basin_list
        basins = [b for b in project.Basin_list if b.ID in basin_ids]

        # ✅ filter TS_list
        transitions = [ts for ts in project.TS_list if ts.ID in ts_ids]

        # ✅ open dialog
        dialog = TunnelInfoDialog(name, basins, transitions, project.grid, project.tunnel_info)
        dialog.setAttribute(Qt.WA_DeleteOnClose)
        dialog.show()

        self.open_dialogs.append(dialog)
        dialog.destroyed.connect(lambda: self.open_dialogs.remove(dialog))

    def on_isolated_clusters_selected(self):

        project = self.current_project()
        if project is None:
            return

        clusters = project.isolated_clusters
        if not clusters:
            return

        dialog = IsolatedClustersDialog(project, clusters)
        dialog.setAttribute(Qt.WA_DeleteOnClose)
        dialog.show()

        self.open_dialogs.append(dialog)
        dialog.destroyed.connect(lambda: self.open_dialogs.remove(dialog))

    def on_tunnel_visibility_changed(self, tunnel_name, checked):

        project = self.current_project()

        project.set_tunnel_visibility(
            tunnel_name,
            checked
        )
    
    def on_isolated_cluster_visibility_changed(self, checked):

        project = self.current_project()

        project.set_isolated_cluster_visibility(
            checked
        )

    def on_plot_tunnel_system(
        self,
        tunnel_name
    ):

        project = self.current_project()

        widgets = self.controls.tunnel_widgets[tunnel_name]

        button = widgets["color"]

        color = getattr(
            button,
            "selected_color",
            (0.059, 0.322, 0.729)
        )

        try:
            opacity = float(
                widgets["opacity"].text()
            )
        except ValueError:
            opacity = 0.5

        project.update_tunnel_plotting(
            tunnel_name,
            color,
            opacity
        )

    def on_plot_isolated_clusters(self):

        project = self.current_project()

        widgets = self.controls.isolated_cluster_widgets

        color = getattr(
            widgets["color"],
            "selected_color",
            (0.839, 0.255, 0.243)
        )

        try:
            opacity = float(
                widgets["opacity"].text()
            )
        except ValueError:
            opacity = 0.5

        project.update_isolated_cluster_plotting(
            color,
            opacity
        )

    def choose_tunnel_color(self, button):

        color = QColorDialog.getColor()

        if not color.isValid():
            return

        button.setStyleSheet(
            f"""
            border: 1px solid #888;
            border-radius: 3px;
            background-color: rgb(
                {color.red()},
                {color.green()},
                {color.blue()}
            );
            """
        )

        # store temporarily on the button
        button.selected_color = (
            color.redF(),
            color.greenF(),
            color.blueF()
        )


    def on_mep_toggle_changed(
        self,
        tunnel_name,
        checked
    ):

        project = self.current_project()

        project.tunnel_systems_plotting[
            tunnel_name
        ]["show_MEP"] = checked

        project.refresh_MEP_styling()


    def on_mep_direction_changed(
        self,
        tunnel_name,
        direction,
        checked
    ):

        if not checked:
            return

        project = self.current_project()

        project.tunnel_systems_plotting[
            tunnel_name
        ]["MEP_direction"] = direction

        project.refresh_MEP_styling()

    #Create a new project, open a tab for it
    def create_new_project(self, name="Untitled"):

        project = Project()
        index = self.tabs.addTab(project, name)
        self.tabs.setCurrentIndex(index)

        self.controls.viewer = project.viewer

        return project
    
    #Simple tab closer
    def close_tab(self, index):
        widget = self.tabs.widget(index)
        widget.deleteLater()
        self.tabs.removeTab(index)

    #Loading function, starts a new project and stores information about it into the project
    def load_file(self):

        filename, _ = QFileDialog.getOpenFileName(
            self,
            "Open structure file",
            "",
            "Log files (*.log);;All files (*)"
        )

        if not filename:
            return

        print("Loading:", filename)
        directory = Path(filename).parent
        files = get_them_files(directory)
        
        #atoms = read_atom_info(directory)
        name = os.path.basename(files["cube_file"])
        project = self.create_new_project(name)

        project.load_structure(directory, files)

        lv_Max = project.E_levels["N_steps"]
        E_step = project.E_levels["E_step"]
        self.controls.iso_slider.setMaximum(lv_Max)
        self.controls.iso_label.setText(
            f"Energy: {self.controls.iso_slider.value()*E_step} / {lv_Max*E_step} kJ/mol")
        self.controls.update_iso_range(lv_Max, E_step)

        self.controls.update_tunnel_systems(project.tunnel_system_layout, project.isolated_clusters, self.current_project())
        self.update_info_menu()
        self.connect_tunnel_widgets()
        

    def on_plot_MEP_energy_diagram(self, tunnel_name):

        project = self.current_project()

        plot_MEP_energy_diagram(
            project,
            tunnel_name
        )

    def save_file(self):

        project = self.current_project()

        if project is None:
            return

        filename, _ = QFileDialog.getSaveFileName(
            self,
            "Save screenshot",
            "",
            "PNG Images (*.png)"
        )

        if not filename:
            return

        project.viewer.plotter.screenshot(filename)#, scale = 4)

        print(f"Saved screenshot to {filename}")






