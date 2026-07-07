import numpy as np
from PySide6.QtWidgets import QWidget, QVBoxLayout
from gui.pyvista_view import PyVistaView
from pathlib import Path
from data.parsers import read_grid_info, read_atom_info, read_E_levels, tunnel_summary, read_tunnel_info
from plotting.Basins_transition_states import Basin, Transition_state, Organize_TuTraSt

class Project(QWidget):

    def __init__(self):
        super().__init__()

        #Information
        self.directory = None
        self.files = None
        
        #Layer visibility state
        self.show_unit_cell = True
        self.show_uc_params = True
        self.show_atoms = True
        self.show_bonds = True
        self.view_mode = "volume"

        self.show_isosurface = False


        self.show_labels_basins = False
        self.show_labels_TS = False

        self.show_energies_basins = False
        self.show_energies_TS = False


        self.isolated_clusters = None
        self.tunnel_systems = None
        self.tunnel_info = None

        self.E_levels = {"N_steps": 1,
                         "E_step": 1,
                         "E_cutoff": 1}

        self.isosurface_params = {
            "level": 1,
            "opacity": 0.5,
            "color": (0.529, 0.808, 0.922)
            }
        
        self.tunnel_system_layout = {}

        layout = QVBoxLayout(self)

        self.viewer = PyVistaView()
        self.viewer.set_project(self)
        layout.addWidget(self.viewer)

    def load_structure(self, directory, files):

        #Save file locations
        self.directory = directory
        self.files = files

        #Load matrices
        self.Level_matrix = np.load(f"{files["NumPy_matrices"]["Level_matrix"]}")
        self.Basin_matrix = np.load(f"{files["NumPy_matrices"]["Basin_matrix"]}")
        self.TS_matrix = np.load(f"{files["NumPy_matrices"]["TS_matrix"]}")
        self.Tunnel_matrix = np.load(f"{files["NumPy_matrices"]["Tunnel_matrix"]}")

        #Load UC info
        grid = read_grid_info(files["output"])
        self.grid = grid

        #Load atoms info
        atoms = read_atom_info(files["cube_file"])
        self.atoms = atoms

        #Load E_step info
        E_levels = read_E_levels(files["output"])
        self.E_levels = E_levels

        #Load Tunnel info
        self.tunnel_system_layout = tunnel_summary(files["Tunnels"])
        self.Basin_list, self.TS_list, self.tunnel_systems, self.isolated_clusters, self.tunnel_systems_plotting, self.isolated_clusters_plotting = Organize_TuTraSt(self.files)     

        self.basins_by_id = {
            basin.ID: basin
            for basin in self.Basin_list
        }
        
        self.TS_by_id = {
            ts.ID: ts
            for ts in self.TS_list
        } 

        self.tunnel_info = read_tunnel_info(files["output"])

        self.viewer.clear()

        #Initially open with plotting UC and UC_parameters
        self.viewer.draw_unit_cell(grid)
        self.viewer.display_UC_params(grid)

        self.viewer.plotter.hide()


        self.viewer.draw_atoms(atoms)

        self.viewer.draw_bonds(atoms)

        self.viewer.draw_isosurface(self.grid, self.Level_matrix, self.isosurface_params)

        #Apply visibility upon loading
        self.viewer.set_visibility("unit_cell_params", self.show_uc_params)
        self.viewer.set_visibility("unit_cell", self.show_unit_cell)

        self.viewer.set_visibility("atoms", self.show_atoms)
        self.viewer.set_visibility("bonds", self.show_bonds)

        self.viewer.set_visibility("isosurface", self.show_isosurface)
        self.viewer.update_TuTraSt_plots(self)
        self.viewer.initialize_labels(self)

        self.viewer.plotter.show()
        self.viewer.plotter.render()


    
    #Visibility control methods
    def set_unit_cell_params_visibility(self, visible):
        self.show_uc_params = visible
        self.viewer.set_visibility("unit_cell_params", visible)

    def set_unit_cell_visibility(self, visible):
        self.show_unit_cell = visible
        self.viewer.set_visibility("unit_cell", visible)

    def set_atoms_visibility(self, visible):
        self.show_atoms = visible
        if not hasattr(self.viewer, "atom_actors"):
            return
        if visible:
            self.viewer.update_visible_atoms(self.viewer.visible_elements)
        else:
            for actor in self.viewer.atom_actors.values():
                actor.SetVisibility(False)
                actor.pickable = False
            self.viewer.plotter.render()

    def set_bonds_visibility(self, visible):
        self.show_bonds = visible
        if not hasattr(self.viewer, "bond_actors"):
            return
        if visible:
            self.viewer.update_visible_bonds(self.viewer.visible_bonds)
        else:
            for actor in self.viewer.bond_actors.values():
                actor.SetVisibility(False)
            self.viewer.plotter.render()

    def set_isosurface_visibility(self, visible):
        self.show_isosurface = visible
        self.viewer.set_visibility("isosurface", visible)

    def update_isosurface_params(self, level, opacity, color):
        self.isosurface_params["level"] = level
        self.isosurface_params["opacity"] = opacity
        self.isosurface_params["color"] = color


    def set_tunnel_visibility(self, tunnel_name, visible):

        tunnel = self.tunnel_systems[tunnel_name]

        basin_ids = set(tunnel["basins"])
        ts_ids = set(tunnel["transitions"])

        for basin in self.Basin_list:
            if basin.ID in basin_ids:
                basin.visible = visible

        for ts in self.TS_list:
            if ts.ID in ts_ids:
                ts.visible = visible
        
        self.tunnel_systems_plotting[tunnel_name]["basins"]["visible"] = visible
        self.tunnel_systems_plotting[tunnel_name]["TS"]["visible"] = visible
        self.viewer.update_TuTraSt_plots(self)

    def set_isolated_cluster_visibility(self, visible):

        for cluster in self.isolated_clusters.values():

            basin_ids = set(cluster["basins"])
            ts_ids = set(cluster["transitions"])

            for basin in self.Basin_list:
                if basin.ID in basin_ids:
                    basin.visible = visible

            for ts in self.TS_list:
                if ts.ID in ts_ids:
                    ts.visible = visible

        self.isolated_clusters_plotting["basins"]["visible"] = visible
        self.isolated_clusters_plotting["TS"]["visible"] = visible
        self.viewer.update_TuTraSt_plots(self)

    
    def update_isosurface(self):
        self.viewer.draw_isosurface(self.grid, self.Level_matrix, self.isosurface_params)


    def set_view_mode(self, mode):
        self.view_mode = mode

    def update_tunnel_plotting(
        self,
        tunnel_name,
        color,
        opacity
    ):

        self.tunnel_systems_plotting[tunnel_name]["basins"]["color"] = color
        self.tunnel_systems_plotting[tunnel_name]["TS"]["color"] = color

        self.tunnel_systems_plotting[tunnel_name]["basins"]["opacity"] = opacity
        self.tunnel_systems_plotting[tunnel_name]["TS"]["opacity"] = opacity

        tunnel = self.tunnel_systems[tunnel_name]

        basin_ids = set(tunnel["basins"])
        ts_ids = set(tunnel["transitions"])

        for basin in self.Basin_list:

            if basin.ID in basin_ids:

                basin.color = color
                basin.opacity = opacity

        for ts in self.TS_list:

            if ts.ID in ts_ids:

                ts.color = color
                ts.opacity = opacity

        self.viewer.update_TuTraSt_plots(self)
        self.refresh_MEP_styling()

    def refresh_MEP_styling(self):

        # -------------------------
        # Restore default tunnel colors
        # -------------------------

        for tunnel_name, tunnel in self.tunnel_systems.items():

            plotting = self.tunnel_systems_plotting[tunnel_name]

            color = plotting["basins"]["color"]
            opacity = plotting["basins"]["opacity"]

            for basin_id in tunnel["basins"]:

                basin = self.basins_by_id[basin_id]

                basin.color = color
                basin.opacity = opacity

            for ts_id in tunnel["transitions"]:

                ts = self.TS_by_id[ts_id]

                ts.color = color
                ts.opacity = opacity

        # -------------------------
        # Apply active MEP overlays
        # -------------------------

        for tunnel_name, plotting in (
            self.tunnel_systems_plotting.items()
        ):

            if plotting["show_MEP"]:

                self.apply_MEP_styling(
                    tunnel_name
                )

        self.viewer.update_TuTraSt_plots(self)

    def update_isolated_cluster_plotting(
        self,
        color,
        opacity
    ):

        self.isolated_clusters_plotting["basins"]["color"] = color
        self.isolated_clusters_plotting["basins"]["opacity"] = opacity

        self.isolated_clusters_plotting["TS"]["color"] = color
        self.isolated_clusters_plotting["TS"]["opacity"] = opacity

        for cluster in self.isolated_clusters.values():

            for basin_id in cluster["basins"]:

                self.basins_by_id[basin_id].color = color
                self.basins_by_id[basin_id].opacity = opacity

            for ts_id in cluster["transitions"]:

                self.TS_by_id[ts_id].color = color
                self.TS_by_id[ts_id].opacity = opacity

        self.viewer.update_TuTraSt_plots(self)

    def apply_MEP_styling(self, tunnel_name):

        plotting = self.tunnel_systems_plotting[tunnel_name]

        if not plotting["show_MEP"]:
            return

        direction = plotting["MEP_direction"]

        mep = self.tunnel_systems[tunnel_name]["MEPs"][direction]

        if mep is None:
            return

        # Basin colors
        for basin_id, color in zip(
            mep["basin_ids"],
            mep["basin_colors"]
        ):

            basin = self.basins_by_id.get(basin_id)

            if basin is not None:

                basin.color = color
                basin.opacity = mep["opacity"]

        # TS colors
        for ts_id in mep["ts_ids"]:

            ts = self.TS_by_id.get(ts_id)

            if ts is not None:

                ts.color = mep["TS_color"]
                ts.opacity = mep["opacity"]