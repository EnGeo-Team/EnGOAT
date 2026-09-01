from PySide6.QtWidgets import QWidget, QVBoxLayout
from PySide6.QtCore import Signal

import numpy as np
import matplotlib.pyplot as plt
import json

from gui.pyvista_view import PyVistaView
from plotting.helpers import element_color, covalent_radius, create_bond_data


class Project(QWidget):
    selected_basin_changed = Signal()
    atom_visibility_changed = Signal()

    def __init__(self):
        super().__init__()

        #Directory where the input files are stored
        self.location = None                                

        #All the information loaded from EnGOAT
        self.metadata = None
        self.unit_cell = None
        self.atoms = None
        self.bonds = None     #created on the fly
        self.basin_data = None
        self.TS_data = None
        self.tunnel_systems = None
        self.isolated_groups = None
        self.kMC_data = None

        self.level_matrix = None
        self.basin_matrix = None
        self.TS_matrix = None

        #Actor visibility states
        self.visibility = {
        
            # Unit cell
            "UC_outline": True,
            "UC_grid": False,
            "UC_parameters": False,

            # Atoms
            "atoms": True,
            "individual_atoms": None,
            "atom_labels": None,

            # Bonds
            "bonds": True,
            "individual_bonds": None,
            "bond_labels": None,

            # Isosurface
            "isosurface": False,

            # Basins
            "individual_basins": None,
            "basin_labels": None,
            "basin_energies": None,

            # Transition states
            "individual_TS": None,
            "TS_labels": None,
            "TS_energies": None,

            # Tunnel systems
            "tunnel_systems": None,

            # MEPs
            "MEPs": None,

            # Isolated groups
            "isolated_groups": None,
        }


        #Actor plotting states
        self.plotting_data = {
            "atoms": {},
            "isosurface": None,
            "basins": None,
            "TS": None,
            "tunnel_sysyems": None,
            "isolated_groups": None
        }

        self.invisible_bonds = []   #Some bonds are (partially) outside the supercell - they are not displayed! This list of bond_ID is filled when bond actors are created in PyVistaView

        #View mode
        self.view_mode = "volume"
        self.supercell = np.array([1, 1, 1])

        #Set up viewer
        layout = QVBoxLayout(self)
        self.viewer = PyVistaView(self)
        #self.viewer.set_project(self)
        layout.addWidget(self.viewer)

        #Additional states for the plotting dialog commands
        self.selected_basin = None
        self.select_radius = 4

        self.show_visible_bonds = False #To show bonds only between visible atoms

    #-----------------------------------------------------
    # lOADING THE STRUCTURE
    #-----------------------------------------------------

    def load_structure(self, location):

        self.location = location
        
        #Load EnGOAT_data
        filepath = location / "EnGOAT_data.json"
        with open(filepath, "r") as f:
            EnGOAT_data = json.load(f)
        
        self.metadata = EnGOAT_data["metadata"]
        self.metadata["grid_vectors"] = np.asarray(self.metadata["grid_vectors"], dtype=float)  #convert vectors to np arrays

        self.unit_cell = EnGOAT_data["unit_cell"]
        self.atoms = EnGOAT_data["atoms"]
        self.basin_data = EnGOAT_data["basin_data"]
        self.TS_data = EnGOAT_data["TS_data"]
        self.tunnel_systems = EnGOAT_data["tunnel_systems"]
        self.isolated_groups = EnGOAT_data["isolated_groups"]
        self.kMC_data = EnGOAT_data["kMC_data"]

        atom_types = sorted({atom["type"] for atom in self.atoms.values()})             #store colors and radii of atoms in dictionaries
        self.plotting_data["atoms"]["colors"] = {atom_type: element_color(atom_type) for atom_type in atom_types}
        self.plotting_data["atoms"]["radii"] = {atom_type: covalent_radius(atom_type) for atom_type in atom_types}

        self.bonds = create_bond_data(self.atoms, self.metadata["grid_shape"], self.metadata["grid_vectors"], self.plotting_data["atoms"]["radii"])

        #Load NumPy matrices
        self.level_matrix = np.load(location / "NumPy_matrices" / "Level_matrix.npy")
        self.basin_matrix = np.load(location / "NumPy_matrices" / "Basin_matrix.npy")
        self.TS_matrix = np.load(location / "NumPy_matrices" / "TS_matrix.npy")

        #Initialize actor visibilities and plotting states
        self.Initialize_actor_states()

        #Plot the initial state
        self.viewer.plotter.hide()
        self.viewer.create_all_actors()
        self.viewer.update_all_actors()
        self.viewer.plotter.show()

    def Initialize_actor_states(self):

        #Initialize missing Actor visibilities
        self.visibility["individual_atoms"] = {
            atom_ID: True
            for atom_ID in self.atoms
        }

        self.visibility["atom_labels"] = {
            atom_ID: False
            for atom_ID in self.atoms
        }

        self.visibility["individual_bonds"] = {
            bond_ID: True
            for bond_ID in self.bonds
        }

        self.visibility["bond_labels"] = {
            bond_ID: False
            for bond_ID in self.bonds
        }
        
        self.visibility["individual_basins"] = {
            basin_ID: False
            for basin_ID in self.basin_data
        }

        self.visibility["basin_labels"] = {
            basin_ID: False
            for basin_ID in self.basin_data
        }

        self.visibility["basin_energies"] = {
            basin_ID: False
            for basin_ID in self.basin_data
        }
        
        self.visibility["individual_TS"] = {
            TS_ID: False
            for TS_ID in self.TS_data
        }

        self.visibility["TS_labels"] = {
            TS_ID: False
            for TS_ID in self.TS_data
        }

        self.visibility["TS_energies"] = {
            TS_ID: False
            for TS_ID in self.TS_data
        }
        
        self.visibility["tunnel_systems"] = {
            tunnel_ID: False
            for tunnel_ID in self.tunnel_systems
        }
        
        self.visibility["MEPs"] = {
            tunnel_ID: False
            for tunnel_ID in self.tunnel_systems
        }
        
        self.visibility["isolated_groups"] = {
            group_ID: False
            for group_ID in self.isolated_groups
        }

        #Initialize actor plotting states

        self.plotting_data["isosurface"] = {"level": 1, "opacity": 0.5, "color": (0.529, 0.808, 0.922)}

        basin_plotting_data = {}
        TS_plotting_data = {}
        tunnel_system_plotting_data = {}
        isolated_groups_plotting_data = {}

        for tunnel_system_ID, data in self.tunnel_systems.items():
            first_direction = next(d for d, mep in data["MEPs"].items() if mep is not None)
            tunnel_system_plotting_data[tunnel_system_ID] = {"opacity": 0.5, "color": (0.059, 0.322, 0.729), "MEP": first_direction}
            for B_ID in data["basin_list"]:
                basin_plotting_data[str(B_ID)] = {"opacity": 0.5, "color": (0.059, 0.322, 0.729)}
            for TS_ID in data["TS_list"]:
                TS_plotting_data[str(TS_ID)] = {"opacity": 0.5, "color": (0.059, 0.322, 0.729)}

        for isolated_group_ID, data in self.isolated_groups.items():
            isolated_groups_plotting_data[isolated_group_ID] = {"opacity": 0.5, "color": (0.839, 0.255, 0.243)}
            for B_ID in data["basin_list"]:
                basin_plotting_data[str(B_ID)] = {"opacity": 0.5, "color": (0.839, 0.255, 0.243)}
            for TS_ID in data["TS_list"]:
                TS_plotting_data[str(TS_ID)] = {"opacity": 0.5, "color": (0.839, 0.255, 0.243)}

        self.plotting_data["basins"] = basin_plotting_data
        self.plotting_data["TS"] = TS_plotting_data
        self.plotting_data["tunnel_systems"] = tunnel_system_plotting_data
        self.plotting_data["isolated_groups"] = isolated_groups_plotting_data

        return None

    #-----------------------------------------------------
    # VISIBILITY CONTROL FUNCTIONS
    #-----------------------------------------------------

    def set_visibility(self, attribute, state):             #For attributes containing complex states (dictionary), the input state has to be complex too
        self.visibility[attribute] = state
        if attribute == "individual_atoms":
            self.update_bonds_from_atom_visibility()
            self.atom_visibility_changed.emit()

        self.viewer.update_all_actors()
        return None
    
    def set_group_visibility(self, group_type, ID, state):  #group type = 'tunnel_systems'/'isolated_grups'   /also MEPs?
        self.visibility[group_type][str(ID)] = state

        basin_list = []
        TS_list = []
        if group_type == "tunnel_systems":
            basin_list = self.tunnel_systems[str(ID)]["basin_list"]
            TS_list = self.tunnel_systems[str(ID)]["TS_list"]
        elif group_type == "isolated_groups":
            basin_list = self.isolated_groups[str(ID)]["basin_list"]
            TS_list = self.isolated_groups[str(ID)]["TS_list"]
        elif group_type == "MEPs":
            self.apply_MEP(ID, state)

        for B_ID in basin_list:
            self.visibility["individual_basins"][str(B_ID)] = state
        for TS_ID in TS_list:
            self.visibility["individual_TS"][str(TS_ID)] = state

        self.viewer.update_all_actors()
        return None

    #-----------------------------------------------------
    # PLOTTING STATE CONTROL FUNCTIONS
    #-----------------------------------------------------

    def set_view_mode(self, mode):
        self.view_mode = mode

        self.viewer.update_all_actors()
        return None
    
    def set_plotting_data(self, attribute, data):
        self.plotting_data[attribute] = data

        self.viewer.update_all_actors()
        return None

    def set_group_plotting_data(self, group_type, ID, data):
        self.plotting_data[group_type][str(ID)] = data

        color = data["color"]
        opacity = data["opacity"]

        if group_type == "tunnel_systems":
            basin_list = self.tunnel_systems[str(ID)]["basin_list"]
            TS_list = self.tunnel_systems[str(ID)]["TS_list"]
        elif group_type == "isolated_groups":
            basin_list = self.isolated_groups[str(ID)]["basin_list"]
            TS_list = self.isolated_groups[str(ID)]["TS_list"]

        for B_ID in basin_list:
            self.plotting_data["basins"][str(B_ID)] = {"color": color, "opacity": opacity}
        for TS_ID in TS_list:
            self.plotting_data["TS"][str(TS_ID)] = {"color": color, "opacity": opacity}

        if group_type == "tunnel_systems":
            if self.visibility["MEPs"][str(ID)]:
                self.apply_MEP(ID, True)

        self.viewer.update_all_actors()
        return None
    
    def apply_MEP(self, tunnel_ID, state):

        self.visibility["MEPs"][str(tunnel_ID)] = state

        direction = self.plotting_data["tunnel_systems"][str(tunnel_ID)]["MEP"]
        basin_list = [step["start_basin"] for step in self.tunnel_systems[str(tunnel_ID)]["MEPs"][direction]["path"]]
        TS_list = [step["transition_state"] for step in self.tunnel_systems[str(tunnel_ID)]["MEPs"][direction]["path"]]

        if state:
            cmap = plt.get_cmap("tab10")
            basin_plotting_data = {basin: {"color": cmap(i % 10)[:3], "opacity": 0.8} for i, basin in enumerate(basin_list)}
            TS_plotting_data = {TS: {"color": (0.0, 0.0, 0.0), "opacity": 0.8} for TS in TS_list}

        else:
            color = self.plotting_data["tunnel_systems"][str(tunnel_ID)]["color"]
            opacity = self.plotting_data["tunnel_systems"][str(tunnel_ID)]["opacity"]
            basin_plotting_data = {basin: {"color": color, "opacity": opacity} for basin in basin_list}
            TS_plotting_data = {TS: {"color": color, "opacity": opacity} for TS in TS_list}

        for B_ID , data in basin_plotting_data.items():
            self.plotting_data["basins"][str(B_ID)] = data
        
        for TS_ID , data in TS_plotting_data.items():
            self.plotting_data["TS"][str(TS_ID)] = data
        
        self.viewer.update_all_actors()
        return None

    def create_supercell(self, supercell):
        self.supercell = np.array(supercell)

        self.viewer.plotter.hide()
        self.viewer.clear()
        self.viewer.create_all_actors()
        self.viewer.update_all_actors()
        self.viewer.plotter.show()

        return None

    #For the clicking update
    def set_selected_basin(self, basin_id):
        self.selected_basin = basin_id
        self.selected_basin_changed.emit()


    # Trigger function to change bond visibility to show only bonds between visible atoms and change the state
    def set_show_visible_bonds(self, state):
        self.show_visible_bonds = bool(state)

        if self.show_visible_bonds:
            self.update_bonds_from_atom_visibility()

    # Change bond visibility to show only bonds between visible atoms
    def update_bonds_from_atom_visibility(self):
        if not self.show_visible_bonds:
            return

        atom_visibility = self.visibility["individual_atoms"]
        bond_visibility = self.visibility["individual_bonds"].copy()

        for bond_id, bond in self.bonds.items():
            atom1, atom2 = bond["atoms"]

            bond_visibility[bond_id] = (
                bool(atom_visibility[atom1])
                and bool(atom_visibility[atom2])
            )

        self.set_visibility(
            "individual_bonds",
            bond_visibility
        )