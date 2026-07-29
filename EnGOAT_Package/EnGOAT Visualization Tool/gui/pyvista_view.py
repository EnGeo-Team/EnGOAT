from PySide6.QtWidgets import QWidget, QVBoxLayout, QTextEdit, QHBoxLayout, QPushButton, QButtonGroup
from PySide6.QtCore import Qt
from pyvistaqt import QtInteractor
from PySide6.QtGui import QTextCursor

import numpy as np
import pyvista as pv
from plotting.helpers import pbc_images
from plotting.helpers import get_pbc_ts_geometry, voxel_surface
from plotting.helpers import all_bonds_inside_UC
from plotting.helpers import cartesian_to_fractional

from datetime import datetime

#Clicking
from enum import Enum
from dataclasses import dataclass

class ClickMode(Enum):
    SELECT = 0
    DISTANCE = 1
    ANGLE = 2

@dataclass
class PickedObject:
    kind: str
    ID: str
    center: np.ndarray
    copies: list[np.ndarray]


class PyVistaView(QWidget):

    def __init__(self, project):
        super().__init__()

        #Set the correspomding project that contains data and states
        self.project = project 

        #Plotting widget
        layout = QVBoxLayout(self)
        self.plotter = QtInteractor(self)

        self.plotter.setFocusPolicy(Qt.StrongFocus) #?
        self.plotter.setFocus() #?

        self.plotter.enable_depth_peeling() #?
        
        self.plotter.show_axes()
        layout.addWidget(self.plotter, stretch=1)

        self.plotter.enable_mesh_picking(
            callback=self.mouse_clicked,
            use_actor=True,
            show=False,
            show_message=False
        )

        #Info box widget and clicking buttons
        self.info_box = QTextEdit()
        self.info_box.setReadOnly(True)
        self.info_box.setPlaceholderText(f"EnGOAT Visualization Tool")#ADD LATER
        self.info_box.clear()
        self.info_box.setStyleSheet("""
            QTextEdit {
                background-color: #f8f8f8;
                border: 1px solid #aaa;
                padding: 8px;
                font-family: Consolas;
                font-size: 12px;
            }
        """)
        self.info_box.setFixedHeight(100) 

        self.select_button = QPushButton("•")
        self.distance_button = QPushButton("↔")
        self.angle_button = QPushButton("∠")
        #self.angle_button = QPushButton("∡")

        for button in (self.select_button,
               self.distance_button,
               self.angle_button):
            button.setCheckable(True)
            button.setStyleSheet("""
                QPushButton {
                    font-size: 20px;
                    font-weight: bold;
                }
            """)

        self.tool_group = QButtonGroup(self)
        self.tool_group.setExclusive(True)
        for button in (
            self.select_button,
            self.distance_button,
            self.angle_button,
        ):
            self.tool_group.addButton(button)

        self.select_button.setChecked(True)

        bottom = QWidget()
        bottom_layout = QHBoxLayout(bottom)
        bottom_layout.setContentsMargins(0, 0, 0, 0)

        bottom_layout.addWidget(self.info_box, stretch=1)

        # Button container
        button_widget = QWidget()
        button_widget.setFixedHeight(self.info_box.height())

        buttons = QVBoxLayout(button_widget)
        buttons.setContentsMargins(0, 0, 0, 0)
        buttons.setSpacing(5)  # small gap between buttons

        button_size = 30

        for button in (
            self.select_button,
            self.distance_button,
            self.angle_button,
        ):
            button.setFixedSize(button_size, button_size)

        buttons.addWidget(self.select_button)
        buttons.addWidget(self.distance_button)
        buttons.addWidget(self.angle_button)

        bottom_layout.addWidget(button_widget)

        layout.addWidget(bottom)

        self.select_button.clicked.connect(
            lambda: self.set_click_mode(ClickMode.SELECT)
        )

        self.distance_button.clicked.connect(
            lambda: self.set_click_mode(ClickMode.DISTANCE)
        )

        self.angle_button.clicked.connect(
            lambda: self.set_click_mode(ClickMode.ANGLE)
        )

        #Actor storage

        self.actors = {

            "UC": {"outline": [], "UC_grid": [], "parameters": []},

            "atoms": {"volume": {}, "outline": {}, "labels": {}},

            "bonds": {"volume": {}, "labels": {}},

            "isosurface": {},

            "basins": {"volume": {}, "graph": {}, "labels": {}, "energies": {}},

            "TS": {"volume": {}, "graph": {}, "labels": {}, "energies": {}},
        }

        #PBC atom and basin positions (for clicking mostly)
        self.PBC_atom_pos = {}
        self.PBC_basin_pos = {}

        # Clicking
        self.actor_metadata = {}    # store ALL clickable actors as keys and {"kind": "atom", "ID": "Li1"} as values
        self.click_mode = ClickMode.SELECT
        self.selected_objects = []
        self.click_actors = []

    #-----------------------------------------------------
    # MASTER FUNCTIONS
    #-----------------------------------------------------

    # Reset button function
    def clear(self):

        def remove_actors(obj):

            if isinstance(obj, dict):
                for value in obj.values():
                    remove_actors(value)

            elif isinstance(obj, (list, tuple)):
                for value in obj:
                    remove_actors(value)

            elif obj is not None:
                self.plotter.remove_actor(obj)

        remove_actors(self.actors)
        self.clear_click_state()

        self.plotter.enable_depth_peeling()
        self.actors = {

            "UC": {"outline": [], "UC_grid": [], "parameters": []},

            "atoms": {"volume": {}, "outline": {}, "labels": {}},

            "bonds": {"volume": {}, "labels": {}},

            "isosurface": {},

            "basins": {"volume": {}, "graph": {}, "labels": {}, "energies": {}},

            "TS": {"volume": {}, "graph": {}, "labels": {}, "energies": {}},
        }
        self.plotter.show_axes()
        self.plotter.reset_camera()

    def create_all_actors(self):
        
        self.create_cartesian_grid()
        
        self.create_UC_actors()
        print("0    %")
        self.create_atom_actors()
        print("20    %")
        self.create_bond_actors()
        print("40    %")
        self.create_isosurface_actors()
        print("60    %")
        self.create_basin_actors()
        print("80    %")
        self.create_TS_actors()
        print("100    %")

        return None

    def update_all_actors(self):

        self.update_UC_actors()
        self.update_atom_actors()
        self.update_bond_actors()
        self.update_isosurface_actors()
        self.update_basin_actors()
        self.update_TS_actors()

        self.plotter.render()

        return None

    def create_cartesian_grid(self):
    
        Na, Nb, Nc = self.project.metadata["grid_shape"]
        Sa, Sb, Sc = self.project.supercell
    
        a_vec = np.asarray(self.project.metadata["grid_vectors"][0], dtype=float)
        b_vec = np.asarray(self.project.metadata["grid_vectors"][1], dtype=float)
        c_vec = np.asarray(self.project.metadata["grid_vectors"][2], dtype=float)
    
        # Indices over the entire supercell
        i = np.arange(Na * Sa)
        j = np.arange(Nb * Sb)
        k = np.arange(Nc * Sc)
    
        I, J, K = np.meshgrid(i, j, k, indexing="ij")
    
        self.cartesian_points = (
            I[..., None] * a_vec +
            J[..., None] * b_vec +
            K[..., None] * c_vec
        )
    
        return None

    #-----------------------------------------------------
    # UC ACTORS
    #-----------------------------------------------------

    def create_UC_actors(self):

        Na, Nb, Nc = self.project.metadata["grid_shape"]

        a_vec = self.project.metadata["grid_vectors"][0]
        b_vec = self.project.metadata["grid_vectors"][1]
        c_vec = self.project.metadata["grid_vectors"][2]

        sx, sy, sz = map(int, self.project.supercell)

        A = Na * a_vec * sx
        B = Nb * b_vec * sy
        C = Nc * c_vec * sz

        def angle(v1, v2):
            return np.degrees(
                np.arccos(
                    np.clip(
                        np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)),
                        -1.0,
                        1.0
                    )
                )
            )

        alpha = angle(b_vec, c_vec)
        beta  = angle(a_vec, c_vec)
        gamma = angle(a_vec, b_vec)

        corners = np.array([
            [0, 0, 0],
            A,
            B,
            C,
            A + B,
            A + C,
            B + C,
            A + B + C
        ])

        edges = [
            (0,1),(0,2),(0,3),
            (1,4),(1,5),
            (2,4),(2,6),
            (3,5),(3,6),
            (4,7),(5,7),(6,7)
        ]

        # Main supercell outline
        for i, j in edges:
            actor = self.plotter.add_lines(
                np.vstack((corners[i], corners[j])),
                color="black",
                width=3
            )
            actor.SetPickable(False)
            self.actors["UC"]["outline"].append(actor)

        # Internal unit-cell grid (no outer outline)
        if (sx, sy, sz) != (1, 1, 1):
        
            ua = Na * a_vec
            ub = Nb * b_vec
            uc = Nc * c_vec
        
            def add_grid_line(p0, p1):
                actor = self.plotter.add_lines(
                    np.vstack((p0, p1)),
                    color="grey",
                    width=1
                )
                actor.GetProperty().SetLineStipplePattern(0xF0F0)
                actor.GetProperty().SetLineStippleRepeatFactor(2)
                actor.SetPickable(False)
                self.actors["UC"]["UC_grid"].append(actor)
        
            # Internal planes perpendicular to a
            for i in range(1, sx):
                for j in range(sy + 1):
                    for k in range(sz + 1):
                        p0 = i * ua + j * ub + k * uc
                        p1 = i * ua + (j + 1) * ub + k * uc if j < sy else None
        
                        if p1 is not None:
                            add_grid_line(p0, p1)
        
                        p0 = i * ua + j * ub + k * uc
                        p1 = i * ua + j * ub + (k + 1) * uc if k < sz else None
        
                        if p1 is not None:
                            add_grid_line(p0, p1)
        
        
            # Internal planes perpendicular to b
            for j in range(1, sy):
                for i in range(sx + 1):
                    for k in range(sz + 1):
                    
                        p0 = i * ua + j * ub + k * uc
        
                        if i < sx:
                            add_grid_line(
                                p0,
                                (i + 1) * ua + j * ub + k * uc
                            )
        
                        if k < sz:
                            add_grid_line(
                                p0,
                                i * ua + j * ub + (k + 1) * uc
                            )
        
        
            # Internal planes perpendicular to c
            for k in range(1, sz):
                for i in range(sx + 1):
                    for j in range(sy + 1):
                    
                        p0 = i * ua + j * ub + k * uc
        
                        if i < sx:
                            add_grid_line(
                                p0,
                                (i + 1) * ua + j * ub + k * uc
                            )
        
                        if j < sy:
                            add_grid_line(
                                p0,
                                i * ua + (j + 1) * ub + k * uc
                            )
        
        # Labels
        labels = [
            (A/2, f"a = {np.linalg.norm(A):.2f}Å"),
            (B/2, f"b = {np.linalg.norm(B):.2f}Å"),
            (C/2, f"c = {np.linalg.norm(C):.2f}Å"),
            (0.2*B + 0.2*C, rf"$\alpha$ = {alpha:.1f}°"),
            (0.2*A + 0.2*C, rf"$\beta$ = {beta:.1f}°"),
            (0.2*A + 0.2*B, rf"$\gamma$ = {gamma:.1f}°"),
        ]

        for pos, text in labels:
            actor = self.plotter.add_point_labels(
                [pos],
                [text],
                font_size=14,
                text_color="black",
                shape=None,
                always_visible=True,
                show_points=False
            )
            actor.SetPickable(False)
            self.actors["UC"]["parameters"].append(actor)

        return None
    
    def update_UC_actors(self):

        for actor in self.actors["UC"]["outline"]:
            actor.SetVisibility(self.project.visibility["UC_outline"])

        for actor in self.actors["UC"]["parameters"]:
            actor.SetVisibility(self.project.visibility["UC_parameters"])

        for actor in self.actors["UC"]["UC_grid"]:
            actor.SetVisibility(self.project.visibility["UC_grid"])

        return None
    
    #-----------------------------------------------------
    # ATOM ACTORS
    #-----------------------------------------------------
    
    def create_atom_actors(self):

        for atom_ID, data in self.project.atoms.items():
            color = self.project.plotting_data["atoms"]["colors"][data["type"]]
            radius = self.project.plotting_data["atoms"]["radii"][data["type"]]

            positions = pbc_images(data["center"], self.project.metadata["grid_shape"], self.project.metadata["grid_vectors"], self.project.supercell)

            self.PBC_atom_pos[atom_ID] = positions

            poly = pv.PolyData(positions)
            poly["radius"] = np.full(len(positions), radius)
            sphere = pv.Sphere(theta_resolution=24, phi_resolution=24)
            glyphs = poly.glyph(scale="radius", geom=sphere, orient=False)

            actor = self.plotter.add_mesh(
                glyphs,
                color = color,
                lighting=True,
                ambient=0.4,
                diffuse=0.8,
                specular=0.0,
                opacity=1.0
            )
            actor.SetPickable(True)
            self.actor_metadata[actor] = {"kind": "atoms", "ID": atom_ID}

            self.actors["atoms"]["volume"][atom_ID] = actor

            outline = self.plotter.add_silhouette(glyphs, color = "black", line_width = 2)
            outline.SetPickable(False)

            self.actors["atoms"]["outline"][atom_ID] = outline

            label = self.plotter.add_point_labels(
                positions,
                [atom_ID] * len(positions),   # repeat the atom ID for every image
                font_size=12,
                text_color="black",
                shape=None,
                show_points=False,
                always_visible=True
            )

            label.SetPickable(False)
            self.actors["atoms"]["labels"][atom_ID] = label
        
        return None
    
    def update_atom_actors(self):

        atoms_visibility = self.project.visibility["atoms"]

        for atom_ID in self.project.atoms:

            label_visibility = self.project.visibility["atom_labels"][atom_ID]

            individual_atom_visibility = self.project.visibility["individual_atoms"][atom_ID]
            self.actors["atoms"]["volume"][atom_ID].SetVisibility(individual_atom_visibility and atoms_visibility)
            self.actors["atoms"]["volume"][atom_ID].SetPickable(individual_atom_visibility and atoms_visibility)
            self.actors["atoms"]["outline"][atom_ID].SetVisibility(individual_atom_visibility and atoms_visibility)
            self.actors["atoms"]["labels"][atom_ID].SetVisibility(individual_atom_visibility and atoms_visibility and label_visibility)

        return None
    
    #-----------------------------------------------------
    # BOND ACTORS
    #-----------------------------------------------------
    
    def create_bond_actors(self):

        invisible_bonds = []

        # Helper function to make half of a bond
        def make_half_bond(p_start, p_end, color):
        
            cyl = pv.Cylinder(
                center=0.5 * (p_start + p_end),
                direction=p_end - p_start,
                radius=0.08,
                height=np.linalg.norm(p_end - p_start)
            )

            cyl.cell_data["colors"] = np.tile(
                color,
                (cyl.n_cells, 1)
            )

            return cyl


        for bond_ID, data in self.project.bonds.items():

            copies = all_bonds_inside_UC(data, self.project.metadata, self.project.atoms, self.project.supercell)

            if copies:

                atom1, atom2 = data["atoms"]
                color1 = self.project.plotting_data["atoms"]["colors"][self.project.atoms[atom1]["type"]]
                color2 = self.project.plotting_data["atoms"]["colors"][self.project.atoms[atom2]["type"]]

                bond_pieces = []
                centers = []

                for copy in copies:

                    p1, p2 = copy

                    mid = 0.5 * (p1 + p2)

                    bond_pieces.append(
                        make_half_bond(p1, mid, color1)
                    )

                    bond_pieces.append(
                        make_half_bond(mid, p2, color2)
                    )

                    centers.append(mid)

                bond_mesh = bond_pieces[0]

                for piece in bond_pieces[1:]:
                    bond_mesh = bond_mesh + piece

                
                bond_actor = self.plotter.add_mesh(
                    bond_mesh,
                    scalars="colors",
                    rgb=True,
                    preference="cell",
                    lighting=True,
                    ambient=0.3,
                    diffuse=0.8,
                    specular=0.1,
                    smooth_shading=True
                )
                bond_actor.SetPickable(False)

                self.actors["bonds"]["volume"][bond_ID] = bond_actor

                label = self.plotter.add_point_labels(
                    centers,
                    [bond_ID]*len(centers),
                    font_size=12,
                    text_color="black",
                    shape=None,
                    show_points=False,
                    always_visible=True
                )

                label.SetPickable(False)
                self.actors["bonds"]["labels"][bond_ID] = label
            
            else:   #For the cases where bond is not wholy in the unit/supercell
                self.actors["bonds"]["volume"][bond_ID] = None
                self.actors["bonds"]["labels"][bond_ID] = None
                invisible_bonds.append(bond_ID)

        for bond_ID in invisible_bonds:
            self.project.visibility["individual_bonds"][bond_ID] = True
        
        self.project.invisible_bonds = invisible_bonds
    
        return None

    def update_bond_actors(self):

        bond_visibility = self.project.visibility["bonds"]

        for bond_ID in self.project.bonds:
            label_visibility = self.project.visibility["bond_labels"][bond_ID]

            if self.actors["bonds"]["volume"][bond_ID] is not None: # Check if the bond actor exists - if the bond is within the unit/supercell
                individual_bond_visibility = self.project.visibility["individual_bonds"][bond_ID]
                self.actors["bonds"]["volume"][bond_ID].SetVisibility(individual_bond_visibility and bond_visibility)
                self.actors["bonds"]["labels"][bond_ID].SetVisibility(individual_bond_visibility and bond_visibility and label_visibility)

        return None

    #-----------------------------------------------------
    # ISOSURFACE ACTOR
    #-----------------------------------------------------

    def create_isosurface_actors(self):

        N_levels = self.project.metadata["N_levels"]
        color = self.project.plotting_data["isosurface"]["color"]
        opacity = self.project.plotting_data["isosurface"]["opacity"]
        
        a_vec = self.project.metadata["grid_vectors"][0]
        b_vec = self.project.metadata["grid_vectors"][1]
        c_vec = self.project.metadata["grid_vectors"][2]

        for level in range(1, N_levels+1):

            mask = np.tile(self.project.level_matrix <= level, self.project.supercell)

            surface = voxel_surface(mask,  self.cartesian_points, a_vec, b_vec, c_vec)

            actor = self.plotter.add_mesh(
                surface,
                color=color,
                opacity=opacity,
                lighting = True,
                smooth_shading=True,
                ambient=0.3,
                diffuse=0.8,
                specular=0.0,
                reset_camera=False,
                show_edges=False
            )
            actor.SetPickable(False)

            self.actors["isosurface"][str(level)] = actor

        return None
    
    def update_isosurface_actors(self):

        level = self.project.plotting_data["isosurface"]["level"]
        color = self.project.plotting_data["isosurface"]["color"]
        opacity = self.project.plotting_data["isosurface"]["opacity"]
        visible = self.project.visibility["isosurface"]

        for iso_level, isosurface_actor in self.actors["isosurface"].items():
            
            isosurface_actor.GetProperty().SetColor(*color)
            isosurface_actor.GetProperty().SetOpacity(opacity)
            isosurface_actor.SetVisibility(visible and level == int(iso_level))

        return None

    #-----------------------------------------------------
    # BASIN ACTORS
    #-----------------------------------------------------

    def create_basin_actors(self):

        a_vec = self.project.metadata["grid_vectors"][0]
        b_vec = self.project.metadata["grid_vectors"][1]
        c_vec = self.project.metadata["grid_vectors"][2]

        for B_ID, data in self.project.plotting_data["basins"].items():

            color = data["color"]
            opacity = data["opacity"]

            #Volume actor
            mask = np.tile(self.project.basin_matrix == int(B_ID), self.project.supercell)

            surface = voxel_surface(mask,  self.cartesian_points, a_vec, b_vec, c_vec)

            volume_actor = self.plotter.add_mesh(
                surface,
                color=color,
                opacity=opacity,
                lighting = True,
                smooth_shading=True,
                ambient=0.3,
                diffuse=0.8,
                specular=0.0,
                reset_camera=False,
                show_edges=False
            )
            volume_actor.SetPickable(False)
            self.actors["basins"]["volume"][B_ID] = volume_actor

            #Graph actor
            i, j, k = self.project.basin_data[B_ID]["center"]
            positions = pbc_images(self.cartesian_points[i, j, k], self.project.metadata["grid_shape"], self.project.metadata["grid_vectors"], self.project.supercell)

            self.PBC_basin_pos[B_ID] = positions

            poly = pv.PolyData(positions)
            graph_actor = self.plotter.add_mesh(
                poly,
                color=color,
                opacity=opacity,
                lighting = True,
                smooth_shading=True,
                ambient=0.3,
                diffuse=0.8,
                specular=0.0,
                point_size=20,
                render_points_as_spheres=True
            )
            graph_actor.SetPickable(True)
            self.actor_metadata[graph_actor] = {"kind": "basins", "ID": B_ID}
            self.actors["basins"]["graph"][B_ID] = graph_actor

            #Label actors
            labels = [f"B{B_ID}"] * len(positions)
            label_actor = self.plotter.add_point_labels(
                positions + np.array([0, 0, 0.15]),
                labels,
                font_size=10,
                show_points=False,
                shape=None,
                always_visible=True
            )
            label_actor.SetPickable(False)
            self.actors["basins"]["labels"][B_ID] = label_actor

            energies = [f"E={self.project.basin_data[B_ID]["E_min"]:.2f}kJ/mol"] * len(positions)
            energy_actor = self.plotter.add_point_labels(
                positions + np.array([0, 0, -0.15]),
                energies,
                font_size=9,
                text_color = "darkgray",
                show_points=False,
                shape=None,
                always_visible=True
            )
            energy_actor.SetPickable(False)
            self.actors["basins"]["energies"][B_ID] = energy_actor

        return None
    
    def update_basin_actors(self):

        view_mode = self.project.view_mode

        for B_ID in self.project.basin_data:

            color = self.project.plotting_data["basins"][B_ID]["color"]
            opacity = self.project.plotting_data["basins"][B_ID]["opacity"]
            visibility = self.project.visibility["individual_basins"][B_ID]

            label_visibility = self.project.visibility["basin_labels"][B_ID]
            energy_visibility = self.project.visibility["basin_energies"][B_ID]

            self.actors["basins"]["volume"][B_ID].GetProperty().SetColor(*color)
            self.actors["basins"]["volume"][B_ID].GetProperty().SetOpacity(opacity)
            self.actors["basins"]["volume"][B_ID].SetVisibility(visibility and (view_mode == "volume"))
                
            self.actors["basins"]["graph"][B_ID].GetProperty().SetColor(*color)
            self.actors["basins"]["graph"][B_ID].GetProperty().SetOpacity(opacity)
            self.actors["basins"]["graph"][B_ID].SetVisibility(visibility and (view_mode == "graph"))
            self.actors["basins"]["graph"][B_ID].SetPickable(visibility and (view_mode == "graph"))

            self.actors["basins"]["labels"][B_ID].SetVisibility(visibility and label_visibility)

            self.actors["basins"]["energies"][B_ID].SetVisibility(visibility and energy_visibility)

        return None
    
    #-----------------------------------------------------
    # TS ACTORS
    #-----------------------------------------------------

    def create_TS_actors(self):

        a_vec = np.array(self.project.metadata["grid_vectors"][0])
        b_vec = np.array(self.project.metadata["grid_vectors"][1])
        c_vec = np.array(self.project.metadata["grid_vectors"][2])

        for TS_ID, data in self.project.plotting_data["TS"].items():

            color = data["color"]
            opacity = data["opacity"]

            #Volume actor
            mask = np.tile(self.project.TS_matrix == int(TS_ID), self.project.supercell)

            surface = voxel_surface(mask,  self.cartesian_points, a_vec, b_vec, c_vec)

            volume_actor = self.plotter.add_mesh(
                surface,
                color=color,
                opacity=opacity,
                lighting = True,
                smooth_shading=True,
                ambient=0.3,
                diffuse=0.8,
                specular=0.0,
                reset_camera=False,
                show_edges=False
            )
            volume_actor.SetPickable(False)
            self.actors["TS"]["volume"][TS_ID] = volume_actor

            #Graph actor
            points, lines, centers = get_pbc_ts_geometry(self.project.TS_data[TS_ID], self.project.metadata["grid_shape"], (a_vec, b_vec, c_vec), self.project.supercell)

            poly = pv.PolyData()
            poly.points = np.array(points)
            poly.lines = np.array(lines)
            graph_actor = self.plotter.add_mesh(
                poly,
                color=color,
                lighting = True,
                smooth_shading=True,
                ambient=0.3,
                diffuse=0.8,
                specular=0.0,
                opacity=opacity,
                line_width=3
            )
            graph_actor.SetPickable(False)
            self.actors["TS"]["graph"][TS_ID] = graph_actor

            #Label actors
            labels = [f"TS{TS_ID}"] * len(centers)
            label_actor = self.plotter.add_point_labels(
                centers + np.array([0, 0, 0.15]),
                labels,
                font_size=10,
                show_points=False,
                shape=None, 
                always_visible=True
            )
            label_actor.SetPickable(False)
            self.actors["TS"]["labels"][TS_ID] = label_actor

            energies = [f"E={self.project.TS_data[TS_ID]["E_min"]:.2f}kJ/mol"] * len(centers)
            energy_actor = self.plotter.add_point_labels(
                centers + np.array([0, 0, -0.15]),
                energies,
                font_size=9,
                text_color = "darkgray",
                show_points=False,
                shape=None,
                always_visible=True
            )
            energy_actor.SetPickable(False)
            self.actors["TS"]["energies"][TS_ID] = energy_actor
    
        return None

    def update_TS_actors(self):

        view_mode = self.project.view_mode

        for TS_ID in self.project.TS_data:

            color = self.project.plotting_data["TS"][TS_ID]["color"]
            opacity = self.project.plotting_data["TS"][TS_ID]["opacity"]
            visibility = self.project.visibility["individual_TS"][TS_ID]
            label_visibility = self.project.visibility["TS_labels"][TS_ID]
            energy_visibility = self.project.visibility["TS_energies"][TS_ID]

            self.actors["TS"]["volume"][TS_ID].GetProperty().SetColor(*color)
            self.actors["TS"]["volume"][TS_ID].GetProperty().SetOpacity(opacity)
            self.actors["TS"]["volume"][TS_ID].SetVisibility(visibility and (view_mode == "volume"))
                
            self.actors["TS"]["graph"][TS_ID].GetProperty().SetColor(*color)
            self.actors["TS"]["graph"][TS_ID].GetProperty().SetOpacity(opacity)
            self.actors["TS"]["graph"][TS_ID].SetVisibility(visibility and (view_mode == "graph"))

            self.actors["TS"]["labels"][TS_ID].SetVisibility(visibility and label_visibility)

            self.actors["TS"]["energies"][TS_ID].SetVisibility(visibility and energy_visibility)

        return None


    #-----------------------------------------------------
    #-----------------------------------------------------
    # CLICKING LOGIC
    #-----------------------------------------------------
    #-----------------------------------------------------


    def set_click_mode(self, mode):

        if self.click_mode == ClickMode.SELECT:
            self.info_box.clear()

        self.click_mode = mode
        self.clear_click_state()

        return None

    def clear_click_state(self):

        self.selected_objects.clear()
        for actor in self.click_actors:
            self.plotter.remove_actor(actor)
        self.click_actors.clear()

        return None

    def mouse_clicked(self, actor):

        if actor is None:
            return

        kind = self.actor_metadata[actor]["kind"]
        ID = self.actor_metadata[actor]["ID"]
        if kind == "atoms":
            obj = PickedObject(kind, ID, self.project.atoms[ID]["center"], self.PBC_atom_pos[ID])

        elif kind == "basins":
            obj = PickedObject(kind, ID, self.project.basin_data[ID]["center"], self.PBC_basin_pos[ID])
        
        self.dispatch_click(obj)


    def dispatch_click(self, obj):

        if self.click_mode == ClickMode.SELECT:
            self.selection_mode(obj)

        elif self.click_mode == ClickMode.DISTANCE:
            self.distance_mode(obj)

        elif self.click_mode == ClickMode.ANGLE:
            self.angle_mode(obj)

        return None



    def selection_mode(self, obj):

        self.clear_click_state()
        self.selected_objects = [obj]
        self.highlight(obj)
        self.display_object(obj)

        return None



    def distance_mode(self, obj):

        if len(self.selected_objects) == 2:
            self.clear_click_state()

        self.selected_objects.append(obj)
        self.highlight(obj)

        if len(self.selected_objects) == 2:

            self.measure_distance(
                self.selected_objects[0],
                self.selected_objects[1]
            )

    def measure_distance(self, obj1, obj2):

        p1, p2 = self.closest_pbc_positions(
            obj1.copies,
            obj2.copies
        )

        self.draw_distance(p1, p2)

        distance = np.linalg.norm(p1 - p2)

        self.info_box.append(
            f"Distance ({obj1.ID} → {obj2.ID}) = {distance:.3f} Å"
        )

        self.plotter.render()


    def draw_distance(self, p1, p2):

        poly = pv.PolyData()

        poly.points = np.array([p1, p2])

        poly.lines = np.array([2, 0, 1])

        actor = self.plotter.add_mesh(
            poly,
            color="black",
            line_width=3,
            lighting=False,
            reset_camera=False
        )

        actor.SetPickable(False)

        self.click_actors.append(actor)

    def closest_pbc_positions(self, copies1, copies2):

        best1 = None
        best2 = None

        min_dist = np.inf

        for p1 in copies1:

            for p2 in copies2:

                dist = np.linalg.norm(p1 - p2)

                if dist < min_dist:

                    min_dist = dist

                    best1 = p1
                    best2 = p2

        return best1, best2

    def angle_mode(self, obj):

        if len(self.selected_objects) == 3:

            self.clear_click_state()

        self.selected_objects.append(obj)

        self.highlight(obj)

        if len(self.selected_objects) == 3:

            self.measure_angle(
                self.selected_objects
            )

    def measure_angle(self, objects):

        obj1, obj2, obj3 = objects

        p1, p2 = self.closest_pbc_positions(
            obj1.copies,
            obj2.copies
        )

        p2_1, p3 = self.closest_pbc_positions(
            obj2.copies,
            obj3.copies
        )

        v1 = p1 - p2
        v2 = p3 - p2_1

        v1 /= np.linalg.norm(v1)
        v2 /= np.linalg.norm(v2)

        angle = np.degrees(
            np.arccos(
                np.clip(np.dot(v1, v2), -1.0, 1.0)
            )
        )

        self.draw_distance(p2, p1)
        self.draw_distance(p2_1, p3)

        self.info_box.append(
            f"Angle ({obj1.ID} - {obj2.ID} - {obj3.ID}) = {angle:.2f}°"
        )


    def highlight(self, obj):

        if obj.kind == "atoms":
            radius = (
                self.project.plotting_data["atoms"]["radii"][
                    self.project.atoms[obj.ID]["type"]
                ] * 1.2
            )

            points = pv.PolyData(np.array(obj.copies))

            points["radius"] = np.full(len(obj.copies), radius)

            spheres = points.glyph(
                scale="radius",
                geom=pv.Sphere(theta_resolution=24, phi_resolution=24),
                orient=False
            )

            actor = self.plotter.add_mesh(
                spheres,
                color="yellow",
                opacity=0.4,
                reset_camera=False
            )




        elif obj.kind == "basins":
            points = pv.PolyData(np.array(obj.copies))
            poly = pv.PolyData(points)
            actor = self.plotter.add_mesh(
                poly,
                color = "yellow",
                opacity = 0.4,
                point_size=30,
                render_points_as_spheres=True,
                reset_camera = False
            )

        actor.SetPickable(False)
        self.click_actors.append(actor)


    def display_object(self, obj):

        self.info_box.clear()

        if obj.kind == "atoms":

            atom = self.project.atoms[obj.ID]

            pos = cartesian_to_fractional(obj.center, self.project.metadata["origin"], self.project.metadata["grid_vectors"], self.project.metadata["grid_shape"])

            self.info_box.append(f"Atom: {obj.ID}")
            self.info_box.append(f"Element: {atom['type']}")
            self.info_box.append(
                f"Position: ({pos[0]:.2f}, "
                f"{pos[1]:.2f}, "
                f"{pos[2]:.2f})"
            )

        elif obj.kind == "basins":

            basin = self.project.basin_data[obj.ID]

            pos = [basin["center"][0]/self.project.metadata["grid_shape"][0], basin["center"][1]/self.project.metadata["grid_shape"][1], basin["center"][2]/self.project.metadata["grid_shape"][2]]

            self.info_box.append(f"Basin {obj.ID}")
            self.info_box.append(
                f"Position: ({pos[0]:.2f}, "
                f"{pos[1]:.2f}, "
                f"{pos[2]:.2f})"
            )
            self.info_box.append(f"Emin: {basin['E_min']:.2f} kJ/mol")
            self.info_box.append(f"Volume: {basin['V']:.2f} Å³")
            self.info_box.append(f"Area: {basin['A']:.2f} Å²")

