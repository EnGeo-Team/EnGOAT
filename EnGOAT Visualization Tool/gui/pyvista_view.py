from PySide6.QtWidgets import QWidget, QVBoxLayout
from PySide6.QtCore import Qt
from pyvistaqt import QtInteractor
from PySide6.QtWidgets import QTextEdit
from PySide6.QtGui import QTextCursor

from plotting.unit_cell import build_unit_cell, get_unit_cell_params
import numpy as np
import pyvista as pv
from plotting.atoms import element_color, covalent_radius, get_atom_info, get_bond_info
from plotting.Basins_transition_states import get_pbc_centers, get_pbc_ts_lines, get_pbc_measurement, voxel_surface
from datetime import datetime

from scipy.ndimage import binary_erosion

class PyVistaView(QWidget):

    def __init__(self):
        super().__init__()

        #Plotting widget
        layout = QVBoxLayout(self)
        self.plotter = QtInteractor(self)
        layout.addWidget(self.plotter)

        #Info box widget
        self.info_box = QTextEdit()
        self.info_box.setReadOnly(True)
        self.info_box.setPlaceholderText("")#ADD LATER
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
        layout.addWidget(self.info_box)

        self.plotter.enable_depth_peeling()

        self.actors = {
                        "unit_cell": [],
                        "unit_cell_params": [],
                        "atoms": [],
                        "bonds": [],
                        "isosurface": []
                        }        
        

        self.basin_actors = {
            "volume": {},   # basin_id → actor
            "graph": {}
        }

        self.basin_labels = {}   # basin_id → label_actors

        self.TS_actors = {
            "volume": {},   
            "graph": {}
        }

        self.TS_labels = {}  
        
        self.selected_atoms = []

        self.plotter.setFocusPolicy(Qt.StrongFocus)
        self.plotter.setFocus()
        self.plotter.mousePressEvent = self._plotter_mouse_press
        
        self.plotter.show_axes()


    def set_project(self, project):
        self.project = project


    def set_visibility(self, category, visible):
        for actor in self.actors.get(category, []):
            # Text actors
            if hasattr(actor, "SetVisibility"):
                actor.SetVisibility(visible)
            else:
                try:
                    actor.VisibilityOn() if visible else actor.VisibilityOff()
                except:
                    pass

    #Reset button function
    def clear(self):
        self.plotter.clear()
        for key in self.actors:
            self.actors[key].clear()

        self.plotter.show_axes()
        self.plotter.reset_camera()

#
#UNIT CELL
#
    def draw_unit_cell(self, grid):

        corners, edges = build_unit_cell(grid)
        for i, j in edges:
            actor = self.plotter.add_lines(
                np.vstack((corners[i], corners[j])),
                color="black",
                width=3
            )
            actor.pickable = False
            self.actors["unit_cell"].append(actor)
        self.plotter.reset_camera()


        Na, Nb, Nc = grid["grid_points"]
        a_vec = np.array(grid["a_vector"])
        b_vec = np.array(grid["b_vector"])
        c_vec = np.array(grid["c_vector"])

        i = np.arange(Na)
        j = np.arange(Nb)
        k = np.arange(Nc)

        I, J, K = np.meshgrid(i, j, k, indexing="ij")

        self.grid = grid
        self.cartesian_points = (
            I[..., None] * a_vec +
            J[..., None] * b_vec +
            K[..., None] * c_vec
        )


    def display_UC_params(self, grid):

        for actor in self.actors["unit_cell_params"]:
            self.plotter.remove_actor(actor)
        self.actors["unit_cell_params"].clear()

        corners, _ = build_unit_cell(grid)

        a_origin = corners[2]
        b_origin = corners[1]
        c_origin = corners[4]
        origin = corners[0]
        a_end  = corners[1]
        b_end  = corners[2]
        c_end  = corners[3]

        A, B, C, alpha, beta, gamma = get_unit_cell_params(grid)

        a_vec = a_end - origin
        b_vec = b_end - origin
        c_vec = c_end - origin


        ab_corner = a_vec + b_vec

        alpha_pos = ab_corner - 0.15 * (b_vec + c_vec)+0.3*c_vec
        beta_pos  = ab_corner - 0.15 * (a_vec + c_vec)+0.3*c_vec
        gamma_pos = ab_corner - 0.15 * (a_vec + b_vec)


        def offset_point(origin, target, d=1):
            v = origin + target
            v = v / np.linalg.norm(v)
            return target + d * v

        labels = [
            (offset_point(a_origin, a_end)/2 + b_end, f"a = {A:.2f}"),
            (offset_point(b_origin, b_end)/2 + a_end, f"b = {B:.2f}"),
            (offset_point(c_origin, c_end)/2 + a_end + b_end, f"c = {C:.2f}"),
            (alpha_pos,
             rf"$\alpha$ = {alpha:.1f}°"),

            (beta_pos,
             rf"$\beta$ = {beta:.1f}°"),

            (gamma_pos,
             rf"$\gamma$ = {gamma:.1f}°"),

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
            actor.pickable = False
            self.actors["unit_cell_params"].append(actor)

#
#ATOMS AND BONDS
#
    def draw_atoms(self, atoms):

        self.atom_actors = {}
        self.atom_positions = []
        self.atom_elements = []

        for element, positions in atoms.items():
            color = element_color(element)
            radius = covalent_radius(element)

            positions = np.array(positions)
            if len(positions) == 0:
                continue

            colors = np.tile([color], (len(positions), 1))
            radii = np.array([radius]*len(positions))

            poly = pv.PolyData(positions)
            poly["radius"] = radii
            poly["colors"] = colors

            sphere = pv.Sphere(theta_resolution=24, phi_resolution=24)

            glyphs = poly.glyph(scale="radius", geom=sphere, orient=False)

            actor = self.plotter.add_mesh(
                glyphs,
                scalars="colors",
                rgb=True,
                lighting=False,
                opacity=1.0
            )

            actor.pickable = True

            # ✅ store actor by element
            self.atom_actors[element] = actor
            self.actors["atoms"].append(actor)

            # ✅ store for picking
            self.atom_positions.extend(positions)
            self.atom_elements.extend([element] * len(positions))

        self.atom_positions = np.array(self.atom_positions)
        self.atom_elements = np.array(self.atom_elements)

        # ✅ store which elements are currently visible
        self.visible_elements = set(self.atom_actors.keys())

    def update_visible_atoms(self, selected_elements):
    
        self.visible_elements = set(selected_elements)
    
        for element, actor in self.atom_actors.items():
        
            visible = element in self.visible_elements
    
            actor.SetVisibility(visible)
            actor.pickable = visible
    
        self.plotter.render()

    def highlight_basin(self, basin, color="yellow"):

        points = self.basin_centers.get(basin.ID, [])

        if len(points) == 0:
            return None

        poly = pv.PolyData(points)

        sphere = pv.Sphere(radius=0.4)

        glyphs = poly.glyph(
            geom=sphere,
            orient=False,
            scale=False
        )

        actor = self.plotter.add_mesh(
            glyphs,
            color=color,
            opacity=0.6,
            reset_camera=False
        )

        return actor

    def clear_selections(self):
    
        self.selected_atoms.clear()
    
        for name in [
            "highlight_actor_1",
            "highlight_actor_2",
            "distance_actor"
        ]:
            if hasattr(self, name):
                self.plotter.remove_actor(getattr(self, name))
    
    def pick_object(self, point):
    
        # -------------------------
        # Atoms
        # -------------------------
        if hasattr(self, "atom_positions") and self.project.show_atoms:
        
            mask = np.isin(
                self.atom_elements,
                list(self.visible_elements)
            )
    
            if np.any(mask):
            
                visible_positions = self.atom_positions[mask]
                visible_elements = self.atom_elements[mask]
    
                distances = np.linalg.norm(
                    visible_positions - point,
                    axis=1
                )
    
                idx = np.argmin(distances)
    
                if distances[idx] <= 0.8:
                
                    return (
                        visible_elements[idx],
                        visible_positions[idx],
                        "atom"
                    )
    
        # -------------------------
        # Basins
        # -------------------------
        basin_result = self.detect_clicked_basin(
            self.project,
            point
        )
    
        if basin_result is not None:
        
            basin, basin_pos = basin_result
    
            return (
                f"B{basin.ID}",
                basin_pos,
                "basin"
            )
    
        return None
    
    def create_highlight(self, label, pos, obj_type, color):
    
        if obj_type == "basin":
        
            basin_id = int(label.removeprefix("B"))
    
            basin = next(
                b for b in self.project.Basin_list
                if b.ID == basin_id
            )
    
            return self.highlight_basin(
                basin,
                color=color
            )
    
        radius = covalent_radius(label) * 0.5 + 0.1
    
        return self.plotter.add_mesh(
            pv.Sphere(center=pos, radius=radius),
            color=color,
            opacity=0.6,
            reset_camera=False
        )
    
    def get_selection_copies(self, label, pos, obj_type):
    
        if obj_type == "atom":
            return [pos]
    
        basin_id = int(label.removeprefix("B"))
    
        return self.basin_centers[basin_id]

    def closest_pbc_points(self, points1, points2):
        
        min_dist = np.inf
        best_p1 = None
        best_p2 = None

        for p1 in points1:
            for p2 in points2:
            
                dist = np.linalg.norm(p1 - p2)

                if dist < min_dist:
                    min_dist = dist
                    best_p1 = p1
                    best_p2 = p2
        return min_dist, best_p1, best_p2

    def draw_distance_actor(self, segments):
    
        if hasattr(self, "distance_actor"):
            self.plotter.remove_actor(self.distance_actor)
    
        poly = pv.PolyData()
    
        points = []
        lines = []
    
        idx = 0
    
        for s1, s2 in segments:
        
            points.append(s1)
            points.append(s2)
    
            lines.append([2, idx, idx + 1])
    
            idx += 2
    
        poly.points = np.array(points)
        poly.lines = np.array(lines)
    
        self.distance_actor = self.plotter.add_mesh(
            poly,
            color="black",
            line_width=1,
            lighting=False,
            reset_camera=False
        )

    def handle_click(self):

        if not hasattr(self, "atom_positions") and not hasattr(self, "basin_centers"):
            return

        point = self.plotter.pick_mouse_position()

        if point is None:
            return

        picked = self.pick_object(np.array(point))

        if picked is None:
            return

        label, pos, obj_type = picked

        # -------------------------
        # Reset after 2 selections
        # -------------------------
        if len(self.selected_atoms) == 2:
            self.clear_selections()

        # -------------------------
        # Prevent duplicate click
        # -------------------------
        if len(self.selected_atoms) == 1:

            _, prev_pos, _ = self.selected_atoms[0]

            if np.linalg.norm(prev_pos - pos) < 1e-3:
                return

        self.selected_atoms.append(
            (label, pos, obj_type)
        )

        time = datetime.now().strftime("%H:%M:%S")

        # -------------------------
        # First click
        # -------------------------
        if len(self.selected_atoms) == 1:

            self.highlight_actor_1 = self.create_highlight(
                label,
                pos,
                obj_type,
                "yellow"
            )

            info_text = (
                f"[{time}] {label}: "
                f"({pos[0]:.2f}, {pos[1]:.2f}, {pos[2]:.2f})"
            )

        # -------------------------
        # Second click
        # -------------------------
        else:

            self.highlight_actor_2 = self.create_highlight(
                label,
                pos,
                obj_type,
                "orange"
            )

            (l1, p1, type1), (l2, p2, type2) = self.selected_atoms

            copies1 = self.get_selection_copies(
                l1, p1, type1
            )

            copies2 = self.get_selection_copies(
                l2, p2, type2
            )

            _, line_p1, line_p2 = self.closest_pbc_points(
                copies1,
                copies2
            )

            distance, segments = get_pbc_measurement(
                line_p1,
                line_p2,
                self.grid
            )

            self.draw_distance_actor(segments)

            info_text = (
                f"[{time}] {label}: "
                f"({pos[0]:.2f}, {pos[1]:.2f}, {pos[2]:.2f})\n"
                f"[{time}] Distance ({l1} → {l2}) = {distance:.3f} Å"
            )

        self.info_box.append(info_text)

        self.info_box.verticalScrollBar().setValue(
            self.info_box.verticalScrollBar().maximum()
        )

        max_lines = 30

        if self.info_box.document().blockCount() > max_lines:

            cursor = self.info_box.textCursor()

            cursor.movePosition(QTextCursor.Start)
            cursor.select(QTextCursor.LineUnderCursor)
            cursor.removeSelectedText()

        self.plotter.render()


    def _plotter_mouse_press(self, event):

        QtInteractor.mousePressEvent(self.plotter, event)
        if event.button() == Qt.LeftButton:
            self.plotter.render()
            self.handle_click()      
        
        if event.button() == Qt.RightButton:
            self.selected_atoms.clear()

            if hasattr(self, "highlight_actor_1"):
                self.plotter.remove_actor(self.highlight_actor_1)
            if hasattr(self, "highlight_actor_2"):
                self.plotter.remove_actor(self.highlight_actor_2)
            if hasattr(self, "distance_actor"):
                self.plotter.remove_actor(self.distance_actor)

            self.plotter.render()

    def draw_bonds(self, atoms):

        self.bond_actors = {}

        bond_data = get_bond_info(atoms)

        for pair, data in bond_data.items():

            poly = pv.PolyData()
            poly.points = data["points"]
            poly.lines = data["lines"]
            poly["colors"] = data["colors"]

            actor = self.plotter.add_mesh(
                poly,
                scalars="colors",
                rgb=True,
                line_width=3,
                lighting=False
            )

            actor.pickable = False

            self.bond_actors[pair] = actor
            self.actors["bonds"].append(actor)

        # ✅ track visible pairs
        self.visible_bonds = set(self.bond_actors.keys())

    def update_visible_bonds(self, selected_pairs):

        self.visible_bonds = set(selected_pairs)

        for pair, actor in self.bond_actors.items():
            visible = pair in self.visible_bonds
            actor.SetVisibility(visible)

        self.plotter.render()

#
#ISOSURFACE
#
    def draw_isosurface(self, grid, level_matrix, isosurface_params):

        level = isosurface_params["level"]
        color = isosurface_params["color"]
        opacity = isosurface_params["opacity"]

        for actor in self.actors["isosurface"]:
            self.plotter.remove_actor(actor)
        self.actors["isosurface"].clear()

        mask = (level_matrix <= level)

        dx = np.array(grid["a_vector"])
        dy = np.array(grid["b_vector"])
        dz = np.array(grid["c_vector"])

        surface = voxel_surface(
            mask,
            self.cartesian_points,
            dx,
            dy,
            dz
        )

        #actor = self.plotter.add_mesh(
        #    surface,
        #    style="wireframe",
        #    line_width=1,
        #    reset_camera=False,
        #    color=color, 
        #    opacity = opacity
        #)

        #actor = self.plotter.add_silhouette(
        #    surface,
        #    color=color,
        #    line_width=2
        #)

#full volume version

        actor = self.plotter.add_mesh(
            surface,
            color=color,
            opacity=opacity,
            reset_camera=False,
            show_edges=False
        )

        self.actors["isosurface"].append(actor)

        self.plotter.render()

#
#BASINS
#

    def _create_volume_basin_actor(self, b, Basin_matrix):

        mask = (Basin_matrix == b.ID)

        a_vec = np.array(self.grid["a_vector"])
        b_vec = np.array(self.grid["b_vector"])
        c_vec = np.array(self.grid["c_vector"])

        surface = voxel_surface(
            mask,
            self.cartesian_points,
            a_vec,
            b_vec,
            c_vec
        )

        actor = self.plotter.add_mesh(
            surface,
            color=b.color,
            opacity=b.opacity,
            show_edges=False
        )

        return actor

    def _create_graph_basin_actor(self, b):

        points = get_pbc_centers(b, self.grid)

        if len(points) == 0:
            return None

        poly = pv.PolyData(points)

        actor = self.plotter.add_mesh(
            poly,
            color=b.color,
            opacity=b.opacity,
            point_size=20,
            render_points_as_spheres=True
        )
        
        if not hasattr(self, "basin_centers"):
            self.basin_centers = {}

        self.basin_centers[b.ID] = points
    
        actor._basin_id = b.ID

        return actor
    
    def update_basins(self, project):

        mode = project.view_mode
        Basin_matrix = project.Basin_matrix
        grid = project.grid

        # ✅ ensure grid points exist
        if not hasattr(self, "cartesian_points"):
            self.build_cartesian_grid(grid)

        if not hasattr(self, "basin_actors"):
            self.basin_actors = {"volume": {}, "graph": {}}

        for b in project.Basin_list:

            basin_id = b.ID

            # -------------------------
            # ✅ VOLUME
            # -------------------------
            if basin_id not in self.basin_actors["volume"]:

                actor = self._create_volume_basin_actor(b, Basin_matrix)
                if actor:
                    self.basin_actors["volume"][basin_id] = actor

            if basin_id in self.basin_actors["volume"]:
                actor = self.basin_actors["volume"][basin_id]

                actor.GetProperty().SetColor(*b.color)
                actor.GetProperty().SetOpacity(b.opacity)
                actor.SetVisibility(b.visible and mode == "volume")
                actor.pickable = False
                

            # -------------------------
            # ✅ GRAPH
            # -------------------------
            if basin_id not in self.basin_actors["graph"]:

                actor = self._create_graph_basin_actor(b)
                if actor:
                    self.basin_actors["graph"][basin_id] = actor

            if basin_id in self.basin_actors["graph"]:
                actor = self.basin_actors["graph"][basin_id]

                actor.GetProperty().SetColor(*b.color)
                actor.GetProperty().SetOpacity(b.opacity)
                actor.SetVisibility(b.visible and mode == "graph")
                actor.pickable = (b.visible and mode == "graph")

        self.plotter.render()


    def detect_clicked_basin(self, project, point):

        # ✅ Only allow in graph mode
        if project.view_mode != "graph":
            return None

        if not hasattr(self, "basin_centers"):
            return None

        point = np.array(point)

        min_dist = float("inf")
        selected = None
        selected_pos = None

        for b in project.Basin_list:

            # ✅ only visible basins
            if not b.visible:
                continue

            centers = self.basin_centers.get(b.ID, [])

            for c in centers:

                dist = np.linalg.norm(c - point)

                if dist < min_dist:
                    min_dist = dist
                    selected = b
                    selected_pos = c

        # ✅ threshold
        if min_dist < 0.8:
            return selected, selected_pos

        return None

#
#TRANSITIONS
#

    def _create_volume_ts_actor(self, ts, TS_matrix):

        mask = (TS_matrix == ts.ID)

        a_vec = np.array(self.grid["a_vector"])
        b_vec = np.array(self.grid["b_vector"])
        c_vec = np.array(self.grid["c_vector"])

        surface = voxel_surface(
            mask,
            self.cartesian_points,
            a_vec,
            b_vec,
            c_vec
        )

        actor = self.plotter.add_mesh(
            surface,
            color=ts.color,
            opacity=ts.opacity,
            show_edges = False
        )

        return actor


    def _create_graph_ts_actor(self, ts, project):

        segments, _ = get_pbc_ts_lines(ts, project.grid)

        if not segments:
            return None

        points = []
        lines = []
        idx = 0

        for p1, p2 in segments:

            points.append(p1)
            points.append(p2)

            lines.append([2, idx, idx + 1])
            idx += 2

        poly = pv.PolyData()
        poly.points = np.array(points)
        poly.lines = np.array(lines)

        actor = self.plotter.add_mesh(
            poly,
            color=ts.color,
            opacity=ts.opacity,
            line_width=3
        )

        actor.pickable = True
        actor._ts_id = ts.ID

        return actor


    def update_transitions(self, project):

        mode = project.view_mode  # "volume" or "graph"
        TS_matrix = project.TS_matrix
        grid = project.grid

        # ✅ ensure storage exists
        if not hasattr(self, "TS_actors"):
            self.TS_actors = {"volume": {}, "graph": {}}

        if not hasattr(self, "TS_labels"):
            self.TS_labels = {}

        # ✅ build cartesian grid if needed (for volume)
        if not hasattr(self, "cartesian_points"):
            self.build_cartesian_grid(grid)

        for ts in project.TS_list:

            ts_id = ts.ID

            # ==========================================
            # ✅ VOLUMETRIC REPRESENTATION
            # ==========================================
            if ts_id not in self.TS_actors["volume"]:

                actor = self._create_volume_ts_actor(ts, TS_matrix)

                if actor:
                    self.TS_actors["volume"][ts_id] = actor

            if ts_id in self.TS_actors["volume"]:
                actor = self.TS_actors["volume"][ts_id]

                actor.GetProperty().SetColor(*ts.color)
                actor.GetProperty().SetOpacity(ts.opacity)

                actor.SetVisibility(ts.visible and mode == "volume")

            # ==========================================
            # ✅ GRAPH REPRESENTATION (PBC + CLIPPING)
            # ==========================================
            if ts_id not in self.TS_actors["graph"]:

                actor = self._create_graph_ts_actor(ts, project)

                if actor:
                    self.TS_actors["graph"][ts_id] = actor

            if ts_id in self.TS_actors["graph"]:
                actor = self.TS_actors["graph"][ts_id]

                actor.GetProperty().SetColor(*ts.color)
                actor.GetProperty().SetOpacity(ts.opacity)

                actor.SetVisibility(ts.visible and mode == "graph")

        self.plotter.render()

    def initialize_labels(self, project):

        self.basin_labels.clear()
        self.TS_labels.clear()

        a_vec = np.array(self.grid["a_vector"])
        b_vec = np.array(self.grid["b_vector"])
        c_vec = np.array(self.grid["c_vector"])

        # ======================
        # BASINS
        # ======================

        for basin in project.Basin_list:

            points = get_pbc_centers(basin, self.grid)

            labels = [f"B{basin.ID}"] * len(points)

            actor = self.plotter.add_point_labels(
                points + np.array([0, 0, 0.3]),
                labels,
                font_size=10,
                show_points=False,
                shape=None
            )
            actor.SetVisibility(False)
            self.basin_labels[f"{basin.ID}_label"] = actor
            
            labels = [f"E = {basin.E:.2f}"] * len(points)
            
            actor = self.plotter.add_point_labels(
                points + np.array([0, 0, -0.3]),
                labels,
                font_size=9,
                text_color = "gray",
                show_points=False,
                shape=None
            )
            actor.SetVisibility(False)
            self.basin_labels[f"{basin.ID}_energy"] = actor



        for ts in project.TS_list:

            segments, centers = get_pbc_ts_lines(ts, self.grid)

            pos = centers[0]
            
            actor = self.plotter.add_point_labels(
                [pos + np.array([0, 0, 0.3])],
                [f"TS{ts.ID}"],
                font_size=10,
                show_points=False,
                shape=None
            )

            actor.pickable = False
            actor.SetVisibility(False)

            self.TS_labels[f"{ts.ID}_label"] = actor

            actor = self.plotter.add_point_labels(
                [pos + np.array([0, 0, -0.3])],
                [f"E = {ts.E:.2f}"],
                font_size=9,
                text_color = "gray",
                show_points=False,
                shape=None
            )

            actor.pickable = False
            actor.SetVisibility(False)

            self.TS_labels[f"{ts.ID}_energy"] = actor



    def update_basin_labels(self, project):

        for basin in project.Basin_list:

            visible = basin.visible

            label_key = f"{basin.ID}_label"
            energy_key = f"{basin.ID}_energy"

            if label_key in self.basin_labels:
                self.basin_labels[label_key].SetVisibility(
                    visible and project.show_labels_basins
                )

            if energy_key in self.basin_labels:
                self.basin_labels[energy_key].SetVisibility(
                    visible and project.show_energies_basins
                )

    def update_TS_labels(self, project):

        for ts in project.TS_list:

            visible = ts.visible

            label_key = f"{ts.ID}_label"
            energy_key = f"{ts.ID}_energy"

            if label_key in self.TS_labels:
                self.TS_labels[label_key].SetVisibility(
                    visible and project.show_labels_TS
                )

            if energy_key in self.TS_labels:
                self.TS_labels[energy_key].SetVisibility(
                    visible and project.show_energies_TS
                )

    def update_TuTraSt_plots(self, project):

        project.viewer.update_basins(project)
        project.viewer.update_basin_labels(project)
        project.viewer.update_transitions(project)
        project.viewer.update_TS_labels(project)