import numpy as np
import pyvista as pv
import re
from matplotlib import colormaps
from itertools import product
import matplotlib.pyplot as plt
from pathlib import Path

class Basin:
    def __init__(self, ID, tunnel, center, E, V, A):
        self.ID = ID
        self.tunnel = tunnel
        self.center = center
        self.E = E
        self.V = V
        self.A = A

        self.visible = False
        self.color =  (0.059, 0.322, 0.729)
        self.opacity = 0.5

    def __repr__(self):
        return (
            f"Basin(ID={self.ID}, tunnel={self.tunnel}, "
            f"center={self.center}, E={self.E}, V={self.V}, A={self.A})"
        )

class Transition_state:
    def __init__(self, ID, B_start, B_end, E, cross_vector):
        self.ID = ID
        self.B_start = B_start
        self.B_end = B_end
        self.E = E
        self.cross_vector = cross_vector
        self.startpoint = None
        self.endpoint = None
        self.tunnel = None

        self.visible = False
        self.color =  (0.059, 0.322, 0.729)
        self.opacity = 0.5
    
    def __repr__(self):
        return (
            f"TS(ID={self.ID}, "
            f"B_start={self.B_start}, "
            f"B_end={self.B_end}, "
            f"E={self.E}, "
            f"cross_vector={self.cross_vector})"
        )

def Make_Basin_TS_lists(basin_file, TS_file):
    Basin_data = basin_file

    pattern = re.compile(
        r"^\s*(\d+)\s+"                      # Basin ID
        r"(\d+|/)\s+"                        # Tunnel system or /
        r"\(\s*(-?\d+)\s*,\s*(-?\d+)\s*,\s*(-?\d+)\s*\)\s+"  # center
        r"([-+]?\d*\.?\d+)\s+"              # E_min
        r"([-+]?\d*\.?\d+)\s+"              # Volume
        r"([-+]?\d*\.?\d+)"                 # Area
    )
    Basin_list = []
    with open(Basin_data, "r") as f:
        for line in f:
            match = pattern.match(line)
            if match:
                ID = int(match.group(1))

                tunnel_str = match.group(2)
                tunnel = None if tunnel_str == "/" else int(tunnel_str)

                center = (
                    int(match.group(3)),
                    int(match.group(4)),
                    int(match.group(5))
                )

                E = float(match.group(6))
                V = float(match.group(7))
                A = float(match.group(8))

                Basin_list.append(Basin(ID, tunnel, center, E, V, A))

    TS_data = TS_file
    
    pattern = re.compile(
        r"^\s*(\d+)\s+"                      # TS ID
        r"([-+]?\d*\.?\d+)\s+"              # E_min
        r"(\d+)\s+"                         # Basin 1
        r"(-?\d+)\s+"                       # Basin 2
        r"\[\s*(-?\d+)\s+(-?\d+)\s+(-?\d+)\s*\]"  # Cross vector
    )
    TS_list =[]
    with open(TS_data, "r") as f:
        for line in f:
            match = pattern.match(line)

            if match:
                ID = int(match.group(1))
                E = float(match.group(2))
                B_start = int(match.group(3))
                B_end = int(match.group(4))

                cross_vector = (
                    int(match.group(5)),
                    int(match.group(6)),
                    int(match.group(7))
                )

                TS_list.append(
                    Transition_state(
                        ID,
                        B_start,
                        B_end,
                        E,
                        cross_vector
                    )
                )
    for TS in TS_list:
        B_start = TS.B_start
        if int(TS.B_end) == -1:
            B_end = B_start
            TS.B_end = B_start
        else:
            B_end = TS.B_end
        for B in Basin_list:
            if B.ID == B_start:
                TS.startpoint = B.center
                TS.tunnel = B.tunnel
            if B.ID == B_end:
                TS.endpoint = B.center
                TS.tunnel = B.tunnel

    return (Basin_list, TS_list)

def read_iso_groups(filename):
    clusters = {}

    with open(filename, "r") as f:
        lines = f.readlines()

    cluster_num = 1
    family_to_name = {}

    # Skip header
    for line in lines[1:]:
        line = line.strip()
        if not line:
            continue

        parts = line.split()

        # Transition state ID (3rd column)
        ts_id = int(parts[2])

        # Extract basin family
        match = re.search(r"\{([^}]*)\}", line)
        if not match:
            continue

        basins = tuple(sorted(
            int(x.strip())
            for x in match.group(1).split(",")
        ))

        # Give each unique basin family a cluster name
        if basins not in family_to_name:
            name = f"Isolated cluster {cluster_num}"
            family_to_name[basins] = name
            clusters[name] = {
                "basins": list(basins),
                "transitions": set()
            }
            cluster_num += 1

        # Add transition state
        clusters[family_to_name[basins]]["transitions"].add(ts_id)

    # Convert transition sets to sorted lists
    for cluster in clusters.values():
        cluster["transitions"] = sorted(cluster["transitions"])

    return clusters

import re

def read_tunnel_system(filename):
    basins = []
    transitions = set()

    with open(filename, "r") as f:
        lines = f.readlines()

    for line in lines:
        line = line.strip()
        if not line:
            continue

        # Basin family
        if line.startswith("{") and line.endswith("}"):
            basins = sorted(
                int(x.strip())
                for x in line.strip("{}").split(",")
            )
            continue

        # Skip comments/header lines
        if line.startswith("#"):
            continue

        parts = line.split()

        # Process lines start with two basin IDs
        if len(parts) >= 3:
            try:
                ts_id = int(parts[2])
                transitions.add(ts_id)
            except ValueError:
                pass

    return {
        "basins": basins,
        "transitions": sorted(transitions)
    }

def read_MEP(filename):
    data = np.loadtxt(filename, skiprows=2)
    data = np.atleast_2d(data)
    MEP = {}
    MEP["basin_ids"] = data[:, 0].astype(int)
    MEP["ts_ids"]    = data[:, 2].astype(int)
    MEP["basin_E"]   = data[:, 3]
    MEP["TS_E"]        = data[:, 5]
    MEP["basin_coords"] = data[:, 6]
    MEP["PBC_crossing"] = data[:, 8]

    # -------------------------
    # MEP plotting colors
    # -------------------------

    n = len(MEP["basin_ids"])

    cmap = plt.get_cmap("tab10")

    MEP["basin_colors"] = [
        cmap(i % 10)[:3]
        for i in range(n)
    ]

    MEP["TS_color"] = (0.0, 0.0, 0.0)
    MEP["opacity"] = 0.8

    return MEP

def Organize_TuTraSt(files):
    basin_file = files["TuTraSt_data"]["basin_data"]
    TS_file = files["TuTraSt_data"]["TS_data"]
    Basin_list, TS_list = Make_Basin_TS_lists(basin_file, TS_file)

    iso_file = files["TuTraSt_data"]["isolated_processes"]
    if Path(iso_file).exists():
        isolated_clusters = read_iso_groups(iso_file)
    else:
        isolated_clusters = None

    clustered_basins = set()
    for cluster in isolated_clusters.values():
        clustered_basins.update(cluster["basins"])

    # Add completely isolated basins (no tunnel and not already clustered)
    next_cluster = len(isolated_clusters) + 1

    for basin in Basin_list:
        if basin.tunnel is None:
            basin.color = (0.839, 0.255, 0.243)

    for TS in TS_list:
        if TS.tunnel is None:
            TS.color = (0.839, 0.255, 0.243)

    for basin in Basin_list:
        if basin.tunnel is None and basin.ID not in clustered_basins:
            isolated_clusters[f"Isolated cluster {next_cluster}"] = {
                "basins": [basin.ID],
                "transitions": []
            }
            next_cluster += 1


    tunnel_systems = {}
    tunnel_systems_plotting = {}
    for entry in files["Tunnels"]:
        
        tunnel_file = entry.get("tunnel_file")
        tunnel_id = Path(tunnel_file).stem.replace("_data", "")
        tunnel_id = tunnel_id.removeprefix("tunnel")
        key = f"Tunnel system {tunnel_id}"

        info = read_tunnel_system(tunnel_file)
        meps = entry.get("MEP_files")
        MEP_dict = {"a":None, "b": None, "c":None}
        for direction, mep_file in meps.items():
            MEP_dict[direction] = read_MEP(mep_file)
        
        info["MEPs"] = MEP_dict
        tunnel_systems[key] = info

        info_plotting = {}

        basins_plotting = {}
        TS_plotting = {}

        basins_plotting["visible"] = False
        basins_plotting["color"] = (0.059, 0.322, 0.729)
        basins_plotting["opacity"] = 0.50
        info_plotting["basins"] = basins_plotting

        TS_plotting["visible"] = False
        TS_plotting["color"] = (0.059, 0.322, 0.729)
        TS_plotting["opacity"] = 0.50

        info_plotting["show_MEP"] = False

        # default direction
        for direction in ["a", "b", "c"]:
            if MEP_dict[direction]:
                info_plotting["MEP_direction"] = direction
                break

        info_plotting["TS"] = TS_plotting
        tunnel_systems_plotting[key] = info_plotting

    isolated_clusters_plotting = {}

    basins_plotting = {}
    TS_plotting = {}
    basins_plotting["visible"] = False
    basins_plotting["color"] = (0.839, 0.255, 0.243)
    basins_plotting["opacity"] = 0.50
    TS_plotting["visible"] = False
    TS_plotting["color"] = (0.839, 0.255, 0.243)
    TS_plotting["opacity"] = 0.50
    isolated_clusters_plotting["basins"] = basins_plotting
    isolated_clusters_plotting["TS"] = TS_plotting

    return (Basin_list, TS_list, tunnel_systems, isolated_clusters, tunnel_systems_plotting, isolated_clusters_plotting)





















def get_pbc_centers(basin, grid):

    Na, Nb, Nc = grid["grid_points"]

    a_vec = np.array(grid["a_vector"])
    b_vec = np.array(grid["b_vector"])
    c_vec = np.array(grid["c_vector"])

    Na_i = Na
    Nb_i = Nb
    Nc_i = Nc

    center = np.array(basin.center)

    # -------------------------
    # ✅ determine shifts
    # -------------------------
    x_shifts = [0]
    if center[0] == 0:
        x_shifts.append(Na_i)
    elif center[0] == Na_i:
        x_shifts.append(-Na_i)

    y_shifts = [0]
    if center[1] == 0:
        y_shifts.append(Nb_i)
    elif center[1] == Nb_i:
        y_shifts.append(-Nb_i)

    z_shifts = [0]
    if center[2] == 0:
        z_shifts.append(Nc_i)
    elif center[2] == Nc_i:
        z_shifts.append(-Nc_i)

    # ✅ generate all images
    points = []

    for shift in product(x_shifts, y_shifts, z_shifts):

        shifted = center + np.array(shift)

        cart = (
            shifted[0] * a_vec +
            shifted[1] * b_vec +
            shifted[2] * c_vec
        )

        points.append(cart)

    return np.array(points)


def get_pbc_ts_lines(ts, grid):

    Na, Nb, Nc = grid["grid_points"]

    a_vec = np.array(grid["a_vector"])
    b_vec = np.array(grid["b_vector"])
    c_vec = np.array(grid["c_vector"])

    Na_i = Na - 1
    Nb_i = Nb - 1
    Nc_i = Nc - 1

    def grid2cart(p):
        return (
            p[0] * a_vec +
            p[1] * b_vec +
            p[2] * c_vec
        )

    start = np.array(ts.startpoint)
    end = np.array(ts.endpoint)

    lines = []

    # -----------------------------
    # No PBC crossing
    # -----------------------------
    if ts.cross_vector == (0, 0, 0):

        p1 = grid2cart(start)
        p2 = grid2cart(end)

        lines.append([p1, p2])

        midpoint = 0.5 * (p1 + p2)

        return lines, [midpoint]

    # -----------------------------
    # PBC crossing
    # -----------------------------
    shift = np.array([
        ts.cross_vector[0] * Na_i,
        ts.cross_vector[1] * Nb_i,
        ts.cross_vector[2] * Nc_i
    ])

    mins = np.array([0, 0, 0])
    maxs = np.array([Na_i, Nb_i, Nc_i])

    clipped_segments = []

    clipped = clip_to_box(
        start,
        end + shift,
        mins,
        maxs
    )

    if clipped is not None:

        q1, q2 = clipped

        p1 = grid2cart(q1)
        p2 = grid2cart(q2)

        lines.append([p1, p2])

        clipped_segments.append((p1, p2))

    clipped = clip_to_box(
        start - shift,
        end,
        mins,
        maxs
    )

    if clipped is not None:

        q1, q2 = clipped

        p1 = grid2cart(q1)
        p2 = grid2cart(q2)

        lines.append([p1, p2])

        clipped_segments.append((p1, p2))

    # -----------------------------
    # Find longest segment(s)
    # -----------------------------
    centers = []

    if clipped_segments:

        lengths = [
            np.linalg.norm(p2 - p1)
            for p1, p2 in clipped_segments
        ]

        max_length = max(lengths)

        tol = 1e-6

        for (p1, p2), L in zip(clipped_segments, lengths):

            if abs(L - max_length) < tol:
                centers.append(0.5 * (p1 + p2))

    return lines, centers

#Clips the graph vertice to the box dimensions
def clip_to_box(p1, p2, mins, maxs):
    """
    Clip line segment p1->p2 to axis-aligned box.

    Returns:
        (q1, q2) if visible
        None      otherwise
    """

    p1 = np.asarray(p1, dtype=float)
    p2 = np.asarray(p2, dtype=float)

    d = p2 - p1

    t0 = 0.0
    t1 = 1.0

    for i in range(3):

        if abs(d[i]) < 1e-12:

            if p1[i] < mins[i] or p1[i] > maxs[i]:
                return None

        else:

            t_enter = (mins[i] - p1[i]) / d[i]
            t_exit  = (maxs[i] - p1[i]) / d[i]

            if t_enter > t_exit:
                t_enter, t_exit = t_exit, t_enter

            t0 = max(t0, t_enter)
            t1 = min(t1, t_exit)

            if t0 > t1:
                return None

    q1 = p1 + t0*d
    q2 = p1 + t1*d

    return q1, q2


def point_inside_unit_cell(cart_point, grid):
    Na, Nb, Nc = grid["grid_points"]

    a = np.array(grid["a_vector"])*Na
    b = np.array(grid["b_vector"])*Nb
    c = np.array(grid["c_vector"])*Nc

    M = np.column_stack([a, b, c])

    frac = np.linalg.solve(M, cart_point)

    return np.all(frac >= 0.0) and np.all(frac <= 1.0)


#For picking
def get_pbc_measurement(p1_cart, p2_cart, grid):


    Na, Nb, Nc = grid["grid_points"]

    a = np.array(grid["a_vector"])*Na
    b = np.array(grid["b_vector"])*Nb
    c = np.array(grid["c_vector"])*Nc

    M = np.column_stack([a, b, c])

    Na_i = Na - 1
    Nb_i = Nb - 1
    Nc_i = Nc - 1

    # -------------------------
    # Cartesian -> grid coords
    # -------------------------
    p1_frac = np.linalg.solve(M, p1_cart)
    p2_frac = np.linalg.solve(M, p2_cart)

    p1 = np.array([
        p1_frac[0] * Na_i,
        p1_frac[1] * Nb_i,
        p1_frac[2] * Nc_i
    ])

    p2 = np.array([
        p2_frac[0] * Na_i,
        p2_frac[1] * Nb_i,
        p2_frac[2] * Nc_i
    ])

    # -------------------------
    # Find shortest PBC image
    # -------------------------
    delta = p2 - p1

    cross_vector = np.zeros(3, dtype=int)

    for i, size in enumerate([Na_i, Nb_i, Nc_i]):

        if delta[i] > size / 2:
            delta[i] -= size
            cross_vector[i] = -1

        elif delta[i] < -size / 2:
            delta[i] += size
            cross_vector[i] = 1

    p2_image = p1 + delta

    # -------------------------
    # PBC distance
    # -------------------------
    d_frac = np.array([
        delta[0] / Na_i,
        delta[1] / Nb_i,
        delta[2] / Nc_i
    ])

    d_cart = (
        d_frac[0] * a +
        d_frac[1] * b +
        d_frac[2] * c
    )

    distance = np.linalg.norm(d_cart)

    mins = np.array([0, 0, 0])
    maxs = np.array([Na_i, Nb_i, Nc_i])

    segments = []

    if np.all(cross_vector == 0):

        q1 = (
            p1[0]/Na_i*a +
            p1[1]/Nb_i*b +
            p1[2]/Nc_i*c
        )

        q2 = (
            p2[0]/Na_i*a +
            p2[1]/Nb_i*b +
            p2[2]/Nc_i*c
        )

        segments.append([q1, q2])

    else:

        clipped = clip_to_box(
            p1,
            p2_image,
            mins,
            maxs
        )

        if clipped is not None:

            q1, q2 = clipped

            q1_cart = (
                q1[0]/Na_i*a +
                q1[1]/Nb_i*b +
                q1[2]/Nc_i*c
            )

            q2_cart = (
                q2[0]/Na_i*a +
                q2[1]/Nb_i*b +
                q2[2]/Nc_i*c
            )

            segments.append([q1_cart, q2_cart])

        shift = np.array([
            cross_vector[0] * Na_i,
            cross_vector[1] * Nb_i,
            cross_vector[2] * Nc_i
        ])

        clipped = clip_to_box(
            p1 - shift,
            p2,
            mins,
            maxs
        )

        if clipped is not None:

            q1, q2 = clipped

            q1_cart = (
                q1[0]/Na_i*a +
                q1[1]/Nb_i*b +
                q1[2]/Nc_i*c
            )

            q2_cart = (
                q2[0]/Na_i*a +
                q2[1]/Nb_i*b +
                q2[2]/Nc_i*c
            )

            segments.append([q1_cart, q2_cart])

    return distance, segments




def voxel_surface(mask, points, da, db, dc):

    verts = []
    faces = []

    face_dirs = [
        (( 1, 0, 0), np.array([[ 0.5, -0.5, -0.5],
                               [ 0.5,  0.5, -0.5],
                               [ 0.5,  0.5,  0.5],
                               [ 0.5, -0.5,  0.5]])),

        ((-1, 0, 0), np.array([[-0.5, -0.5, -0.5],
                               [-0.5, -0.5,  0.5],
                               [-0.5,  0.5,  0.5],
                               [-0.5,  0.5, -0.5]])),

        (( 0, 1, 0), np.array([[-0.5,  0.5, -0.5],
                               [-0.5,  0.5,  0.5],
                               [ 0.5,  0.5,  0.5],
                               [ 0.5,  0.5, -0.5]])),

        (( 0,-1, 0), np.array([[-0.5, -0.5, -0.5],
                               [ 0.5, -0.5, -0.5],
                               [ 0.5, -0.5,  0.5],
                               [-0.5, -0.5,  0.5]])),

        (( 0, 0, 1), np.array([[-0.5, -0.5,  0.5],
                               [ 0.5, -0.5,  0.5],
                               [ 0.5,  0.5,  0.5],
                               [-0.5,  0.5,  0.5]])),

        (( 0, 0,-1), np.array([[-0.5, -0.5, -0.5],
                               [-0.5,  0.5, -0.5],
                               [ 0.5,  0.5, -0.5],
                               [ 0.5, -0.5, -0.5]])),
    ]

    nx, ny, nz = mask.shape

    for i in range(nx):
        for j in range(ny):
            for k in range(nz):

                if not mask[i, j, k]:
                    continue

                center = points[i, j, k]

                for (di, dj, dk), corners in face_dirs:

                    ii = i + di
                    jj = j + dj
                    kk = k + dk

                    exposed = (
                        ii < 0 or ii >= nx or
                        jj < 0 or jj >= ny or
                        kk < 0 or kk >= nz or
                        not mask[ii, jj, kk]
                    )

                    if not exposed:
                        continue

                    base = len(verts)

                    quad = (
                        center
                        + corners[:, [0]] * da
                        + corners[:, [1]] * db
                        + corners[:, [2]] * dc
                    )

                    verts.extend(quad)

                    faces.extend([
                        4,
                        base,
                        base + 1,
                        base + 2,
                        base + 3
                    ])

    poly = pv.PolyData(
        np.asarray(verts),
        np.asarray(faces)
    )

    poly = poly.clean()

    return poly
