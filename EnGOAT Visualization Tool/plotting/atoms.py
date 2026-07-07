from ase.data import atomic_numbers, vdw_radii
from ase.data.colors import jmol_colors
from ase.data import atomic_numbers
from mendeleev import element
import numpy as np
from scipy.spatial import cKDTree

def covalent_radius(symbol):
    el = element(symbol)

    if el.covalent_radius is None:
        raise ValueError(f"No covalent radius for {symbol}")

    return el.covalent_radius / 100.0  # pm → Å


def element_color(symbol):
    return tuple(jmol_colors[atomic_numbers[symbol]])

def get_atom_info(atoms):

    points = []
    colors = []
    radii = []
    elements_list = []

    for symbol, positions in atoms.items():
        color = element_color(symbol)
        radius = covalent_radius(symbol)

        for pos in positions:
            points.append(pos)
            elements_list.append(symbol)
            colors.append(color)
            radii.append(radius)   # reasonable size

    points = np.array(points)
    colors = np.array(colors)
    radii = np.array(radii)
    return (points, colors, radii, elements_list)

def get_bond_info(atoms):
    positions = []
    elements = []

    color_lib = {}
    radius_lib = {}

    for symbol, element_positions in atoms.items():
        color_lib[symbol] = element_color(symbol)
        radius_lib[symbol] = covalent_radius(symbol)
        for pos in element_positions:
            positions.append(np.array(pos))
            elements.append(symbol)

    positions = np.array(positions)

    tree = cKDTree(positions)
    pairs = tree.query_pairs(2.5)

    bond_data = {}

    idx = 0

    
    for i, j in pairs:
        
        p1 = positions[i]
        p2 = positions[j]

        direction = p2 - p1
        dist_sq = np.dot(direction, direction)

        if dist_sq < 0.04:
            continue

        r1 = radius_lib[elements[i]]
        r2 = radius_lib[elements[j]]

        cutoff_sq = (r1 + r2 + 0.2) ** 2

        if dist_sq > cutoff_sq:
            continue

        pair = tuple(sorted([elements[i], elements[j]]))

        if pair not in bond_data:
            bond_data[pair] = {
                "points": [],
                "lines": [],
                "colors": [],
                "idx": 0
            }

        entry = bond_data[pair]
        idx = entry["idx"]

        mid = (p1 + p2) / 2

        color1 = color_lib[elements[i]]
        color2 = color_lib[elements[j]]

        # ✅ First segment

        entry["points"].extend([p1, mid])
        entry["lines"].extend([2, idx, idx + 1])
        entry["colors"].extend([color1, color1])

        idx += 2

        # ✅ Second segment

        entry["points"].extend([mid, p2])
        entry["lines"].extend([2, idx, idx + 1])
        entry["colors"].extend([color2, color2])

        idx += 2

        entry["idx"] = idx


    for pair in bond_data:
        bond_data[pair]["points"] = np.array(bond_data[pair]["points"])
        bond_data[pair]["lines"] = np.array(bond_data[pair]["lines"])
        bond_data[pair]["colors"] = np.array(bond_data[pair]["colors"])

    return bond_data
