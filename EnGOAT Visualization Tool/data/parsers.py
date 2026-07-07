import numpy as np
import numpy as np
import pyvista as pv
import re
import os
from pathlib import Path
from collections import defaultdict
from ase.data import chemical_symbols

def read_grid_info(logfile):
    """
    Read EnGOAT output.log and return:
      - grid points (Na, Nb, Nc)
      - a-axis vector
      - b-axis vector
      - c-axis vector
    """

    with open(logfile, "r") as f:
        text = f.read()

    grid_match = re.search(
        r"Grid points \(Na,Nb,Nc\):\s*\(\s*(\d+),\s*(\d+),\s*(\d+)\)",
        text
    )

    a_match = re.search(
        r"a-axis vector:\s*\(\s*([-0-9.]+),\s*([-0-9.]+),\s*([-0-9.]+)\)",
        text
    )

    b_match = re.search(
        r"b-axis vector:\s*\(\s*([-0-9.]+),\s*([-0-9.]+),\s*([-0-9.]+)\)",
        text
    )

    c_match = re.search(
        r"c-axis vector:\s*\(\s*([-0-9.]+),\s*([-0-9.]+),\s*([-0-9.]+)\)",
        text
    )

    if not all([grid_match, a_match, b_match, c_match]):
        raise ValueError("Could not find all grid/cell information in file.")

    Na, Nb, Nc = map(int, grid_match.groups())

    a_vec = tuple(map(float, a_match.groups()))
    b_vec = tuple(map(float, b_match.groups()))
    c_vec = tuple(map(float, c_match.groups()))

    return {
        "grid_points": (Na, Nb, Nc),
        "a_vector": a_vec,
        "b_vector": b_vec,
        "c_vector": c_vec,
    }

def read_E_levels(logfile):
    with open(logfile, "r") as f:
        text = f.read()

    N_match = re.search(r"Number of levels:\s*(\d+)", text)

    E_step_match = re.search(
        r"E step:\s*([-+]?\d*\.?\d+)",
        text
    )

    E_cutoff_match = re.search(
        r"E cutoff:\s*([-+]?\d*\.?\d+)",
        text
    )

    if not all([N_match, E_step_match, E_cutoff_match]):
        raise ValueError("Could not find energy level information.")

    return {
        "N_steps": int(N_match.group(1)),
        "E_step": float(E_step_match.group(1)),
        "E_cutoff": float(E_cutoff_match.group(1)),
    }

def get_them_files(main_dir):

    main_dir = Path(main_dir).resolve()

    files = {
        "output": None,
        "cube_file": None,
        "TuTraSt_data": {},
        "Tunnels": [],
        "NumPy_matrices": {}
    }

    # Top-level files
    cube_files = list(main_dir.glob("*.cube"))
    if cube_files:
        files["cube_file"] = str(cube_files[0].resolve())

    output_file = main_dir / "output.log"
    if output_file.exists():
        files["output"] = str(output_file.resolve())

    # TuTraSt_data
    tutrast_dir = main_dir / "TuTraSt_data"
    tunnel_files = {}

    if tutrast_dir.exists():
        for f in tutrast_dir.glob("*.dat"):
            name = f.stem

            m = re.match(r"(tunnel\d+)_data$", name)
            if m:
                tunnel_id = m.group(1)

                tunnel_files[tunnel_id] = {
                    "tunnel_file": str(f.resolve()),
                    "MEP_files": {}   # dict: direction -> file
                }

            else:
                files["TuTraSt_data"][name] = str(f.resolve())

    # Tunnel_data (MEPs)
    tunnel_data_dir = main_dir / "Tunnel_data"

    if tunnel_data_dir.exists():
        for f in tunnel_data_dir.glob("*.dat"):

            # min_E_path_a_tunnel1.dat
            m = re.match(r"min_E_path_([abc])_(tunnel\d+)$", f.stem)

            if m:
                direction = m.group(1)   # a, b, or c
                tunnel_id = m.group(2)

                if tunnel_id not in tunnel_files:
                    tunnel_files[tunnel_id] = {
                        "tunnel_file": None,
                        "MEP_files": {}
                    }

                # store only if exists (no assumption all directions present)
                tunnel_files[tunnel_id]["MEP_files"][direction] = str(f.resolve())

    files["Tunnels"] = list(tunnel_files.values())

    # NumPy matrices
    NumPy_matrices = main_dir / "NumPy_matrices"

    files["NumPy_matrices"] = {
        "Basin_matrix": NumPy_matrices / "Basin_matrix.npy",
        "TS_matrix": NumPy_matrices / "TS_matrix.npy",
        "Level_matrix": NumPy_matrices / "Level_matrix.npy",
        "Tunnel_matrix": NumPy_matrices / "Tunnel_matrix.npy",
    }

    return files

def tunnel_summary(tunnels):

    result = {}

    for tunnel in tunnels:

        # extract tunnel ID from filename
        tunnel_file = tunnel.get("tunnel_file")

        tunnel_id = Path(tunnel_file).stem.replace("_data", "")
        tunnel_id = tunnel_id.removeprefix("tunnel")
        
        key = f"Tunnel system {tunnel_id}"

        mep = tunnel.get("MEP_files", {})

        result[key] = {
            "a": "a" in mep,
            "b": "b" in mep,
            "c": "c" in mep
        }

    return result


def read_atom_info(cube_file):

    atoms = defaultdict(list)

    with open(cube_file, "r") as f:
        lines = f.readlines()

    n_atoms = abs(int(lines[2].split()[0]))

    atom_start = 6
    atom_end = atom_start + n_atoms

    for line in lines[atom_start:atom_end]:
        fields = line.split()

        Z = int(fields[0])  # atomic number
        x, y, z = map(float, fields[2:5])

        symbol = chemical_symbols[Z]
        atoms[symbol].append(np.array((x, y, z))*0.529177210903)

    return dict(atoms)

def read_tunnel_info(filename):

    with open(filename, "r", encoding="utf-8") as f:
        text = f.read()

    # Split into individual tunnel system sections
    blocks = re.split(r"Tunnel system #(\d+)", text)

    tunnel_data = {}

    # blocks = [before, id1, block1, id2, block2, ...]
    for i in range(1, len(blocks), 2):
        tunnel_id = blocks[i]
        block = blocks[i + 1]

        V = re.search(r"Volume:\s+([0-9.+-Ee]+)", block)
        A = re.search(r"Surface area:\s+([0-9.+-Ee]+)", block)
        Emin = re.search(r"Lowest energy point:\s+([0-9.+-Ee]+)", block)
        Dim = re.search(r"Dimensionality:\s+(\d+)", block)

        if V and A and Emin and Dim:
            tunnel_data[f"Tunnel system {tunnel_id}"] = {
                "V": float(V.group(1)),
                "A": float(A.group(1)),
                "E_min": float(Emin.group(1)),
                "Dim": int(Dim.group(1)),
            }

    return tunnel_data