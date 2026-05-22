from pathlib import Path
import re
import ast


class Basin:
    def __init__(self, ID, ts, center, E, V, V_boltz):
        self.ID = ID
        self.ts = ts
        self.center = center
        self.E = E
        self.V = V
        self.V_boltz = V_boltz

    def __repr__(self):
        return (f"Basin(id={self.ID}, energy={self.E}, "
                f"center={self.center})")


def read_basins(file_path):

    basins = []

    with open(file_path, "r") as f:
        lines = f.readlines()

    for line in lines:

        # skip header or empty lines
        if not line.strip():
            continue
        if not line.strip()[0].isdigit():
            continue

        # ---------------------------
        # Extract data
        # ---------------------------

        # ID
        parts = line.split()
        basin_id = int(parts[0])

        # Tunnel system
        tunnel = parts[1]
        if tunnel == "/":
            tunnel = None

        # Center (extract using regex)
        center_match = re.search(r"\((.*?)\)", line)
        center_str = center_match.group(1)
        center = tuple(int(x.strip()) for x in center_str.split(","))

        # Extract floats AFTER the center
        after_center = line.split(")")[-1].split()

        energy = float(after_center[0])
        volume = float(after_center[1])
        v_boltz = float(after_center[2])

        basin = Basin(
            basin_id,
            tunnel,
            center,
            energy,
            volume,
            v_boltz
        )

        basins.append(basin)

    return basins


def read_merges(file_path):

    merges = []

    with open(file_path, "r") as f:
        lines = f.readlines()

    for line in lines:

        # skip header / empty lines
        if not line.strip():
            continue
        if not line.strip()[0].isdigit():
            continue

        parts = line.split()

        # columns:
        # ID | Energy | Basin1 | Basin2
        energy = float(parts[1])
        b1 = int(parts[2])
        b2 = int(parts[3])

        merges.append((str(b1), str(b2), energy))

    return merges

def read_cell_data(filename):
    cell_data = {}

    with open(filename, 'r') as f:
        text = f.read()

    # Grid points
    grid_match = re.search(
        r"Grid points \(Na,Nb,Nc\):\s*\(\s*(\d+),\s*(\d+),\s*(\d+)\)",
        text
    )
    if grid_match:
        cell_data['grid_size'] = tuple(map(int, grid_match.groups()))

    # Axis vectors
    def extract_vector(label):
        match = re.search(
            rf"{label}-axis vector:\s*\(\s*([-\d\.]+),\s*([-\d\.]+),\s*([-\d\.]+)\)",
            text
        )
        if match:
            return tuple(map(float, match.groups()))
        return None
    cell_data['b'] = extract_vector('b')
    cell_data['c'] = extract_vector('c')
    cell_data['a'] = extract_vector('a')
    return cell_data

def read_processes(filepath):
    processes = []

    with open(filepath, "r") as f:
        # skip first 5 lines
        for _ in range(5):
            next(f)

        for line in f:
            parts = re.findall(r"\[.*?\]|\(.*?\)|\S+", line)

            if len(parts) < 7:
                continue

            process = {
                "B1": int(parts[0]),
                "B2": int(parts[1]),
                "TS": int(parts[2]),
                "E": float(parts[3]),
                "start_point": ast.literal_eval(parts[4]),
                "end_point": ast.literal_eval(parts[5]),
                "process_vector": ast.literal_eval(parts[6]),
            }

            processes.append(process)
        
    TS_indices = set()
    unique_processes = []
    for process in processes:
        if process["TS"] not in TS_indices:
            TS_indices.add(process["TS"])
            unique_processes.append(process)

    return unique_processes