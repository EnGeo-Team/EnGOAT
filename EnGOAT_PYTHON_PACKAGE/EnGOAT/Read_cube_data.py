import numpy as np
from ase.data import chemical_symbols
from collections import defaultdict

def Read_cube_data(E_unit, cube_file):
    with open(cube_file) as f:
        f.readline()
        f.readline()

        # Number of atoms and origin
        header = f.readline().split()
        n_atoms = abs(int(float(header[0])))
        origin = np.array(header[1:], dtype=float) * 0.529177249

        # Grid
        grid_shape = np.zeros(3, dtype=int)
        grid_vectors = np.zeros((3, 3))

        for i in range(3):
            line = f.readline().split()
            grid_shape[i] = int(line[0])
            grid_vectors[i] = np.array(line[1:], dtype=float) * 0.529177249

        # Atoms
        atom_data = {}
        atom_counts = defaultdict(int)

        for _ in range(n_atoms):
            line = f.readline().split()

            Z = int(line[0])
            symbol = chemical_symbols[Z]

            atom_counts[symbol] += 1
            atom_ID = f"{symbol}{atom_counts[symbol]}"

            center = tuple(np.array(line[2:5], dtype=float) * 0.529177249)

            atom_data[atom_ID] = {
                "Z": Z,
                "type": symbol,
                "center": center,
            }

    voxel_size = np.linalg.norm(grid_vectors, axis=1)

    grid = {
        "shape": grid_shape,
        "vectors": grid_vectors,
        "voxel_size": voxel_size,
        "origin": origin,
    }

    conversion = {
        1: 1.0,         # kJ/mol
        2: 4.184,       # kcal/mol
        3: 1312.7497,   # Ry
        4: 96.4853,     # eV
        5: 627.5096,    # Hartree
    }[E_unit]

    volumetric_data = np.loadtxt(cube_file, skiprows=6 + n_atoms)

    shift_volumetric_data = (
        volumetric_data - volumetric_data.min()
    ) * conversion

    Energy_matrix = shift_volumetric_data.reshape(grid_shape)

    return atom_data, grid, Energy_matrix