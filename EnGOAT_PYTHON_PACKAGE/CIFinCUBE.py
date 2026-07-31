import numpy as np
from ase.io import read


BOHR_TO_ANGSTROM = 0.529177249


def replace_cube_atoms(cube_file, cif_file, output_cube):
    """
    Replace the atom block of a cube file using atoms from a CIF file.

    The original atom count and original atom coordinates are completely
    ignored.

    The volumetric data are identified from the grid dimensions.
    """

    # -------------------------------------------------------------
    # Read original cube file
    # -------------------------------------------------------------
    with open(cube_file, "r") as f:
        lines = f.readlines()

    comment1 = lines[0]
    comment2 = lines[1]

    # Read only the origin from the third line.
    # The original atom count is intentionally ignored.
    third_line = lines[2].split()

    origin = np.array(
        third_line[1:4],
        dtype=float
    )

    # Grid information
    grid_lines = lines[3:6]

    grid_shape = np.array(
        [int(line.split()[0]) for line in grid_lines],
        dtype=int
    )

    n_grid_points = int(np.prod(grid_shape))

    print(f"Grid shape: {tuple(grid_shape)}")
    print(f"Expected volumetric values: {n_grid_points}")

    # -------------------------------------------------------------
    # Extract volumetric data
    # -------------------------------------------------------------
    #
    # We do NOT use the original atom count.
    #
    # Instead, all remaining tokens in the file are read as numbers.
    #
    all_data = []

    for line in lines[6:]:
        all_data.extend(line.split())

    all_data = np.asarray(all_data, dtype=float)

    if len(all_data) < n_grid_points:
        raise ValueError(
            f"Not enough volumetric data. "
            f"Expected {n_grid_points} values, "
            f"but found only {len(all_data)}."
        )

    # The volumetric data are the final n_grid_points values.
    #
    # Any original atom lines before them are ignored.
    volumetric_data = all_data[-n_grid_points:]

    # -------------------------------------------------------------
    # Read CIF
    # -------------------------------------------------------------
    atoms = read(cif_file)

    # Wrap atoms into the unit cell
    atoms.wrap()

    # -------------------------------------------------------------
    # Remove duplicate periodic images
    # -------------------------------------------------------------
    frac = atoms.get_scaled_positions()

    rounded = np.round(frac, 8)

    _, unique_idx = np.unique(
        rounded,
        axis=0,
        return_index=True
    )

    atoms = atoms[np.sort(unique_idx)]

    numbers = atoms.get_atomic_numbers()
    positions = atoms.get_positions()

    natoms_new = len(atoms)

    # -------------------------------------------------------------
    # Write new cube file
    # -------------------------------------------------------------
    with open(output_cube, "w") as f:

        # Comments
        f.write(comment1)
        f.write(comment2)

        # Atom count comes exclusively from the CIF
        f.write(
            f"{natoms_new:5d}"
            f"{origin[0]:13.6f}"
            f"{origin[1]:13.6f}"
            f"{origin[2]:13.6f}\n"
        )

        # Preserve original grid vectors
        for line in grid_lines:
            f.write(line)

        # Write atom data exclusively from CIF
        for Z, pos in zip(numbers, positions):

            f.write(
                f"{Z:5d}"
                f"{0.0:13.6f}"
                f"{pos[0] / BOHR_TO_ANGSTROM:13.6f}"
                f"{pos[1] / BOHR_TO_ANGSTROM:13.6f}"
                f"{pos[2] / BOHR_TO_ANGSTROM:13.6f}\n"
            )

        # Write exactly the expected number of volumetric values
        for i, value in enumerate(volumetric_data, start=1):

            f.write(f"{value:13.5E}")

            # Standard cube formatting: six values per line
            if i % 6 == 0:
                f.write("\n")

        if len(volumetric_data) % 6 != 0:
            f.write("\n")

    print(f"Done! Wrote {output_cube}")
    print(f"Atoms written: {natoms_new}")
    print(f"Volumetric values written: {len(volumetric_data)}")
