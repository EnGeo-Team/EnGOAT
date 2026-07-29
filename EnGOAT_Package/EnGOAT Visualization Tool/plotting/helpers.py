from ase.data import atomic_numbers
from ase.data.colors import jmol_colors
from ase.data import atomic_numbers
from mendeleev import element
import numpy as np
import pyvista as pv
from itertools import product

def covalent_radius(symbol):
    el = element(symbol)

    if el.covalent_radius is None:
        raise ValueError(f"No covalent radius for {symbol}")

    return el.covalent_radius / 100.0  # pm → Å


def element_color(symbol):
    return tuple(jmol_colors[atomic_numbers[symbol]])

def create_bond_data(atoms, grid_shape, grid_vectors, atom_radius_data):

    Na, Nb, Nc = grid_shape

    A = Na * np.asarray(grid_vectors[0])
    B = Nb * np.asarray(grid_vectors[1])
    C = Nc * np.asarray(grid_vectors[2])

    M = np.column_stack((A, B, C))
    M_inv = np.linalg.inv(M)

    atom_ids = list(atoms.keys())

    positions = np.array([
        atoms[atom_id]["center"]
        for atom_id in atom_ids
    ])

    elements = [
        atoms[atom_id]["type"]
        for atom_id in atom_ids
    ]

    bond_data = {}

    bond_idx = 1
    n_atoms = len(atom_ids)

    for i in range(n_atoms):

        for j in range(i + 1, n_atoms):

            p1 = positions[i]
            p2 = positions[j]

            dr = p2 - p1

            # Cartesian -> fractional
            frac = M_inv @ dr

            # Minimum-image shift
            shift = -np.round(frac).astype(int)

            # Apply minimum image convention
            frac_min = frac + shift
            dr_min = M @ frac_min

            dist = np.linalg.norm(dr_min)

            r1 = atom_radius_data[elements[i]]
            r2 = atom_radius_data[elements[j]]

            if dist > (r1 + r2 + 0.2):
                continue

            bond_id = f"Bond{bond_idx}"

            bond_data[bond_id] = {
                "type": tuple(sorted([
                    elements[i],
                    elements[j]
                ])),
                "atom_types": (
                    elements[i],
                    elements[j]
                ),
                "atoms": (
                    atom_ids[i],
                    atom_ids[j]
                ),
                "fractional_shift": tuple(
                    int(x) for x in shift
                ),
                "bond_vector": dr_min,
                "distance": dist
            }

            bond_idx += 1

    return bond_data



def pbc_images(cart_coord, grid_shape, grid_vectors, supercell, tol=1e-8):
    """
    Return all periodic images of a Cartesian point inside a supercell.

    Parameters
    ----------
    cart_coord : array_like, shape (3,)
        Cartesian coordinates.
    grid_shape : (Na, Nb, Nc)
        Number of voxels along each lattice vector.
    grid_vectors : (3,3)
        Voxel vectors from metadata["grid_vectors"].
    supercell : (Sa, Sb, Sc)
        Number of repeated unit cells along each lattice vector.
        (1,1,1) reproduces the original behavior.
    tol : float
        Tolerance for detecting points on unit-cell boundaries.

    Returns
    -------
    list[np.ndarray]
        Cartesian coordinates of all unique equivalent images inside the
        requested supercell.
    """
    Na, Nb, Nc = grid_shape
    Sa, Sb, Sc = map(int, supercell)

    # Unit-cell lattice vectors
    a = Na * np.asarray(grid_vectors[0])
    b = Nb * np.asarray(grid_vectors[1])
    c = Nc * np.asarray(grid_vectors[2])

    M = np.column_stack((a, b, c))

    # Fractional coordinates in the unit cell
    frac = np.linalg.solve(M, cart_coord)
    frac %= 1.0

    options = []
    for f, S in zip(frac, (Sa, Sb, Sc)):
        if np.isclose(f, 0.0, atol=tol):
            # Boundary point: include all equivalent copies
            vals = [k for k in range(S + 1)]
        else:
            # Interior point: one copy per repeated unit cell
            vals = [f + k for k in range(S)]
        options.append(vals)

    images = []
    seen = set()

    for u, v, w in product(*options):
        frac_super = np.array([u, v, w])
        cart = M @ frac_super

        # Remove duplicates caused by tolerance
        key = tuple(np.round(cart / tol).astype(int))
        if key not in seen:
            seen.add(key)
            images.append(cart)

    return images



def all_bonds_inside_UC(
    data,
    metadata,
    atoms,
    supercell,
    tol=1e-8,
):
    """
    Return all unique periodic copies of a bond whose endpoints lie inside
    the requested supercell.

    Parameters
    ----------
    data : dict
        Bond data.
    metadata : dict
        Grid metadata.
    atoms : dict
        Atomic information.
    supercell : (Sa, Sb, Sc)
        Number of repeated unit cells along each lattice vector.
        (1,1,1) reproduces the original behavior.
    tol : float

    Returns
    -------
    list[(np.ndarray, np.ndarray)]
        Cartesian coordinates of all bond copies.
    """

    atom1, atom2 = data["atoms"]

    p1 = np.asarray(atoms[atom1]["center"], dtype=float)
    p2 = p1 + np.asarray(data["bond_vector"], dtype=float)

    origin = np.asarray(metadata["origin"], dtype=float)
    grid_vectors = np.asarray(metadata["grid_vectors"], dtype=float)
    grid_shape = np.asarray(metadata["grid_shape"], dtype=float)

    Sa, Sb, Sc = map(int, supercell)

    # Unit-cell lattice vectors
    cell_vectors = grid_shape[:, None] * grid_vectors

    # Fractional coordinates
    f1 = np.linalg.solve(cell_vectors.T, p1 - origin)
    f2 = np.linalg.solve(cell_vectors.T, p2 - origin)

    # Require both endpoints to lie inside the supercell:
    #
    #     0 <= f1 + shift <= S
    #     0 <= f2 + shift <= S
    #
    S = np.array([Sa, Sb, Sc], dtype=float)

    lower = np.maximum(-f1, -f2)
    upper = np.minimum(S - f1, S - f2)

    shift_min = np.ceil(lower - tol).astype(int)
    shift_max = np.floor(upper + tol).astype(int)

    if np.any(shift_min > shift_max):
        return []

    bond_copies = []
    seen = set()

    for shift in product(
        range(shift_min[0], shift_max[0] + 1),
        range(shift_min[1], shift_max[1] + 1),
        range(shift_min[2], shift_max[2] + 1),
    ):

        shift = np.asarray(shift, dtype=float)
        translation = shift @ cell_vectors

        p1_copy = p1 + translation
        p2_copy = p2 + translation

        # Remove duplicates caused by tolerance
        key = (
            tuple(np.round(p1_copy / tol).astype(int)),
            tuple(np.round(p2_copy / tol).astype(int)),
        )

        if key not in seen:
            seen.add(key)
            bond_copies.append((p1_copy, p2_copy))

    return bond_copies


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



def get_pbc_ts_geometry(
    TS_data,
    grid_shape,
    grid_vectors,
    supercell,
):

    Na, Nb, Nc = grid_shape
    Sa, Sb, Sc = supercell

    a_vec, b_vec, c_vec = map(np.asarray, grid_vectors)

    Na_i = Na - 1
    Nb_i = Nb - 1
    Nc_i = Nc - 1

    def grid2cart(p):
        return (
            p[0]*a_vec +
            p[1]*b_vec +
            p[2]*c_vec
        )

    start = np.asarray(TS_data["start_center"], float)
    end   = np.asarray(TS_data["end_center"], float)

    cross = np.asarray(TS_data["cross_vector"], int)

    # Put the second endpoint in the neighboring periodic image
    end = end + cross*np.array([Na_i, Nb_i, Nc_i])

    mins = np.zeros(3)

    maxs = np.array([
        Sa*Na_i,
        Sb*Nb_i,
        Sc*Nc_i
    ], float)

    points = []
    lines = []
    centers = []

    idx = 0

    # Search neighboring translated copies
    for ia in range(-1, Sa+1):
        for ib in range(-1, Sb+1):
            for ic in range(-1, Sc+1):

                shift = np.array([
                    ia*Na_i,
                    ib*Nb_i,
                    ic*Nc_i
                ], float)

                clipped = clip_to_box(
                    start + shift,
                    end + shift,
                    mins,
                    maxs,
                )

                if clipped is None:
                    continue

                q1, q2 = clipped

                p1 = grid2cart(q1)
                p2 = grid2cart(q2)

                points.extend((p1, p2))
                lines.extend((2, idx, idx+1))
                centers.append(0.5*(p1+p2))

                idx += 2

    return (
        np.asarray(points),
        np.asarray(lines),
        np.asarray(centers),
    )

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

def cartesian_to_fractional(position, origin, grid_vectors, grid_shape):
    """
    Convert Cartesian coordinates to fractional coordinates.

    Parameters
    ----------
    position : array-like, shape (3,)
        Cartesian coordinate of the atom.
    origin : array-like, shape (3,)
        Origin of the grid/unit cell.
    grid_vectors : array-like, shape (3,3)
        Grid vectors (rows or columns depending on your convention).
    grid_shape : array-like, shape (3,)
        Number of grid points along each lattice direction.

    Returns
    -------
    np.ndarray, shape (3,)
        Fractional coordinates.
    """

    position = np.asarray(position, dtype=float)
    origin = np.asarray(origin, dtype=float)
    grid_vectors = np.asarray(grid_vectors, dtype=float)
    grid_shape = np.asarray(grid_shape, dtype=float)

    # Unit cell vectors
    cell_vectors = grid_shape[:, None] * grid_vectors

    # Cartesian -> fractional
    fractional = np.linalg.solve(
        cell_vectors.T,
        position - origin
    )

    return fractional



