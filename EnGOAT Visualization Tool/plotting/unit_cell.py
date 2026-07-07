import numpy as np

def build_unit_cell(grid):

    Na, Nb, Nc = grid["grid_points"]

    a = np.array(grid["a_vector"])
    b = np.array(grid["b_vector"])
    c = np.array(grid["c_vector"])


    A = Na * a
    B = Nb * b
    C = Nc * c

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

    return corners, edges

import numpy as np


def get_unit_cell_params(grid):

    Na, Nb, Nc = grid["grid_points"]

    a = np.array(grid["a_vector"])
    b = np.array(grid["b_vector"])
    c = np.array(grid["c_vector"])

    A = Na * np.linalg.norm(a)
    B = Nb * np.linalg.norm(b)
    C = Nc * np.linalg.norm(c)

    alpha = np.degrees(
        np.arccos(
            np.clip(
                np.dot(b, c)
                / (np.linalg.norm(b) * np.linalg.norm(c)),
                -1.0,
                1.0
            )
        )
    )

    beta = np.degrees(
        np.arccos(
            np.clip(
                np.dot(a, c)
                / (np.linalg.norm(a) * np.linalg.norm(c)),
                -1.0,
                1.0
            )
        )
    )

    gamma = np.degrees(
        np.arccos(
            np.clip(
                np.dot(a, b)
                / (np.linalg.norm(a) * np.linalg.norm(b)),
                -1.0,
                1.0
            )
        )
    )

    return A, B, C, alpha, beta, gamma