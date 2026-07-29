import numpy as np


def get_process_cross_vector(
    start_point,
    end_point,
    C1_C2_TS,
    minIDmatrix_clusters
):
    """
    Find the periodic-boundary cross vector between start_point and end_point.

    Parameters
    ----------
    start_point : array-like, shape (3,)
        Starting point [i, j, k].

    end_point : array-like, shape (3,)
        Target point [i, j, k].

    C1_C2_TS : array-like, shape (3,)
        [C1, C2, TS].

        TS is currently not used, matching the C++/Fortran implementation.

    minIDmatrix_clusters : np.ndarray, shape (nx, ny, nz)
        Cluster ID matrix.

    Returns
    -------
    np.ndarray, shape (3,)
        Crossing vector [cross_i, cross_j, cross_k].

        Returns [100, 100, 100] if the endpoint cannot be reached.
    """

    # ------------------------------------------------------------
    # Convert inputs
    # ------------------------------------------------------------

    start_point = np.asarray(start_point)
    end_point = np.asarray(end_point)
    C1_C2_TS = np.asarray(C1_C2_TS)

    if start_point.shape != (3,):
        raise ValueError(
            "start_point must have shape (3,)"
        )

    if end_point.shape != (3,):
        raise ValueError(
            "end_point must have shape (3,)"
        )

    if C1_C2_TS.shape != (3,):
        raise ValueError(
            "C1_C2_TS must have shape (3,)"
        )

    # ------------------------------------------------------------
    # Grid dimensions
    # ------------------------------------------------------------

    nx, ny, nz = minIDmatrix_clusters.shape

    # ------------------------------------------------------------
    # Cluster IDs
    # ------------------------------------------------------------

    C1 = C1_C2_TS[0]
    C2 = C1_C2_TS[1]

    # ------------------------------------------------------------
    # Start and end coordinates
    # ------------------------------------------------------------

    i_start, j_start, k_start = start_point
    i_end, j_end, k_end = end_point

    # ------------------------------------------------------------
    # Visited array
    # ------------------------------------------------------------

    checked = np.zeros(
        (nx, ny, nz),
        dtype=np.bool_
    )

    # ------------------------------------------------------------
    # BFS queue
    #
    # Each point:
    #
    # [i, j, k,
    #  cross_i, cross_j, cross_k]
    #
    # ------------------------------------------------------------

    point_list = [
        [
            i_start,
            j_start,
            k_start,
            0,
            0,
            0
        ]
    ]

    checked[
        i_start,
        j_start,
        k_start
    ] = True

    point_index = 0

    # ------------------------------------------------------------
    # Breadth-first search
    # ------------------------------------------------------------

    while point_index < len(point_list):

        (
            i,
            j,
            k,
            current_cross_i,
            current_cross_j,
            current_cross_k
        ) = point_list[point_index]

        # --------------------------------------------------------
        # Periodic neighbours
        # --------------------------------------------------------

        ip = (i + 1) % nx
        im = (i - 1) % nx

        jp = (j + 1) % ny
        jm = (j - 1) % ny

        kp = (k + 1) % nz
        km = (k - 1) % nz

        # --------------------------------------------------------
        # Crossing information
        # --------------------------------------------------------

        cross_x_positive = (
            1 if i == nx - 1 else 0
        )

        cross_x_negative = (
            -1 if i == 0 else 0
        )

        cross_y_positive = (
            1 if j == ny - 1 else 0
        )

        cross_y_negative = (
            -1 if j == 0 else 0
        )

        cross_z_positive = (
            1 if k == nz - 1 else 0
        )

        cross_z_negative = (
            -1 if k == 0 else 0
        )

        # --------------------------------------------------------
        # Six neighbours
        # --------------------------------------------------------

        neighbours = (

            (
                ip,
                j,
                k,
                current_cross_i + cross_x_positive,
                current_cross_j,
                current_cross_k
            ),

            (
                im,
                j,
                k,
                current_cross_i + cross_x_negative,
                current_cross_j,
                current_cross_k
            ),

            (
                i,
                jp,
                k,
                current_cross_i,
                current_cross_j + cross_y_positive,
                current_cross_k
            ),

            (
                i,
                jm,
                k,
                current_cross_i,
                current_cross_j + cross_y_negative,
                current_cross_k
            ),

            (
                i,
                j,
                kp,
                current_cross_i,
                current_cross_j,
                current_cross_k + cross_z_positive
            ),

            (
                i,
                j,
                km,
                current_cross_i,
                current_cross_j,
                current_cross_k + cross_z_negative
            )
        )

        # --------------------------------------------------------
        # Examine neighbours
        # --------------------------------------------------------

        for (
            ni,
            nj,
            nk,
            neighbour_cross_i,
            neighbour_cross_j,
            neighbour_cross_k
        ) in neighbours:

            neighbour_id = (
                minIDmatrix_clusters[
                    ni,
                    nj,
                    nk
                ]
            )

            # Only C1 and C2 are included
            if (
                neighbour_id != C1
                and neighbour_id != C2
            ):
                continue

            # Already visited?
            if checked[
                ni,
                nj,
                nk
            ]:
                continue

            # Mark visited
            checked[
                ni,
                nj,
                nk
            ] = True

            # Add to BFS queue
            point_list.append(
                [
                    ni,
                    nj,
                    nk,
                    neighbour_cross_i,
                    neighbour_cross_j,
                    neighbour_cross_k
                ]
            )

            # Endpoint reached
            if (
                ni == i_end
                and nj == j_end
                and nk == k_end
            ):

                return np.array(
                    [
                        neighbour_cross_i,
                        neighbour_cross_j,
                        neighbour_cross_k
                    ],
                    dtype=np.int64
                )

        point_index += 1

    # ------------------------------------------------------------
    # Endpoint not found
    # ------------------------------------------------------------

    return np.array(
        [100, 100, 100],
        dtype=np.int64
    )