import numpy as np
import data
import logging
#import initiate_cluster                                                                                                         #Fortran file containing the subroutine for finding new clusters at a current level
logger = logging.getLogger(__name__)

def find_new_clusters(level):

    # Find all unexplored points
    unexplored_x, unexplored_y, unexplored_z = np.where(
        data.Cluster_matrix.Levels == 0
    )

    # Grid dimensions
    nx, ny, nz = data.level_matrix.shape

    # ------------------------------------------------------------
    # Loop over all unexplored points
    # ------------------------------------------------------------

    for unexplored_point in range(len(unexplored_x)):

        i_0 = unexplored_x[unexplored_point]
        j_0 = unexplored_y[unexplored_point]
        k_0 = unexplored_z[unexplored_point]

        # Only points belonging to the current level can start a cluster
        if data.level_matrix[i_0, j_0, k_0] != level:
            continue

        # --------------------------------------------------------
        # Periodic boundary conditions for the starting point
        # --------------------------------------------------------

        i = i_0
        j = j_0
        k = k_0

        ip = (i + 1) % nx
        im = (i - 1) % nx

        jp = (j + 1) % ny
        jm = (j - 1) % ny

        kp = (k + 1) % nz
        km = (k - 1) % nz

        # Crossing information:
        #
        # [positive x, negative x,
        #  positive y, negative y,
        #  positive z, negative z]

        cross = np.array([
            1 if i == nx - 1 else 0,
            -1 if i == 0 else 0,
            1 if j == ny - 1 else 0,
            -1 if j == 0 else 0,
            1 if k == nz - 1 else 0,
            -1 if k == 0 else 0
        ], dtype=np.int64)

        # IDs of neighbouring points
        neighbour_ID = [
            data.Cluster_matrix.IDs[ip, j, k],
            data.Cluster_matrix.IDs[im, j, k],
            data.Cluster_matrix.IDs[i, jp, k],
            data.Cluster_matrix.IDs[i, jm, k],
            data.Cluster_matrix.IDs[i, j, kp],
            data.Cluster_matrix.IDs[i, j, km]
        ]

        # Same condition as the Fortran routine:
        # only start a new cluster if all six neighbours have ID = 0
        if sum(neighbour_ID) != 0:
            continue

        # --------------------------------------------------------
        # Create a new cluster
        # --------------------------------------------------------

        data.Cluster.N_clusters[0] += 1

        cluster_id = data.Cluster.N_clusters[0]

        data.Cluster_matrix.Levels[i, j, k] = level
        data.Cluster_matrix.IDs[i, j, k] = cluster_id

        # List of points in this cluster
        #
        # Columns:
        # 0: x
        # 1: y
        # 2: z
        # 3: crossing x
        # 4: crossing y
        # 5: crossing z
        # 6: boundary flag
        # 7: level
        # 8: transition-state information

        cluster_points = [
            [i, j, k, 0, 0, 0, 0, level, 0]
        ]

        # --------------------------------------------------------
        # Breadth-first search through the cluster
        # --------------------------------------------------------

        point_index = 0

        while point_index < len(cluster_points):

            # Current point
            i, j, k = cluster_points[point_index][0:3]

            # Periodic neighbours
            ip = (i + 1) % nx
            im = (i - 1) % nx

            jp = (j + 1) % ny
            jm = (j - 1) % ny

            kp = (k + 1) % nz
            km = (k - 1) % nz

            # Crossing information for current point
            cross = np.array([
                1 if i == nx - 1 else 0,
                -1 if i == 0 else 0,
                1 if j == ny - 1 else 0,
                -1 if j == 0 else 0,
                1 if k == nz - 1 else 0,
                -1 if k == 0 else 0
            ], dtype=np.int64)

            boundary = 0

            # The six neighbours:
            #
            # (neighbour i, neighbour j, neighbour k,
            #  crossing contribution x,
            #  crossing contribution y,
            #  crossing contribution z)

            neighbours = [
                (ip, j,  k,  cross[0], 0,        0),
                (im, j,  k,  cross[1], 0,        0),
                (i,  jp, k,  0,        cross[2], 0),
                (i,  jm, k,  0,        cross[3], 0),
                (i,  j,  kp, 0,        0,        cross[4]),
                (i,  j,  km, 0,        0,        cross[5])
            ]

            # ----------------------------------------------------
            # Examine all six neighbours
            # ----------------------------------------------------

            for ni, nj, nk, cross_x, cross_y, cross_z in neighbours:

                # Only examine points that have not yet been explored
                if data.Cluster_matrix.Levels[ni, nj, nk] != 0:
                    continue

                # If it belongs to the current level, add it
                if data.level_matrix[ni, nj, nk] == level:

                    data.Cluster_matrix.Levels[ni, nj, nk] = level
                    data.Cluster_matrix.IDs[ni, nj, nk] = cluster_id

                    # Calculate crossing information relative to
                    # the current point
                    cross_i = (
                        data.cross_matrix.i[i, j, k]
                        + cross_x
                    )

                    cross_j = (
                        data.cross_matrix.j[i, j, k]
                        + cross_y
                    )

                    cross_k = (
                        data.cross_matrix.k[i, j, k]
                        + cross_z
                    )

                    # Store crossing information for the new point
                    data.cross_matrix.i[ni, nj, nk] = cross_i
                    data.cross_matrix.j[ni, nj, nk] = cross_j
                    data.cross_matrix.k[ni, nj, nk] = cross_k

                    # Add point to cluster
                    cluster_points.append([
                        ni,
                        nj,
                        nk,
                        cross_i,
                        cross_j,
                        cross_k,
                        0,
                        level,
                        0
                    ])

                else:

                    # The current point has a neighbour that is not
                    # part of this cluster
                    boundary = 1

            # Store boundary information for current point
            cluster_points[point_index][6] = boundary

            # Move to the next point in the cluster
            point_index += 1

        # --------------------------------------------------------
        # Convert cluster points to NumPy array
        # --------------------------------------------------------

        Cluster_points = np.asarray(
            cluster_points,
            dtype=np.int64
        )

        # --------------------------------------------------------
        # Get energies directly from Energy_matrix
        # --------------------------------------------------------

        point_energies = data.Energy_matrix[
            Cluster_points[:, 0],
            Cluster_points[:, 1],
            Cluster_points[:, 2]
        ]

        # Minimum energy
        E_min = np.min(point_energies)

        # Point with minimum energy
        Emin_point = np.argmin(point_energies)

        center = (
            Cluster_points[Emin_point, 0],
            Cluster_points[Emin_point, 1],
            Cluster_points[Emin_point, 2]
        )

        # --------------------------------------------------------
        # Add cluster to data structure
        # --------------------------------------------------------

        data.Cluster.Cluster_list.append(
            data.Cluster(
                cluster_id,
                Cluster_points,
                E_min,
                center,
                point_energies
            )
        )

        data.Cluster.Cluster_families.append(
            {cluster_id}
        )

    return None