import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import read_basin_data
import os

# -----------------------------
# Convert grid → Cartesian
# -----------------------------
def grid_to_cartesian(idx, grid):
    i, j, k = idx
    a_vec = np.array(grid["a"])
    b_vec = np.array(grid["b"])
    c_vec = np.array(grid["c"])
    return i * a_vec + j * b_vec + k * c_vec


def plot_graph_network(ts_map, ts_colors, iso_map, iso_colors):

    # -----------------------------
    # Build unit cell (parallelepiped)
    # -----------------------------

    grid = read_basin_data.read_cell_data("output.log")

    A = grid["grid_size"][0]*np.array(grid["a"])
    B = grid["grid_size"][1]*np.array(grid["b"])
    C = grid["grid_size"][2]*np.array(grid["c"])

    origin = np.array([0, 0, 0])

    vertices = np.array([
        origin,
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

    # Plot

    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111, projection='3d')

    for i, j in edges:
        ax.plot(
            [vertices[i][0], vertices[j][0]],
            [vertices[i][1], vertices[j][1]],
            [vertices[i][2], vertices[j][2]],
            color='black',
            linewidth=1
        )


    # -----------------------------
    # Plot the basins
    # -----------------------------

    tunnel_systems = list(ts_map.keys())
    #colors = ["blue", "red"]#plt.cm.tab20(np.linspace(0, 1, len(tunnel_systems)))
    for ts in tunnel_systems:
        color = ts_colors[ts]
        basins = ts_map[ts]
        PBC_points = []
        for b in basins:
            
            shifts = []

            if b.center[0] == 0:
                shifts.append([grid["grid_size"][0], 0, 0])
            elif b.center[0] == grid["grid_size"][0]:
                shifts.append([-grid["grid_size"][0], 0, 0])

            if b.center[1] == 0:
                shifts.append([0, grid["grid_size"][1], 0])
            elif b.center[1] == grid["grid_size"][1]:
                shifts.append([0, -grid["grid_size"][1], 0])

            if b.center[2] == 0:
                shifts.append([0, 0, grid["grid_size"][2]])
            elif b.center[2] == grid["grid_size"][2]:
                shifts.append([0, 0, -grid["grid_size"][2]])

            # always include original
            shifts = [[0, 0, 0]] + shifts

            shifts2 = []
            for s1_idx in range(len(shifts)):
                s1 = np.array(shifts[s1_idx])
                shifts2.append(s1)
                for s2_idx in range(s1_idx+1, len(shifts)):
                    s2 = np.array(shifts[s2_idx])
                    shifts2.append(s1 + s2)


            for shift in shifts2:
                point = np.array(b.center) + np.array(shift)
                PBC_points.append(grid_to_cartesian(point, grid))

        PBC_points = np.array(PBC_points)

        if len(PBC_points)!=0:
            sc = ax.scatter(
            PBC_points[:, 0],
            PBC_points[:, 1],
            PBC_points[:, 2],
            c=color,
            s=80,          # size of spheres
            edgecolors='k' # optional: black outline
            )

        # -----------------------------
        # Plot the transitions
        # -----------------------------

        processes = read_basin_data.read_processes(f"TuTraSt_data/tunnel{ts}_data.dat")

        for process in processes:
            start_point = np.array(process["start_point"])
            end_point = np.array(process["end_point"]) + np.array((process["process_vector"][0]*(grid["grid_size"][0]),
                                                       process["process_vector"][1]*(grid["grid_size"][1]),
                                                       process["process_vector"][2]*(grid["grid_size"][2])))
            if (
            (end_point >= 0).all() and
            (end_point <= np.array(grid["grid_size"])).all()):
                start_point = grid_to_cartesian(start_point, grid)
                end_point = grid_to_cartesian(end_point, grid)
                ax.plot(
                [start_point[0], end_point[0]],
                [start_point[1], end_point[1]],
                [start_point[2], end_point[2]],
                color=color,
                linewidth=1
                )

            start_point = np.array(process["start_point"])- np.array((process["process_vector"][0]*(grid["grid_size"][0]),
                                                       process["process_vector"][1]*(grid["grid_size"][1]),
                                                       process["process_vector"][2]*(grid["grid_size"][2])))
            end_point = np.array(process["end_point"]) 
            if (
            (start_point >= 0).all() and
            (start_point <= np.array(grid["grid_size"])).all()):
                start_point = grid_to_cartesian(start_point, grid)
                end_point = grid_to_cartesian(end_point, grid)
                ax.plot(
                [start_point[0], end_point[0]],
                [start_point[1], end_point[1]],
                [start_point[2], end_point[2]],
                color=color,
                linewidth=1
                )

    processes = read_basin_data.read_processes(f"TuTraSt_data/isolated_processes.dat")
    iso_groups = list(iso_map.keys())
    #colors = ["blue", "red"]#plt.cm.tab20(np.linspace(0, 1, len(tunnel_systems)))
    for iso in iso_groups:
        color = iso_colors[iso]
        basins = iso_map[iso]
        PBC_points = []
        for b in basins:
            

            shifts = []
            corner = [False, False, False]
            if b.center[0] == 0:
                shifts.append([grid["grid_size"][0], 0, 0])
                corner[0] = True
            elif b.center[0] == grid["grid_size"][0]:
                shifts.append([-grid["grid_size"][0], 0, 0])
                corner[0] = True

            if b.center[1] == 0:
                shifts.append([0, grid["grid_size"][1], 0])
                corner[1] = True
            elif b.center[1] == grid["grid_size"][1]:
                shifts.append([0, -grid["grid_size"][1], 0])
                corner[1] = True

            if b.center[2] == 0:
                shifts.append([0, 0, grid["grid_size"][2]])
                corner[2] = True
            elif b.center[2] == grid["grid_size"][2]:
                shifts.append([0, 0, -grid["grid_size"][2]])
                corner[2] = True

            # always include original
            shifts = [[0, 0, 0]] + shifts

            shifts2 = []
            for s1_idx in range(len(shifts)):
                s1 = np.array(shifts[s1_idx])
                shifts2.append(s1)
                for s2_idx in range(s1_idx+1, len(shifts)):
                    s2 = np.array(shifts[s2_idx])
                    for s3_idx in range(s2_idx+1, len(shifts)):
                        s3 = np.array(shifts[s3_idx])
                        shifts2.append(s1 + s2 + s3)
            

                


            for shift in shifts2:
                point = np.array(b.center) + np.array(shift)
                PBC_points.append(grid_to_cartesian(point, grid))

        PBC_points = np.array(PBC_points)

        if len(PBC_points)!=0:
            sc = ax.scatter(
            PBC_points[:, 0],
            PBC_points[:, 1],
            PBC_points[:, 2],
            c=color,
            s=80,          # size of spheres
            edgecolors='k' # optional: black outline
            )

        for process in processes:
            if any(process["B1"] == basin.ID for basin in basins):
                
                start_point = np.array(process["start_point"])
                end_point = np.array(process["end_point"]) + np.array((process["process_vector"][0]*(grid["grid_size"][0]),
                                                           process["process_vector"][1]*(grid["grid_size"][1]),
                                                           process["process_vector"][2]*(grid["grid_size"][2])))
                if (
                (end_point >= 0).all() and
                (end_point <= np.array(grid["grid_size"])).all()):
                    start_point = grid_to_cartesian(start_point, grid)
                    end_point = grid_to_cartesian(end_point, grid)
                    ax.plot(
                    [start_point[0], end_point[0]],
                    [start_point[1], end_point[1]],
                    [start_point[2], end_point[2]],
                    color=color,
                    linewidth=1
                    )

                start_point = np.array(process["start_point"])- np.array((process["process_vector"][0]*(grid["grid_size"][0]),
                                                           process["process_vector"][1]*(grid["grid_size"][1]),
                                                           process["process_vector"][2]*(grid["grid_size"][2])))
                end_point = np.array(process["end_point"]) 
                if (
                (start_point >= 0).all() and
                (start_point <= np.array(grid["grid_size"])).all()):
                    start_point = grid_to_cartesian(start_point, grid)
                    end_point = grid_to_cartesian(end_point, grid)
                    ax.plot(
                    [start_point[0], end_point[0]],
                    [start_point[1], end_point[1]],
                    [start_point[2], end_point[2]],
                    color=color,
                    linewidth=1
                    )

    # Labels
    ax.set_xlabel("X (Å)")
    ax.set_ylabel("Y (Å)")
    ax.set_zlabel("Z (Å)")

    # Set symmetric limits
    xmin, ymin, zmin = vertices.min(axis=0)
    xmax, ymax, zmax = vertices.max(axis=0)
    xmid = (xmin + xmax) / 2
    ymid = (ymin + ymax) / 2
    zmid = (zmin + zmax) / 2
    max_range = max(xmax - xmin, ymax - ymin, zmax - zmin) / 2
    ax.set_xlim(xmid - max_range, xmid + max_range)
    ax.set_ylim(ymid - max_range, ymid + max_range)
    ax.set_zlim(zmid - max_range, zmid + max_range)
    # Force equal aspect
    ax.set_box_aspect([1, 1, 1])

    plt.title("3D Basin Visualization")
    plt.tight_layout()
    plt.savefig(os.path.join("Graph_network_plots", "graph_network.png"), dpi = 300)
    plt.show()
