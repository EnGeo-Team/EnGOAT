import numpy as np
import matplotlib.pyplot as plt
import os
import re
from collections import defaultdict

os.makedirs("Min_E_path_plots", exist_ok=True)

def plot_min_E_path(tunnel, direction, basin_IDs, colors):
    folder = "NumPy_matrices"
    basin_filename = "Basin_matrix.npy"
    basin_filepath = os.path.join(folder, basin_filename)
    Cluster_matrix = np.load(basin_filepath)
    TS_filename = "TS_matrix.npy"
    TS_filepath = os.path.join(folder, TS_filename)
    TS_matrix = np.load(TS_filepath)
    Level_filename = "Level_matrix.npy"
    Level_filepath = os.path.join(folder, Level_filename)
    Level_matrix = np.load(Level_filepath)
    Tunnel_filename = "Tunnel_matrix.npy"
    Tunnel_filepath = os.path.join(folder, Tunnel_filename)
    Tunnel_matrix = np.load(Tunnel_filepath)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.grid(False)

    denom = min(Level_matrix.shape)
    ax.set_box_aspect([Level_matrix.shape[0]/denom,Level_matrix.shape[1]/denom,Level_matrix.shape[2]/denom])
    
    ax.set_xlabel("a")
    ax.set_ylabel("b")
    ax.set_zlabel("c")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])

    ax.set_title(f"Minimum E Path in {direction} Direction")
    

    mask = Tunnel_matrix == tunnel
    ax.voxels(mask, facecolors="skyblue", edgecolor='none', alpha=0.3)

    for i, basin_ID in enumerate(basin_IDs):
        mask = Cluster_matrix == basin_ID
        ax.voxels(mask, facecolors=colors[i], edgecolor='none', alpha=0.8)

    plt.savefig(os.path.join("Min_E_path_plots", f"Min_E_Path_{direction}_tunnel{tunnel}_Visualized.png"), dpi = 300)
    #plt.show()
    plt.close()
    return None

def bezier_quad(P0, P1, P2, t):
    """Quadratic Bézier curve from P0 -> P1 -> P2"""
    return (1 - t)**2 * P0 + 2 * (1 - t) * t * P1 + t**2 * P2    

def plot_min_E_pathway(tunnel, direction):
    #PLOT MINIMUM E PATHWAY
    filename = f"min_E_path_{direction}_tunnel{tunnel}.dat"
    folder = "Tunnel_data"
    filepath = os.path.join(folder, filename)
    data = np.loadtxt(filepath, skiprows=2)
    if data.ndim == 1:
        basin_ids = [data[0].astype(int)]
        end_ids   = [data[1].astype(int)]
        ts_ids    = [data[2].astype(int)]
        basin_E   = [data[3]]
        end_E     = [data[4]]
        TS_E        = [data[5]]
        basin_coords = [data[6]]
        end_coords = [data[7]]
        PBC_crossing = [data[8]]
    else:
        basin_ids = data[:, 0].astype(int)
        end_ids   = data[:, 1].astype(int)
        ts_ids    = data[:, 2].astype(int)
        basin_E   = data[:, 3]
        end_E     = data[:, 4]
        TS_E        = data[:, 5]
        basin_coords = data[:, 6]
        end_coords = data[:, 7]
        PBC_crossing = data[:, 8]

    plt.figure(figsize=(15, 5))

    n_basins = len(basin_ids)
    basin_width = 1
    for i in range(len(basin_coords)-1):
        if abs(basin_coords[i]-basin_coords[i+1]) != 0:
            basin_width = min(basin_width, abs(basin_coords[i]-basin_coords[i+1]))
    basin_width = 0.25/n_basins

    y_max = np.max(TS_E)*1.1
    y_min = min(basin_E) - np.max(TS_E)*0.1
    y_range = y_max-y_min

    # Plot basin plateaus and labels
    colors = plt.cm.tab20(np.linspace(0, 1, len(basin_ids)))
    for i, (bid, E, x) in enumerate(zip(basin_ids, basin_E, basin_coords)):
        color = colors[i]
        plt.hlines(E, x - basin_width, x + basin_width, linewidth=3, color=color)
        if x - basin_width < 0:
            plt.hlines(E, x - basin_width + 1, x + basin_width + 1, linewidth=3, color=color)
            plt.text(x + basin_width, E - 0.05*y_range, f"E = {E} kJ/mol", ha='center', va='bottom', fontsize=9, color = color)
            plt.text(x + basin_width/2, E + 0.015*y_range, f"B{bid}", ha='center', va='bottom', fontsize=9)
        elif x + basin_width > 1:
            plt.hlines(E, x - basin_width - 1, x + basin_width - 1, linewidth=3, color=color)
            plt.text(x - basin_width, E - 0.05*y_range, f"E = {E} kJ/mol", ha='center', va='bottom', fontsize=9, color = color)
            plt.text(x - basin_width/2, E + 0.015*y_range, f"B{bid}", ha='center', va='bottom', fontsize=9)
        else:
            plt.text(x, E - 0.05*y_range, f"E = {E} kJ/mol", ha='center', va='bottom', fontsize=9, color = color)
            plt.text(x, E + 0.015*y_range, f"B{bid}", ha='center', va='bottom', fontsize=9)

    # Plot Bezier Curves
    for i in range(len(data)):
        j = i+1
        if i == len(data)-1:
            j = 0
        if data.ndim == 1:
            i = 0
            j = 0
        x0 = basin_coords[i] + basin_width
        x2 = basin_coords[j] - basin_width

        y0 = basin_E[i]
        y2 = basin_E[j]
        y_ts = TS_E[i]
        y1 = 2 * y_ts - 0.5 * y0 - 0.5 * y2

        if PBC_crossing[i] == 0:
            xm = 0.5 * (x0 + x2)
            if x0 < x2:
                P0 = np.array([x0, y0])
                P1 = np.array([xm, y1])  # transition state peak
                P2 = np.array([x2, y2])
                t = np.linspace(0, 1, 100)
                curve = np.array([bezier_quad(P0, P1, P2, ti) for ti in t])
                plt.plot(curve[:, 0], curve[:, 1], color='black')
                
                if x2 > 1:
                    P0_ = np.array([x0-1, y0])
                    P1_ = np.array([xm-1, y1])  # transition state peak
                    P2_ = np.array([x2-1, y2])
                    t = np.linspace(0, 1, 100)
                    curve = np.array([bezier_quad(P0_, P1_, P2_, ti) for ti in t])
                    plt.plot(curve[:, 0], curve[:, 1], color='black')

                plt.text(
                xm,
                y_ts + y_range*0.01,  # small offset above peak
                f"$E$ = {y_ts:.2f} kJ/mol",
                ha='center',
                va='bottom',
                fontsize=8
                )
            elif x0 > x2:
                P0 = np.array([x0, y0])
                P1 = np.array([xm, y1])  # transition state peak
                P2 = np.array([x2, y2])
                t = np.linspace(0, 1, 100)
                curve = np.array([bezier_quad(P2, P1, P0, ti) for ti in t])
                plt.plot(curve[:, 0], curve[:, 1], color='black')

                if x2 < 0:
                    P0_ = np.array([x0+1, y0])
                    P1_ = np.array([xm+1, y1])  # transition state peak
                    P2_ = np.array([x2+1, y2])
                    t = np.linspace(0, 1, 100)
                    curve = np.array([bezier_quad(P0_, P1_, P2_, ti) for ti in t])
                    plt.plot(curve[:, 0], curve[:, 1], color='black')

                plt.text(
                xm,
                y_ts + y_range*0.01,  # small offset above peak
                f"$E$ = {y_ts:.2f} kJ/mol",
                ha='center',
                va='bottom',
                fontsize=8
                )
            elif x0 == x2:
                P0 = np.array([x0, y0])
                P1 = np.array([xm, y1])  # transition state peak
                P2 = np.array([x2, y2])
                plt.plot([x0, x0], [y0, y_ts], color='black')
                plt.plot([x0, x0], [y2, y_ts], color='black')

                plt.text(
                x0,
                y_ts + y_range*0.01,  # small offset above peak
                f"$E$ = {y_ts:.2f} kJ/mol",
                ha='center',
                va='bottom',
                fontsize=8
                )
            
        elif PBC_crossing[i] == 1:
            x2 += 1
            xm = 0.5 * (x0 + x2)  # midpoint for label (not needed for quadratic)
            P0 = np.array([x0, y0])
            P1 = np.array([xm, y1])  # transition state peak
            P2 = np.array([x2, y2])

            t = np.linspace(0, 1, 100)
            curve = np.array([bezier_quad(P0, P1, P2, ti) for ti in t])
            plt.plot(curve[:, 0], curve[:, 1], color='black')

            P0[0] -= 1
            P1[0] -= 1
            P2[0] -= 1
            curve = np.array([bezier_quad(P0, P1, P2, ti) for ti in t])
            plt.plot(curve[:, 0], curve[:, 1], color='black')

            if xm < 1:
                plt.text(
                    xm,
                    y_ts + y_range*0.01,  # small offset above peak
                    f"$E$ = {y_ts:.2f} kJ/mol",
                    ha='center',
                    va='bottom',
                    fontsize=8
                )
            else:
                plt.text(
                    xm - 1,
                    y_ts + y_range*0.01,  # small offset above peak
                    f"$E$ = {y_ts:.2f} kJ/mol",
                    ha='center',
                    va='bottom',
                    fontsize=8
                )
        
        elif PBC_crossing[i] == -1:
            x2 -= 1
            xm = 0.5 * (x0 + x2)  # midpoint for label (not needed for quadratic)
            P0 = np.array([x0, y0])
            P1 = np.array([xm, y1])  # transition state peak
            P2 = np.array([x2, y2])

            t = np.linspace(0, 1, 100)
            curve = np.array([bezier_quad(P2, P1, P0, ti) for ti in t])
            plt.plot(curve[:, 0], curve[:, 1], color='black')

            P0[0] += 1
            P1[0] += 1
            P2[0] += 1
            curve = np.array([bezier_quad(P2, P1, P0, ti) for ti in t])
            plt.plot(curve[:, 0], curve[:, 1], color='black')

            if xm > 0:
                plt.text(
                    xm,
                    y_ts + y_range*0.01,  # small offset above peak
                    f"$E$ = {y_ts:.2f} kJ/mol",
                    ha='center',
                    va='bottom',
                    fontsize=8
                )
            else:
                plt.text(
                    xm + 1,
                    y_ts + y_range*0.01,  # small offset above peak
                    f"$E$ = {y_ts:.2f} kJ/mol",
                    ha='center',
                    va='bottom',
                    fontsize=8
                )

    title = f"Minimum Energy Path in {direction} direction"
    plt.xlabel(f"Fractional Coordinate {direction}")
    plt.ylabel("Energy [kJ/mol]")
    plt.title(title)

    plt.xlim(0, 1)
    plt.ylim(y_min, y_max )

    plt.xticks([0, 0.2, 0.4, 0.6, 0.8, 1])  # remove xticks

    # Ensure all spines are visible (replace plt.box(True))
    for spine in plt.gca().spines.values():
        spine.set_visible(True)

    plt.tight_layout()
    plt.savefig(os.path.join("Min_E_path_plots", f"Min_E_Path_{direction}_tunnel{tunnel}.png"), dpi = 300)
    plt.close()

    plot_min_E_path(tunnel, direction, basin_ids, colors)

    return None


folder = "Tunnel_data"

# pattern: min_E_path_a_tunnel1
pattern = re.compile(r"min_E_path_([abc])_tunnel(\d+)")

paths = defaultdict(list)

for filename in os.listdir(folder):
    match = pattern.match(filename)
    if match:
        direction = match.group(1)
        tunnel = int(match.group(2))
        paths[tunnel].append(direction)

# Print results
for tunnel in sorted(paths):
    print(f"Tunnel {tunnel}: directions {sorted(paths[tunnel])}")
    for direction in sorted(paths[tunnel]):
        print(f"Plotting direction {direction}...")
        plot_min_E_pathway(tunnel, direction)
        
