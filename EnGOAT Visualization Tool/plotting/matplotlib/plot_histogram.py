import numpy as np
import matplotlib.pyplot as plt


def plot_energy_histogram(level_matrix, mask, grid, E_step, title=None, color = "steelblue"):

    a_vec = np.array(grid["a_vector"])
    b_vec = np.array(grid["b_vector"])
    c_vec = np.array(grid["c_vector"])

    V_vox = abs(np.dot(a_vec, np.cross(b_vec, c_vec)))

    levels = level_matrix[mask]

    if len(levels) == 0:
        print("No points found.")
        return

    counts = np.bincount(levels)
    volumes = counts * V_vox

    # Shift bins one interval to the left
    edges = (np.arange(len(counts) + 1) - 1) * E_step

    fig, ax = plt.subplots(figsize=(7, 4))

    ax.bar(
        edges[:-1],
        volumes,
        width=E_step,
        align="edge",
        color=color,
        edgecolor="black",
        linewidth=0.7,
        alpha=0.85
    )

    ax.set_xlabel("Energy [kJ/mol]")
    ax.set_ylabel(r"Volume [$\AA^3$]")

    if title is not None:
        ax.set_title(title)

    ax.grid(
        axis="y",
        linestyle="--",
        alpha=0.3
    )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    plt.show(block=False)