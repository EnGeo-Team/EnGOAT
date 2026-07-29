import matplotlib.pyplot as plt


def plot_energy_histogram(histogram, title=None, color="steelblue"):

    if not histogram:
        print("No histogram data found.")
        return

    energies = [float(energy) for energy in histogram.keys()]
    volumes = list(histogram.values())

    if len(energies) == 1:
        width = 1.0
    else:
        width = min(b - a for a, b in zip(energies, energies[1:]))

    fig, ax = plt.subplots(figsize=(7, 4))

    ax.bar(
        energies,
        volumes,
        width=width,
        align="edge",
        color=color,
        edgecolor="black",
        linewidth=0.7,
        alpha=0.85,
    )

    ax.set_xlabel("Energy [kJ/mol]")
    ax.set_ylabel(r"Volume [$\AA^3$]")

    if title is not None:
        ax.set_title(title)

    ax.grid(axis="y", linestyle="--", alpha=0.3)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    plt.show(block=False)