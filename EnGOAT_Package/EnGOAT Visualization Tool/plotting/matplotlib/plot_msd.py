import numpy as np
import matplotlib.pyplot as plt


def plot_tunnel_msd(
    kMC_data,
    tunnel_system_id,
    directions,
    direction_colors=None,
    conditions = None
):
    """
    Plot MSD curves for multiple directions of one tunnel system.

    Individual runs are shown transparently, with the mean MSD shown
    as a thicker line and ±1 standard deviation shown as a shaded region.
    """

    tunnel_system_id = str(tunnel_system_id)

    tunnel_system = kMC_data["tunnel_systems"][tunnel_system_id]

    fig, ax = plt.subplots(figsize=(8, 6))

    for direction in directions:

        direction_data = tunnel_system["directions"][direction]

        color = None
        if direction_colors is not None:
            color = direction_colors.get(direction)

        runs = direction_data["MSD"]

        # Plot individual runs
        for run_data in runs.values():
            time = np.asarray(run_data["time"])
            msd = np.asarray(run_data["MSD"])

            ax.plot(
                time,
                msd,
                color=color,
                alpha=0.15,
                linewidth=0.8,
            )

        # Calculate mean and standard deviation across runs
        msd_values = np.array([
            run_data["MSD"]
            for run_data in runs.values()
        ])

        time = np.asarray(next(iter(runs.values()))["time"])

        mean_msd = np.mean(msd_values, axis=0)
        std_msd = np.std(msd_values, axis=0)

        # Mean curve
        ax.plot(
            time,
            mean_msd,
            color=color,
            linewidth=2.5,
            label=direction,
        )

        # Standard deviation region
        ax.fill_between(
            time,
            mean_msd - std_msd,
            mean_msd + std_msd,
            color=color,
            alpha=0.15,
        )

        # Plot fitting range thresholds
        if conditions is not None:
            key = "iso" if direction == "3D" else direction

            if key in conditions:
                start = conditions[key]["start"]
                end = conditions[key]["end"]

                # Highlight fitting region
                ax.axhspan(
                    start,
                    end,
                    color=color,
                    alpha=0.08,
                    zorder=0,
                )

                # Faint boundary lines
                ax.axhline(
                    start,
                    color=color,
                    linestyle="--",
                    linewidth=1.0,
                    alpha=0.35,
                )

                ax.axhline(
                    end,
                    color=color,
                    linestyle="--",
                    linewidth=1.0,
                    alpha=0.35,
                )

    # Logarithmic axes
    ax.set_xscale("log")
    ax.set_yscale("log")

    # Labels
    ax.set_xlabel("Time (s)", fontsize=12)
    ax.set_ylabel(r"MSD ($\mathrm{\AA^2}$)", fontsize=12)

    ax.set_title(
        f"Mean Squared Displacement — Tunnel System {tunnel_system_id}",
        fontsize=14,
        pad=12,
    )

    # Grid
    ax.grid(
        True,
        which="major",
        linestyle="--",
        linewidth=0.7,
        alpha=0.5,
    )

    ax.grid(
        True,
        which="minor",
        linestyle=":",
        linewidth=0.5,
        alpha=0.25,
    )

    # Legend
    ax.legend(
        title="Direction",
        frameon=True,
        fontsize=10,
        title_fontsize=10,
    )

    # Remove unnecessary top/right spines
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()

    plt.show(block=False)

    return fig, ax