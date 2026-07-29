import numpy as np
import matplotlib.pyplot as plt


def bezier_quad(P0, P1, P2, t):
    return (
        (1 - t) ** 2 * P0
        + 2 * (1 - t) * t * P1
        + t ** 2 * P2
    )


def plot_MEP_energy_diagram(
    MEP,
    direction,
    plotting_data,
    basins,
    TS_data,
    grid_shape,
):
    """
    Plot a periodic minimum-energy pathway.

    Parameters
    ----------
    MEP : dict
        {
            "path": [
                {
                    "start_basin": int,
                    "end_basin": int,
                    "transition_state": int,
                    "crossing": -1 | 0 | 1
                },
                ...
            ]
        }

    direction : str
        "a", "b", or "c"

    plotting_data : dict
        Project plotting data.

    basins : dict
        Mapping basin ID -> basin data dictionary.

    TS_data : dict
        Mapping transition-state ID -> TS data dictionary.

    grid_shape : tuple/list
        [Na, Nb, Nc]
    """

    path = MEP.get("path", [])

    if not path:
        return

    dir_idx = {"a": 0, "b": 1, "c": 2}[direction]

    # --------------------------------------------------
    # Ordered basin sequence
    # --------------------------------------------------

    basin_ids = [str(path[0]["start_basin"])]

    for step in path[:-1]:
        basin_ids.append(str(step["end_basin"]))

    ts_ids = [
        str(step["transition_state"])
        for step in path
    ]

    crossings = [
        step["crossing"]
        for step in path
    ]

    # --------------------------------------------------
    # Basin properties
    # --------------------------------------------------

    basin_E = [
        basins[B]["E_min"]
        for B in basin_ids
    ]

    basin_x = [
        basins[B]["center"][dir_idx] / grid_shape[dir_idx]
        for B in basin_ids
    ]

    basin_colors = [
        plotting_data["basins"]
        .get(str(B), {})
        .get("color", "black")
        for B in basin_ids
    ]

    # --------------------------------------------------
    # Transition-state properties
    # --------------------------------------------------

    TS_E = [
        TS_data[TS]["E_min"]
        for TS in ts_ids
    ]

    TS_colors = [
        plotting_data["TS"]
        .get(str(TS), {})
        .get("color", "black")
        for TS in ts_ids
    ]

    # --------------------------------------------------
    # Plot setup
    # --------------------------------------------------

    fig, ax = plt.subplots(figsize=(15, 5))

    basin_width = 0.25 / max(len(basin_ids), 1)

    y_max = max(TS_E) * 1.1
    y_min = min(basin_E) - max(TS_E) * 0.1
    y_range = y_max - y_min

    # --------------------------------------------------
    # Basin plateaus
    # --------------------------------------------------

    for B, E, x, color in zip(
        basin_ids,
        basin_E,
        basin_x,
        basin_colors,
    ):

        ax.hlines(
            E,
            x - basin_width,
            x + basin_width,
            linewidth=3,
            color=color,
        )

        # Periodic copies at the boundaries
        if x - basin_width < 0:
            ax.hlines(
                E,
                x - basin_width + 1,
                x + basin_width + 1,
                linewidth=3,
                color=color,
            )
            label_x = x + 0.01

        elif x + basin_width > 1:
            ax.hlines(
                E,
                x - basin_width - 1,
                x + basin_width - 1,
                linewidth=3,
                color=color,
            )
            label_x = x - 0.01

        else:
            label_x = x

        ha = (
            "left" if label_x < 0.08
            else "right" if label_x > 0.92
            else "center"
        )

        ax.text(
            label_x,
            E + 0.015 * y_range,
            f"B{B}",
            ha=ha,
            fontsize=9,
        )

        ax.text(
            label_x,
            E - 0.05 * y_range,
            f"E={E:.2f}",
            ha=ha,
            fontsize=8,
            color=color,
        )

    # --------------------------------------------------
    # MEP barriers
    # --------------------------------------------------

    n_basins = len(basin_ids)

    for i, (TS, crossing, y_ts, color) in enumerate(
        zip(ts_ids, crossings, TS_E, TS_colors)
    ):

        j = (i + 1) % n_basins

        x0 = basin_x[i] + basin_width
        x2 = basin_x[j] - basin_width

        y0 = basin_E[i]
        y2 = basin_E[j]

        y1 = 2 * y_ts - 0.5 * y0 - 0.5 * y2

        # ------------------------------------------
        # No PBC crossing
        # ------------------------------------------

        if crossing == 0:

            if x0 > x2:
                x0, x2 = x2, x0
                y0, y2 = y2, y0

            xm = 0.5 * (x0 + x2)

            P0 = np.array([x0, y0])
            P1 = np.array([xm, y1])
            P2 = np.array([x2, y2])

            curve = np.array([
                bezier_quad(P0, P1, P2, t)
                for t in np.linspace(0, 1, 100)
            ])

            ax.plot(
                curve[:, 0],
                curve[:, 1],
                color=color,
            )

        # ------------------------------------------
        # +PBC crossing
        # ------------------------------------------

        elif crossing == 1:

            x2_shift = x2 + 1
            xm = 0.5 * (x0 + x2_shift)

            curve = np.array([
                bezier_quad(
                    np.array([x0, y0]),
                    np.array([xm, y1]),
                    np.array([x2_shift, y2]),
                    t,
                )
                for t in np.linspace(0, 1, 100)
            ])

            ax.plot(
                curve[:, 0],
                curve[:, 1],
                color=color,
            )

            ax.plot(
                curve[:, 0] - 1,
                curve[:, 1],
                color=color,
            )

        # ------------------------------------------
        # -PBC crossing
        # ------------------------------------------

        elif crossing == -1:

            x2_shift = x2 - 1
            xm = 0.5 * (x0 + x2_shift)

            curve = np.array([
                bezier_quad(
                    np.array([x2_shift, y2]),
                    np.array([xm, y1]),
                    np.array([x0, y0]),
                    t,
                )
                for t in np.linspace(0, 1, 100)
            ])

            ax.plot(
                curve[:, 0],
                curve[:, 1],
                color=color,
            )

            ax.plot(
                curve[:, 0] + 1,
                curve[:, 1],
                color=color,
            )

        else:
            raise ValueError(
                f"Invalid PBC crossing value: {crossing}"
            )

        ax.text(
            xm % 1,
            y_ts + 0.01 * y_range,
            f"TS{TS}",
            ha="center",
            fontsize=8,
            color=color,
        )

    # --------------------------------------------------
    # Formatting
    # --------------------------------------------------

    ax.set_xlabel(
        f"Fractional coordinate ({direction})"
    )

    ax.set_ylabel(
        "Energy [kJ/mol]"
    )

    ax.set_title(
        f"Minimum Energy Path ({direction}-direction)"
    )

    ax.set_xlim(0, 1)
    ax.set_ylim(y_min, y_max)

    ax.grid(
        axis="y",
        linestyle="--",
        alpha=0.3,
    )

    fig.tight_layout()
    plt.show(block=False)