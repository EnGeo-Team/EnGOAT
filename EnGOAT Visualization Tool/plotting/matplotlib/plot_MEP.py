import numpy as np
import matplotlib.pyplot as plt


def bezier_quad(P0, P1, P2, t):
    return (
        (1 - t) ** 2 * P0
        + 2 * (1 - t) * t * P1
        + t**2 * P2
    )


def plot_MEP_energy_diagram(project, tunnel_name):

    direction = (
        project.tunnel_systems_plotting
        [tunnel_name]
        ["MEP_direction"]
    )

    mep = (
        project.tunnel_systems
        [tunnel_name]
        ["MEPs"]
        [direction]
    )

    if mep is None:
        return

    basin_ids = mep["basin_ids"]
    ts_ids = mep["ts_ids"]

    basin_E = mep["basin_E"]
    TS_E = mep["TS_E"]

    basin_coords = mep["basin_coords"]
    PBC_crossing = mep["PBC_crossing"]

    basin_color_map = {
    basin.ID: basin.color
    for basin in project.Basin_list
    }

    colors = [
    basin_color_map.get(
        basin_id,
        mep["basin_colors"][i]
    )
    for i, basin_id in enumerate(basin_ids)
    ]

    plt.figure(figsize=(15, 5))

    n_basins = len(basin_ids)

    basin_width = 0.25 / max(n_basins, 1)

    y_max = np.max(TS_E) * 1.1
    y_min = np.min(basin_E) - np.max(TS_E) * 0.1
    y_range = y_max - y_min

    # -------------------------
    # Basin plateaus
    # -------------------------

    for bid, E, x, color in zip(
        basin_ids,
        basin_E,
        basin_coords,
        colors
    ):

        plt.hlines(
            E,
            x - basin_width,
            x + basin_width,
            linewidth=3,
            color=color
        )

        if x - basin_width < 0:

            plt.hlines(
                E,
                x - basin_width + 1,
                x + basin_width + 1,
                linewidth=3,
                color=color
            )

            label_x = x + 0.01

        elif x + basin_width > 1:

            plt.hlines(
                E,
                x - basin_width - 1,
                x + basin_width - 1,
                linewidth=3,
                color=color
            )

            label_x = x - 0.01

        else:

            label_x = x

        if label_x < 0.08:
            ha = "left"
        elif label_x > 0.92:
            ha = "right"
        else:
            ha = "center"

        plt.text(
            label_x,
            E - 0.05 * y_range,
            f"E = {E:.2f} kJ/mol",
            ha=ha,
            va="bottom",
            fontsize=9,
            color=color
        )

        plt.text(
            label_x,
            E + 0.015 * y_range,
            f"B{bid}",
            ha=ha,
            va="bottom",
            fontsize=9
        )

    # -------------------------
    # TS color lookup
    # -------------------------

    ts_color_map = {
        ts.ID: ts.color
        for ts in project.TS_list
    }

    # -------------------------
    # Bezier barriers
    # -------------------------

    n = len(basin_ids)

    for i in range(n):

        j = (i + 1) % n

        x0 = basin_coords[i] + basin_width
        x2 = basin_coords[j] - basin_width

        y0 = basin_E[i]
        y2 = basin_E[j]

        y_ts = TS_E[i]

        y1 = (
            2 * y_ts
            - 0.5 * y0
            - 0.5 * y2
        )

        ts_color = ts_color_map.get(
            ts_ids[i],
            "black"
        )

        crossing = PBC_crossing[i]

        # -------------------------
        # No PBC crossing
        # -------------------------

        if crossing == 0:

            xm = 0.5 * (x0 + x2)

            if x0 <= x2:

                P0 = np.array([x0, y0])
                P1 = np.array([xm, y1])
                P2 = np.array([x2, y2])

                t = np.linspace(0, 1, 100)

                curve = np.array([
                    bezier_quad(P0, P1, P2, ti)
                    for ti in t
                ])

                plt.plot(
                    curve[:, 0],
                    curve[:, 1],
                    color=ts_color
                )

            else:

                P0 = np.array([x2, y2])
                P1 = np.array([xm, y1])
                P2 = np.array([x0, y0])

                t = np.linspace(0, 1, 100)

                curve = np.array([
                    bezier_quad(P0, P1, P2, ti)
                    for ti in t
                ])

                plt.plot(
                    curve[:, 0],
                    curve[:, 1],
                    color=ts_color
                )


            label_x = xm % 1

            if label_x < 0.08:
                ha = "left"
            elif label_x > 0.92:
                ha = "right"
            else:
                ha = "center"

            plt.text(
                label_x,
                y_ts + y_range * 0.01,
                f"E = {y_ts:.2f}",
                ha=ha,
                va="bottom",
                fontsize=8,
                color=ts_color
            )

        # -------------------------
        # Positive PBC crossing
        # -------------------------

        elif crossing == 1:

            x2 += 1

            xm = 0.5 * (x0 + x2)

            P0 = np.array([x0, y0])
            P1 = np.array([xm, y1])
            P2 = np.array([x2, y2])

            t = np.linspace(0, 1, 100)

            curve = np.array([
                bezier_quad(P0, P1, P2, ti)
                for ti in t
            ])

            plt.plot(
                curve[:, 0],
                curve[:, 1],
                color=ts_color
            )

            curve[:, 0] -= 1

            plt.plot(
                curve[:, 0],
                curve[:, 1],
                color=ts_color
            )


            label_x = xm % 1

            if label_x < 0.08:
                ha = "left"
            elif label_x > 0.92:
                ha = "right"
            else:
                ha = "center"


            plt.text(
                label_x,
                y_ts + y_range * 0.01,
                f"E = {y_ts:.2f}",
                ha=ha,
                va="bottom",
                fontsize=8,
                color=ts_color
            )

        # -------------------------
        # Negative PBC crossing
        # -------------------------

        elif crossing == -1:

            x2 -= 1

            xm = 0.5 * (x0 + x2)

            P0 = np.array([x0, y0])
            P1 = np.array([xm, y1])
            P2 = np.array([x2, y2])

            t = np.linspace(0, 1, 100)

            curve = np.array([
                bezier_quad(P2, P1, P0, ti)
                for ti in t
            ])

            plt.plot(
                curve[:, 0],
                curve[:, 1],
                color=ts_color
            )

            curve[:, 0] += 1

            plt.plot(
                curve[:, 0],
                curve[:, 1],
                color=ts_color
            )


            label_x = xm % 1

            if label_x < 0.08:
                ha = "left"
            elif label_x > 0.92:
                ha = "right"
            else:
                ha = "center"


            plt.text(
                label_x,
                y_ts + y_range * 0.01,
                f"E = {y_ts:.2f}",
                ha=ha,
                va="bottom",
                fontsize=8,
                color=ts_color
            )

    plt.xlabel(
        f"Fractional coordinate {direction}"
    )

    plt.ylabel(
        "Energy [kJ/mol]"
    )

    plt.title(
        f"{tunnel_name} — MEP in {direction} Direction"
    )

    plt.xlim(0, 1)
    plt.ylim(y_min, y_max)

    plt.grid(
        axis="y",
        linestyle="--",
        alpha=0.3
    )

    for spine in plt.gca().spines.values():
        spine.set_visible(True)

    plt.tight_layout()
    plt.show(block=False)