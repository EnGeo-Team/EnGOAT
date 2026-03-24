#!/usr/bin/env python3
import os
import re
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def parse_files(folder):
    """Scan MSD files and group by temperature and direction."""
    pattern = re.compile(
        r"msd_T(?P<T>[0-9.]+)_tunnel\d+_run\d+(?:_direction(?P<dir>\w+)|_(?P<iso>isotropic))\.dat"
    )

    data = {}

    for filename in os.listdir(folder):
        match = pattern.match(filename)
        if not match:
            continue

        T = float(match.group("T"))
        direction = "isotropic" if match.group("iso") else match.group("dir")

        data.setdefault(T, {}).setdefault(direction, []).append(
            os.path.join(folder, filename)
        )

    return data


def get_labels(direction):
    """Return axis label and title string."""
    if direction == "isotropic":
        return r"$\langle r \rangle^2$", "isotropic"
    else:
        axis = direction[0]
        return fr"$\langle {axis} \rangle^2$", axis


def plot_group(temperature, direction, files, output_folder):
    """Plot MSD curves for a given temperature and direction."""
    fig, ax = plt.subplots()

    # Inset (log-log)
    ax_inset = inset_axes(
        ax,
        width="30%",
        height="30%",
        loc="lower right",
        bbox_to_anchor=(-0.05, 0.075, 1, 1),
        bbox_transform=ax.transAxes,
    )

    colors = plt.cm.tab20(np.linspace(0, 1, len(files)))
    slopes = []

    for idx, filepath in enumerate(sorted(files)):
        data = np.loadtxt(filepath)
        time, msd = data[:, 0], data[:, 1]

        mask = (time > 0) & (msd > 0)
        if not np.any(mask):
            continue

        color = colors[idx]

        # Main plot
        ax.plot(time, msd, color=color, label=f"run {idx+1}")

        # Log-log inset
        log_t = np.log(time[mask])
        log_msd = np.log(msd[mask])
        ax_inset.plot(log_t, log_msd, color=color)

        # Slope
        if len(log_t) > 2:
            slope, _ = np.polyfit(log_t, log_msd, 1)
            slopes.append(slope)

    # Labels
    ylabel, title_dir = get_labels(direction)
    ax.set_xlabel("t [s]")
    ax.set_ylabel(f"{ylabel} [$Å^2$]")
    ax.set_title(f"T = {temperature}, direction = {title_dir}")
    ax.legend()

    ax_inset.set_title("log-log", fontsize=8)
    ax_inset.tick_params(labelsize=8)

    # Output filename
    if direction == "isotropic":
        filename = f"msd_T{temperature}_isotropic.png"
    else:
        filename = f"msd_T{temperature}_direction{direction}.png"

    plt.savefig(os.path.join(output_folder, filename), dpi=300)
    plt.close()


def main():
    input_folder = "MSD_files"
    output_folder = "MSD_plots"

    os.makedirs(output_folder, exist_ok=True)

    data = parse_files(input_folder)

    for T in sorted(data):
        for direction, files in data[T].items():
            plot_group(T, direction, files, output_folder)


if __name__ == "__main__":
    main()