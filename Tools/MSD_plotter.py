#!/usr/bin/env python3
import os
import re
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

def main():

    os.makedirs("MSD_plots", exist_ok=True)

    folder = "MSD_files"

    # Regex pattern for MSD files
    pattern = re.compile(
        r"msd_T(?P<T>[0-9.]+)_tunnel\d+_run(?P<run>\d+)(_direction(?P<dir>\w+))?\.dat"
    )

    # --- Scan files ---
    data = {}
    for filename in os.listdir(folder):
        match = pattern.match(filename)
        if match:
            T = float(match.group("T"))
            direction = match.group("dir") or "none"
            data.setdefault(T, {}).setdefault(direction, []).append(
                os.path.join(folder, filename)
            )

    # --- Plot each temperature and direction ---
    for temperature in sorted(data.keys()):
        for direction in data[temperature]:
            fig, ax = plt.subplots()

            # Inset for log-log plot
            ax_inset = inset_axes(
                ax,
                width="30%",
                height="30%",
                loc="lower right",
                bbox_to_anchor=(-0.05, 0.075, 1, 1),
                bbox_transform=ax.transAxes
            )

            files = sorted(data[temperature][direction])
            colors = plt.cm.tab20(np.linspace(0, 1, len(files)))

            slopes = []

            # Plot each run
            for idx, msd_file in enumerate(files):
                msd_data = np.loadtxt(msd_file)
                time = msd_data[:, 0]
                msd = msd_data[:, 1]

                mask = (time > 0) & (msd > 0)
                color = colors[idx]

                # Linear main plot
                ax.plot(time, msd, color=color, label=f"run {idx+1}")

                # Log-log inset plot
                log_t = np.log(time[mask])
                log_msd = np.log(msd[mask])
                ax_inset.plot(log_t, log_msd, color=color)

                # Fit slope in log-log
                if len(log_t) > 2:
                    slope, _ = np.polyfit(log_t, log_msd, 1)
                    slopes.append(slope)

            # --- Average slope + guide line in inset ---
            if slopes:
                k_avg = np.mean(slopes)

                # Use last dataset for reference
                t_ref = time[mask]
                log_t_ref = np.log(t_ref)
                log_msd_ref = np.log(msd[mask])

                mid = len(log_t_ref) // 2
                t0 = log_t_ref[mid]
                msd0 = log_msd_ref[mid]

                log_t_line = np.linspace(log_t_ref.min(), log_t_ref.max(), 100)
                log_msd_line = msd0 + k_avg * (log_t_line - t0)

            # Main plot labels
            ax.set_xlabel("t [s]")
            ax.set_ylabel(fr"$\langle {direction[0]} \rangle^2$ [$Å^2$]")
            ax.set_title(f"T = {temperature}, direction = {direction[0]}")
            ax.legend()

            # Inset settings
            ax_inset.set_title("log-log plot", fontsize=8)
            ax_inset.tick_params(axis='both', which='major', labelsize=8)

            # Save figure
            figpath = os.path.join(
                "MSD_plots", f"msd_T{temperature}_direction{direction}.png"
            )
            plt.savefig(figpath, dpi=300)
            plt.close()

if __name__ == "__main__":
    main()