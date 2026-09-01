import numpy as np
import matplotlib.pyplot as plt


def plot_arrhenius(kMC_data, selected_directions):
    R = 8.314462618
    directions = ["a", "b", "c", "x", "y", "z", "3D"]

    temperatures = []
    for temperature in kMC_data:
        try:
            temperatures.append(float(temperature))
        except (TypeError, ValueError):
            continue

    temperatures = np.array(sorted(temperatures))

    if len(temperatures) < 2:
        raise ValueError("At least two temperatures are required for an Arrhenius fit.")

    for direction in selected_directions:
        if direction not in directions:
            raise ValueError(f"Unknown direction '{direction}'. Valid directions are {directions}.")

    fig, ax = plt.subplots(figsize=(8.5, 6))
    tab10 = plt.get_cmap("tab10")
    direction_colors = {direction: tab10(i) for i, direction in enumerate(directions)}
    all_used_temperatures = []

    for direction in selected_directions:
        D_values = []
        sd_values = []
        valid_temperatures = []

        for temperature in temperatures:
            temperature_key = str(int(temperature))

            if temperature_key in kMC_data:
                data = kMC_data[temperature_key]
            elif temperature in kMC_data:
                data = kMC_data[temperature]
            else:
                continue

            try:
                direction_data = data["D_tot"][direction]
                D = float(direction_data["D"])
                sd = float(direction_data["sd"])
            except (KeyError, TypeError, ValueError):
                continue

            if D <= 0:
                continue

            valid_temperatures.append(temperature)
            D_values.append(D)
            sd_values.append(sd)

        if len(D_values) < 2:
            print(f"Skipping {direction}: not enough valid data points for a fit.")
            continue

        valid_temperatures = np.asarray(valid_temperatures)
        D_values = np.asarray(D_values)
        sd_values = np.asarray(sd_values)
        all_used_temperatures.extend(valid_temperatures)

        x = 1000.0 / valid_temperatures
        y = np.log10(D_values)
        yerr = sd_values / (D_values * np.log(10))

        slope, intercept = np.polyfit(x, y, 1)
        Ea = -slope * 2.303 * R

        x_fit = np.linspace(x.min(), x.max(), 200)
        y_fit = slope * x_fit + intercept
        color = direction_colors[direction]

        ax.errorbar(
            x,
            y,
            yerr=yerr,
            fmt="o",
            markersize=7,
            markeredgewidth=1.2,
            markeredgecolor="white",
            color=color,
            ecolor=color,
            elinewidth=1.4,
            capsize=4,
            capthick=1.4,
            label=f"{direction}  ($E_a$ = {Ea:.2f} kJ/mol)",
            zorder=3
        )

        ax.plot(
            x_fit,
            y_fit,
            color=color,
            linewidth=2.0,
            linestyle="--",
            zorder=2
        )

    ax.set_xlabel(r"$1000/T$ [K$^{-1}$]", fontsize=12)
    ax.set_ylabel(r"$\log_{10}(D)$", fontsize=12)

    all_used_temperatures = np.unique(np.asarray(all_used_temperatures))

    if len(all_used_temperatures) > 0:
        ax_top = ax.twiny()
        ax_top.set_xlim(ax.get_xlim())

        temperature_positions = 1000.0 / all_used_temperatures
        ax_top.set_xticks(temperature_positions)
        ax_top.set_xticklabels([f"{temperature:g}" for temperature in all_used_temperatures])
        ax_top.set_xlabel("Temperature [K]", fontsize=12, labelpad=8)
        ax_top.tick_params(axis="x", which="major", labelsize=10)
        ax_top.spines["bottom"].set_visible(False)

    ax.grid(True, which="major", linestyle="--", linewidth=0.7, alpha=0.3)
    ax.set_axisbelow(True)

    ax.set_title("Arrhenius Plot", fontsize=15, fontweight="bold", pad=15)

    ax.legend(
        title="Direction",
        frameon=True,
        fancybox=True,
        framealpha=0.95,
        loc="best"
    )

    ax.tick_params(axis="both", which="major", labelsize=10)
    ax.spines["top"].set_visible(False)

    fig.tight_layout()
    plt.show(block = False)

    return fig, ax