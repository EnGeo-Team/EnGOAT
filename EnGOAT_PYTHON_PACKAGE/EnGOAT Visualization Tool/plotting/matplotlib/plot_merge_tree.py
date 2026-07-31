import matplotlib.pyplot as plt


class Node:

    def __init__(
        self,
        name,
        energy,
        left=None,
        right=None,
        merge_color="black",
    ):

        self.name = name
        self.energy = energy
        self.left = left
        self.right = right
        self.merge_color = merge_color
        self.x = None


# ============================================================
# Remove redundant/cyclic merges
# ============================================================

def filter_redundant_merges(merges):

    parent = {}

    def find(x):

        if parent[x] != x:
            parent[x] = find(parent[x])

        return parent[x]

    def union(x, y):

        root_x = find(x)
        root_y = find(y)

        if root_x == root_y:
            return False

        parent[root_y] = root_x

        return True

    # Collect all basin IDs
    nodes = set()

    for b1, b2, _, _ in merges:

        nodes.add(b1)
        nodes.add(b2)

    # Initialize union-find
    for node in nodes:

        parent[node] = node

    filtered = []

    # Process lowest-energy merges first
    for b1, b2, energy, color in merges:

        if union(b1, b2):

            filtered.append(
                (
                    b1,
                    b2,
                    energy,
                    color,
                )
            )

    return filtered


# ============================================================
# Build merge tree
# ============================================================

def build_tree(births, merges):

    # Initial leaf nodes
    nodes = {
        basin_id: Node(
            basin_id,
            energy,
        )
        for basin_id, energy in births.items()
    }

    # Each basin initially belongs to its own cluster
    cluster_of = {
        basin_id: basin_id
        for basin_id in births
    }

    cluster_counter = 0

    for b1, b2, energy, color in merges:

        # Find the current cluster containing each basin
        cluster_1 = cluster_of[b1]
        cluster_2 = cluster_of[b2]

        # If both basins already belong to the same cluster,
        # this merge is redundant
        if cluster_1 == cluster_2:

            continue

        node_1 = nodes[cluster_1]
        node_2 = nodes[cluster_2]

        parent_name = f"C{cluster_counter}"

        parent = Node(
            name=parent_name,
            energy=energy,
            left=node_1,
            right=node_2,
            merge_color=color,
        )

        cluster_counter += 1

        # Store the new cluster
        nodes[parent_name] = parent

        # Every basin in either old cluster now belongs
        # to the new cluster
        for basin_id in cluster_of:

            if cluster_of[basin_id] in (
                cluster_1,
                cluster_2,
            ):

                cluster_of[basin_id] = parent_name

    # Find all remaining root clusters
    roots = set(
        cluster_of.values()
    )

    # If there is only one root, return it directly
    if len(roots) == 1:

        return nodes[next(iter(roots))]

    # If multiple disconnected components remain,
    # connect them with an artificial root
    roots = list(roots)

    current = nodes[roots[0]]

    for root_name in roots[1:]:

        next_node = nodes[root_name]

        current = Node(
            name=f"C{cluster_counter}",
            energy=max(
                current.energy,
                next_node.energy,
            ),
            left=current,
            right=next_node,
            merge_color="black",
        )

        cluster_counter += 1

    return current


# ============================================================
# Assign x coordinates
# ============================================================

def assign_x(node, x=0):

    # Leaf node
    if node.left is None:

        node.x = x

        return x + 1

    x = assign_x(
        node.left,
        x,
    )

    x = assign_x(
        node.right,
        x,
    )

    node.x = (
        node.left.x
        + node.right.x
    ) / 2

    return x


# ============================================================
# Plot tree
# ============================================================

def plot_tree(
    node,
    basin_colors,
):

    # Leaf
    if node.left is None:

        color = basin_colors.get(
            node.name,
            "black",
        )

        plt.scatter(
            node.x,
            node.energy,
            s=40,
            color=color,
            zorder=3,
        )

        plt.text(
            node.x,
            node.energy - 0.5,
            str(node.name),
            ha="center",
            va="top",
            fontsize=8,
        )

        return

    # Left vertical branch
    plt.plot(
        [
            node.left.x,
            node.left.x,
        ],
        [
            node.left.energy,
            node.energy,
        ],
        color=node.merge_color,
        linewidth=1.5,
    )

    # Right vertical branch
    plt.plot(
        [
            node.right.x,
            node.right.x,
        ],
        [
            node.right.energy,
            node.energy,
        ],
        color=node.merge_color,
        linewidth=1.5,
    )

    # Horizontal merge line
    plt.plot(
        [
            node.left.x,
            node.right.x,
        ],
        [
            node.energy,
            node.energy,
        ],
        color=node.merge_color,
        linewidth=1.5,
    )

    # Recursively plot children
    plot_tree(
        node.left,
        basin_colors,
    )

    plot_tree(
        node.right,
        basin_colors,
    )

import matplotlib.pyplot as plt
import textwrap

#MAIN
#MAIN
#MAIN
#MAIN
#MAIN

def create_merge_trees(project):

    plt.figure(
        figsize=(10, 6),
        dpi=150,
    )

    # --------------------------------------------------------
    # Layout parameters
    # --------------------------------------------------------

    x = 0

    group_gap = 5.0

    title_pad_fraction = 0.08

    max_title_width = 22

    groups = []

    # --------------------------------------------------------
    # Visible tunnel systems
    # --------------------------------------------------------

    for group_id, visible in project.visibility[
        "tunnel_systems"
    ].items():

        if not visible:

            continue

        group = project.tunnel_systems[
            group_id
        ]

        groups.append(
            (
                f"Tunnel system {group_id}",
                group["basin_list"],
                group["TS_list"],
            )
        )

    # --------------------------------------------------------
    # Visible isolated groups
    # --------------------------------------------------------

    for group_id, visible in project.visibility[
        "isolated_groups"
    ].items():

        if not visible:

            continue

        group = project.isolated_groups[
            group_id
        ]

        groups.append(
            (
                f"Isolated group {group_id}",
                group["basin_list"],
                group["TS_list"],
            )
        )

    # --------------------------------------------------------
    # Plot each visible group
    # --------------------------------------------------------

    for label, basin_ids, TS_ids in groups:

        basin_ids = set(
            basin_ids
        )

        TS_ids = set(
            TS_ids
        )

        # ----------------------------------------------------
        # Basin birth energies and colors
        # ----------------------------------------------------

        births = {}

        basin_colors = {}

        for basin_id in basin_ids:

            basin_data = project.basin_data[
                str(basin_id)
            ]

            births[basin_id] = basin_data[
                "E_min"
            ]

            basin_colors[basin_id] = (
                project.plotting_data[
                    "basins"
                ]
                .get(
                    str(basin_id),
                    {},
                )
                .get(
                    "color",
                    "black",
                )
            )

        if not births:

            continue

        # ----------------------------------------------------
        # Extract merges
        # ----------------------------------------------------

        merges = []

        for TS_id in TS_ids:

            TS = project.TS_data[
                str(TS_id)
            ]

            TS_basins = TS.get(
                "basins",
                [],
            )

            if len(TS_basins) != 2:

                continue

            basin_1, basin_2 = TS_basins

            if (
                basin_1 not in basin_ids
                or basin_2 not in basin_ids
            ):

                continue

            TS_color = (
                project.plotting_data[
                    "TS"
                ]
                .get(
                    str(TS_id),
                    {},
                )
                .get(
                    "color",
                    "black",
                )
            )

            merges.append(
                (
                    basin_1,
                    basin_2,
                    TS["E_min"],
                    TS_color,
                )
            )

        # ----------------------------------------------------
        # Sort and filter merges
        # ----------------------------------------------------

        merges.sort(
            key=lambda entry: entry[2]
        )

        merges = filter_redundant_merges(
            merges
        )

        # ----------------------------------------------------
        # Plot a merged tree
        # ----------------------------------------------------

        if merges:

            root = build_tree(
                births,
                merges,
            )

            # Remember where this group starts
            group_start_x = x

            # Assign x coordinates
            next_x = assign_x(
                root,
                x,
            )

            # Determine actual group boundaries
            group_left_x = x

            group_right_x = get_max_x(
                root
            )

            # Plot the tree
            plot_tree(
                root,
                basin_colors,
            )

            # ------------------------------------------------
            # Title placement
            # ------------------------------------------------

            energy_values = list(
                births.values()
            )

            energy_values.extend(
                merge[2]
                for merge in merges
            )

            energy_min = min(
                energy_values
            )

            energy_max = max(
                energy_values
            )

            energy_range = (
                energy_max
                - energy_min
            )

            title_pad = max(
                1.0,
                energy_range
                * title_pad_fraction,
            )

            # Wrap long labels
            wrapped_label = "\n".join(
                textwrap.wrap(
                    label,
                    width=max_title_width,
                )
            )

            group_center_x = (
                group_left_x
                + group_right_x
            ) / 2

            plt.text(
                group_center_x,
                energy_max + title_pad,
                wrapped_label,
                ha="center",
                va="bottom",
                fontsize=10,
                fontweight="bold",
            )

            # Leave enough space before next group
            x = (
                group_right_x
                + group_gap
            )

        # ----------------------------------------------------
        # Plot isolated/unmerged basins
        # ----------------------------------------------------

        else:

            basin_items = list(
                births.items()
            )

            group_left_x = x

            group_right_x = (
                x
                + len(basin_items)
                - 1
            )

            group_center_x = (
                group_left_x
                + group_right_x
            ) / 2

            energy_max = max(
                births.values()
            )

            energy_min = min(
                births.values()
            )

            energy_range = (
                energy_max
                - energy_min
            )

            title_pad = max(
                1.0,
                energy_range
                * title_pad_fraction,
            )

            wrapped_label = "\n".join(
                textwrap.wrap(
                    label,
                    width=max_title_width,
                )
            )

            plt.text(
                group_center_x,
                energy_max + title_pad,
                wrapped_label,
                ha="center",
                va="bottom",
                fontsize=10,
                fontweight="bold",
            )

            for basin_id, energy in basin_items:

                plt.scatter(
                    x,
                    energy,
                    color=basin_colors[
                        basin_id
                    ],
                    s=40,
                )

                plt.text(
                    x,
                    energy - 0.5,
                    str(basin_id),
                    ha="center",
                    va="top",
                    fontsize=8,
                )

                x += 1

            # Space before next group
            x += group_gap

    # --------------------------------------------------------
    # Formatting
    # --------------------------------------------------------

    plt.ylabel(
        "Energy (kJ/mol)"
    )

    plt.xticks([])

    plt.xlabel("")

    plt.grid(
        axis="y",
        linestyle="--",
        alpha=0.3,
    )

    ax = plt.gca()

    ax.spines[
        "top"
    ].set_visible(
        False
    )

    ax.spines[
        "right"
    ].set_visible(
        False
    )

    plt.margins(
        x=0.05
    )

    plt.tight_layout()

    plt.show(
        block=False
    )


def get_max_x(node):

    if node.left is None:

        return node.x

    return max(
        get_max_x(
            node.left
        ),
        get_max_x(
            node.right
        ),
    )