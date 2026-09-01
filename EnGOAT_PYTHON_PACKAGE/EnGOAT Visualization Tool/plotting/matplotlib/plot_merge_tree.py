import matplotlib.pyplot as plt
import textwrap


# ============================================================
# Tree node
# ============================================================

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

        # Assigned later
        self.x = None


# ============================================================
# Remove redundant/cyclic merges
# ============================================================

def filter_redundant_merges(merges):

    parent = {}

    def find(x):

        if x not in parent:
            parent[x] = x

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

    # Process merges in energy order.
    # Sorting here makes this function deterministic even if
    # the input list was not sorted.
    merges = sorted(
        merges,
        key=lambda entry: entry[2],
    )

    filtered = []

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

    # --------------------------------------------------------
    # Create all leaf nodes first
    # --------------------------------------------------------

    nodes = {
        basin_id: Node(
            name=basin_id,
            energy=energy,
        )
        for basin_id, energy in births.items()
    }

    # Every basin initially belongs to itself
    cluster_of = {
        basin_id: basin_id
        for basin_id in births
    }

    cluster_counter = 0

    # --------------------------------------------------------
    # Process merges from low to high energy
    # --------------------------------------------------------

    merges = sorted(
        merges,
        key=lambda entry: entry[2],
    )

    for b1, b2, energy, color in merges:

        # Ignore merges referring to missing basins
        if b1 not in cluster_of or b2 not in cluster_of:
            continue

        cluster_1 = cluster_of[b1]
        cluster_2 = cluster_of[b2]

        # Already merged
        if cluster_1 == cluster_2:
            continue

        node_1 = nodes[cluster_1]
        node_2 = nodes[cluster_2]

        parent_name = f"C{cluster_counter}"
        cluster_counter += 1

        parent = Node(
            name=parent_name,
            energy=energy,
            left=node_1,
            right=node_2,
            merge_color=color,
        )

        nodes[parent_name] = parent

        # Update every basin belonging to either cluster
        for basin_id in cluster_of:

            if cluster_of[basin_id] == cluster_1:

                cluster_of[basin_id] = parent_name

            elif cluster_of[basin_id] == cluster_2:

                cluster_of[basin_id] = parent_name

    # --------------------------------------------------------
    # Find remaining disconnected roots
    # --------------------------------------------------------

    roots = sorted(
        set(cluster_of.values()),
        key=str,
    )

    if not roots:
        return None

    # One connected component
    if len(roots) == 1:
        return nodes[roots[0]]

    # --------------------------------------------------------
    # Connect disconnected components with artificial roots
    # --------------------------------------------------------

    current = nodes[roots[0]]

    for root_name in roots[1:]:

        next_node = nodes[root_name]

        # Use an energy slightly above both components.
        #
        # This avoids zero-height artificial branches when
        # both components happen to have exactly the same energy.
        highest_energy = max(
            current.energy,
            next_node.energy,
        )

        energy_offset = max(
            abs(highest_energy) * 1e-6,
            1e-6,
        )

        artificial_energy = (
            highest_energy
            + energy_offset
        )

        current = Node(
            name=f"C{cluster_counter}",
            energy=artificial_energy,
            left=current,
            right=next_node,
            merge_color="black",
        )

        cluster_counter += 1

    return current


# ============================================================
# Count leaves
# ============================================================

def count_leaves(node):

    if node is None:
        return 0

    if node.left is None:
        return 1

    return (
        count_leaves(node.left)
        + count_leaves(node.right)
    )


# ============================================================
# Assign x coordinates
# ============================================================

def assign_x(
    node,
    x=0,
):

    if node is None:
        return x

    # --------------------------------------------------------
    # Leaf
    # --------------------------------------------------------

    if node.left is None:

        node.x = x

        return x + 1

    # --------------------------------------------------------
    # Children
    # --------------------------------------------------------

    x = assign_x(
        node.left,
        x,
    )

    x = assign_x(
        node.right,
        x,
    )

    # Parent sits exactly between children
    node.x = (
        node.left.x
        + node.right.x
    ) / 2.0

    return x


# ============================================================
# Get maximum x coordinate
# ============================================================

def get_max_x(node):

    if node is None:
        return None

    if node.left is None:
        return node.x

    return max(
        get_max_x(node.left),
        get_max_x(node.right),
    )


# ============================================================
# Get minimum x coordinate
# ============================================================

def get_min_x(node):

    if node is None:
        return None

    if node.left is None:
        return node.x

    return min(
        get_min_x(node.left),
        get_min_x(node.right),
    )


# ============================================================
# Get all leaf nodes
# ============================================================

def get_leaf_nodes(node):

    if node is None:
        return []

    if node.left is None:
        return [node]

    return (
        get_leaf_nodes(node.left)
        + get_leaf_nodes(node.right)
    )


# ============================================================
# Plot tree
# ============================================================

def plot_tree(
    node,
    basin_colors,
    label_offset,
):

    if node is None:
        return

    # --------------------------------------------------------
    # Leaf
    # --------------------------------------------------------

    if node.left is None:

        color = basin_colors.get(
            node.name,
            "black",
        )

        # Basin point
        plt.scatter(
            node.x,
            node.energy,
            s=45,
            color=color,
            edgecolor="black",
            linewidth=0.4,
            zorder=5,
        )

        # Basin ID
        plt.text(
            node.x,
            node.energy - label_offset,
            f"B{node.name}",
            ha="center",
            va="top",
            fontsize=8,
            zorder=6,
        )

        return

    # --------------------------------------------------------
    # Vertical branches
    # --------------------------------------------------------

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
        linewidth=1.6,
        solid_capstyle="round",
        zorder=2,
    )

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
        linewidth=1.6,
        solid_capstyle="round",
        zorder=2,
    )

    # --------------------------------------------------------
    # Horizontal merge line
    # --------------------------------------------------------

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
        linewidth=1.6,
        solid_capstyle="round",
        zorder=2,
    )

    # --------------------------------------------------------
    # Children
    # --------------------------------------------------------

    plot_tree(
        node.left,
        basin_colors,
        label_offset,
    )

    plot_tree(
        node.right,
        basin_colors,
        label_offset,
    )


# ============================================================
# Collect energy values from a tree
# ============================================================

def get_tree_energies(node):

    if node is None:
        return []

    values = [node.energy]

    if node.left is not None:
        values.extend(
            get_tree_energies(node.left)
        )

    if node.right is not None:
        values.extend(
            get_tree_energies(node.right)
        )

    return values


# ============================================================
# Main plotting function
# ============================================================

def create_merge_trees(project):

    # ========================================================
    # Collect visible groups
    # ========================================================

    groups = []

    # --------------------------------------------------------
    # Tunnel systems
    # --------------------------------------------------------

    visible_tunnel_systems = project.visibility.get(
        "tunnel_systems",
        {},
    )

    for group_id, visible in sorted(
        visible_tunnel_systems.items(),
        key=lambda item: str(item[0]),
    ):

        if not visible:
            continue

        if group_id not in project.tunnel_systems:
            continue

        group = project.tunnel_systems[group_id]

        groups.append(
            (
                f"Tunnel system {group_id}",
                group.get("basin_list", []),
                group.get("TS_list", []),
            )
        )

    # --------------------------------------------------------
    # Isolated groups
    # --------------------------------------------------------

    visible_isolated_groups = project.visibility.get(
        "isolated_groups",
        {},
    )

    for group_id, visible in sorted(
        visible_isolated_groups.items(),
        key=lambda item: str(item[0]),
    ):

        if not visible:
            continue

        if group_id not in project.isolated_groups:
            continue

        group = project.isolated_groups[group_id]

        groups.append(
            (
                f"Isolated group {group_id}",
                group.get("basin_list", []),
                group.get("TS_list", []),
            )
        )

    # --------------------------------------------------------
    # Nothing to plot
    # --------------------------------------------------------

    if not groups:

        print("No visible tunnel systems or isolated groups.")

        return None

    # ========================================================
    # Prepare data for every group FIRST
    #
    # This is important because it allows us to calculate
    # global energy limits before positioning titles.
    # ========================================================

    prepared_groups = []

    all_energy_values = []

    for label, basin_ids, TS_ids in groups:

        # ----------------------------------------------------
        # Normalize and sort IDs
        # ----------------------------------------------------

        basin_ids = sorted(
            set(basin_ids),
            key=str,
        )

        TS_ids = sorted(
            set(TS_ids),
            key=str,
        )

        # ----------------------------------------------------
        # Basin birth energies and colors
        # ----------------------------------------------------

        births = {}
        basin_colors = {}

        for basin_id in basin_ids:

            basin_key = str(basin_id)

            basin_data = project.basin_data.get(
                basin_key
            )

            if basin_data is None:
                print(
                    f"Warning: basin {basin_id} "
                    f"not found in project.basin_data."
                )
                continue

            if "E_min" not in basin_data:
                print(
                    f"Warning: basin {basin_id} "
                    f"has no E_min."
                )
                continue

            births[basin_id] = basin_data["E_min"]

            basin_plot_data = (
                project.plotting_data
                .get("basins", {})
                .get(basin_key, {})
            )

            basin_colors[basin_id] = (
                basin_plot_data.get(
                    "color",
                    "black",
                )
            )

        # ----------------------------------------------------
        # Nothing valid in this group
        # ----------------------------------------------------

        if not births:
            continue

        # ----------------------------------------------------
        # Extract merges
        # ----------------------------------------------------

        merges = []

        for TS_id in TS_ids:

            TS_key = str(TS_id)

            TS = project.TS_data.get(
                TS_key
            )

            if TS is None:
                print(
                    f"Warning: TS {TS_id} "
                    f"not found in project.TS_data."
                )
                continue

            TS_basins = TS.get(
                "basins",
                [],
            )

            # Only two-basin TSs can be represented
            # by this binary tree.
            if len(TS_basins) != 2:
                continue

            basin_1, basin_2 = TS_basins

            if (
                basin_1 not in births
                or basin_2 not in births
            ):
                continue

            if "E_min" not in TS:
                print(
                    f"Warning: TS {TS_id} "
                    f"has no E_min."
                )
                continue

            TS_plot_data = (
                project.plotting_data
                .get("TS", {})
                .get(TS_key, {})
            )

            TS_color = TS_plot_data.get(
                "color",
                "black",
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
        # Remove cyclic/redundant merges
        # ----------------------------------------------------

        merges = filter_redundant_merges(
            merges
        )

        # ----------------------------------------------------
        # Create tree if possible
        # ----------------------------------------------------

        root = None

        if len(births) > 1:

            root = build_tree(
                births,
                merges,
            )

        # ----------------------------------------------------
        # Collect energies
        # ----------------------------------------------------

        energy_values = list(
            births.values()
        )

        energy_values.extend(
            merge[2]
            for merge in merges
        )

        if root is not None:
            energy_values.extend(
                get_tree_energies(root)
            )

        all_energy_values.extend(
            energy_values
        )

        prepared_groups.append(
            {
                "label": label,
                "births": births,
                "basin_colors": basin_colors,
                "merges": merges,
                "root": root,
            }
        )

    # ========================================================
    # Check whether anything survived
    # ========================================================

    if not prepared_groups:

        print("No valid groups to plot.")

        return None

    # ========================================================
    # Global energy range
    # ========================================================

    global_energy_min = min(
        all_energy_values
    )

    global_energy_max = max(
        all_energy_values
    )

    global_energy_range = (
        global_energy_max
        - global_energy_min
    )

    # Protect against all energies being identical
    if global_energy_range <= 0:

        global_energy_range = max(
            abs(global_energy_max) * 0.05,
            1.0,
        )

    # --------------------------------------------------------
    # Vertical spacing
    # --------------------------------------------------------

    label_offset = max(
        global_energy_range * 0.025,
        0.05,
    )

    title_offset = max(
        global_energy_range * 0.055,
        0.5,
    )

    top_padding = max(
        global_energy_range * 0.14,
        2.0,
    )

    bottom_padding = max(
        global_energy_range * 0.08,
        1.0,
    )

    # ========================================================
    # Figure width
    # ========================================================

    total_leaves = 0

    for group in prepared_groups:

        if group["root"] is not None:

            total_leaves += count_leaves(
                group["root"]
            )

        else:

            total_leaves += len(
                group["births"]
            )

    # Dynamic width.
    #
    # A minimum width keeps small plots readable,
    # while larger datasets get more horizontal space.
    group_width = max(
        2.2,
        total_leaves * 0.65,
    )

    figure_width = max(
        10,
        group_width,
    )

    figure_height = 6.5

    fig, ax = plt.subplots(
        figsize=(
            figure_width,
            figure_height,
        ),
        dpi=150,
    )

    # ========================================================
    # Plot groups
    # ========================================================

    x = 0.0

    group_gap = max(
        2.0,
        total_leaves * 0.03,
    )

    title_y = (
        global_energy_max
        + title_offset
    )

    for group in prepared_groups:

        label = group["label"]
        births = group["births"]
        basin_colors = group["basin_colors"]
        merges = group["merges"]
        root = group["root"]

        # ----------------------------------------------------
        # Merged tree
        # ----------------------------------------------------

        if root is not None:

            # Assign x positions starting at current x
            next_x = assign_x(
                root,
                x,
            )

            group_left_x = get_min_x(
                root
            )

            group_right_x = get_max_x(
                root
            )

            # Plot tree
            plot_tree(
                root,
                basin_colors,
                label_offset,
            )

            # ------------------------------------------------
            # Title
            # ------------------------------------------------

            wrapped_label = "\n".join(
                textwrap.wrap(
                    str(label),
                    width=24,
                    break_long_words=False,
                    break_on_hyphens=False,
                )
            )

            group_center_x = (
                group_left_x
                + group_right_x
            ) / 2.0

            ax.text(
                group_center_x,
                title_y,
                wrapped_label,
                ha="center",
                va="bottom",
                fontsize=10,
                fontweight="bold",
                clip_on=False,
                zorder=10,
            )

            # ------------------------------------------------
            # Move x position for next group
            # ------------------------------------------------

            x = (
                group_right_x
                + group_gap
            )

        # ----------------------------------------------------
        # No merge tree
        #
        # Plot all basins individually.
        # ----------------------------------------------------

        else:

            basin_items = sorted(
                births.items(),
                key=lambda item: str(item[0]),
            )

            if not basin_items:
                continue

            group_left_x = x

            for basin_id, energy in basin_items:

                color = basin_colors.get(
                    basin_id,
                    "black",
                )

                ax.scatter(
                    x,
                    energy,
                    s=45,
                    color=color,
                    edgecolor="black",
                    linewidth=0.4,
                    zorder=5,
                )

                ax.text(
                    x,
                    energy - label_offset,
                    f"B{basin_id}",
                    ha="center",
                    va="top",
                    fontsize=8,
                    zorder=6,
                )

                x += 1.0

            group_right_x = x - 1.0

            group_center_x = (
                group_left_x
                + group_right_x
            ) / 2.0

            wrapped_label = "\n".join(
                textwrap.wrap(
                    str(label),
                    width=24,
                    break_long_words=False,
                    break_on_hyphens=False,
                )
            )

            ax.text(
                group_center_x,
                title_y,
                wrapped_label,
                ha="center",
                va="bottom",
                fontsize=10,
                fontweight="bold",
                clip_on=False,
                zorder=10,
            )

            x += group_gap

    # ========================================================
    # Formatting
    # ========================================================

    ax.set_ylabel(
        "Energy (kJ/mol)"
    )

    ax.set_xlabel("")

    # No x tick labels
    ax.set_xticks([])

    # --------------------------------------------------------
    # Y limits
    #
    # Explicit limits ensure titles and basin labels have
    # room and cannot be accidentally clipped.
    # --------------------------------------------------------

    ax.set_ylim(
        global_energy_min - bottom_padding,
        global_energy_max + top_padding,
    )

    # --------------------------------------------------------
    # Grid
    # --------------------------------------------------------

    ax.grid(
        axis="y",
        linestyle="--",
        alpha=0.3,
        linewidth=0.7,
        zorder=0,
    )

    # --------------------------------------------------------
    # Spines
    # --------------------------------------------------------

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # --------------------------------------------------------
    # Horizontal margins
    # --------------------------------------------------------

    ax.margins(
        x=0.03
    )

    # --------------------------------------------------------
    # Make sure manually positioned titles have room
    # --------------------------------------------------------

    fig.subplots_adjust(
        left=0.08,
        right=0.98,
        bottom=0.10,
        top=0.84,
    )

    # --------------------------------------------------------
    # Display
    # --------------------------------------------------------

    plt.show(
        block=False
    )

    return fig, ax