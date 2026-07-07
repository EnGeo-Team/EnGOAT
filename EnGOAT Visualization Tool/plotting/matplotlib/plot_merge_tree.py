import matplotlib.pyplot as plt

class Node:
    def __init__(
        self,
        name,
        energy,
        left=None,
        right=None,
        merge_color="black"
    ):
        self.name = name
        self.energy = energy
        self.left = left
        self.right = right
        self.merge_color = merge_color
        self.x = None

def filter_redundant_merges(merges):

    parent = {}

    def find(x):

        if parent[x] != x:
            parent[x] = find(parent[x])

        return parent[x]

    def union(x, y):

        parent[find(y)] = find(x)

    nodes = set()

    for b1, b2, _, _ in merges:
        nodes.update([b1, b2])

    for n in nodes:
        parent[n] = n

    filtered = []

    for b1, b2, E, color in merges:

        if b2 == -1:
            continue

        if find(b1) != find(b2):

            filtered.append((b1, b2, E, color))
            union(b1, b2)

    return filtered


def extract_data(project, basin_ids):

    births = {}

    basin_colors = {}

    for basin in project.Basin_list:

        if basin.ID in basin_ids:

            births[basin.ID] = basin.E
            basin_colors[basin.ID] = basin.color

    merges = []

    for ts in project.TS_list:

        if (
            ts.B_start in basin_ids
            and ts.B_end in basin_ids
        ):

            merges.append(
                (
                    ts.B_start,
                    ts.B_end,
                    ts.E,
                    ts.color
                )
            )

    return births, basin_colors, merges

def extract_tree_data(project, basin_ids, tunnel_name=None):

    births = {}
    basin_colors = {}

    # -------------------------
    # Default colors
    # -------------------------
    for basin in project.Basin_list:

        if basin.ID in basin_ids:

            births[basin.ID] = basin.E
            basin_colors[basin.ID] = basin.color

    # -------------------------
    # Override with MEP colors
    # -------------------------
    if tunnel_name is not None:

        plotting = project.tunnel_systems_plotting[tunnel_name]

        if plotting["show_MEP"]:

            direction = plotting["MEP_direction"]

            mep = project.tunnel_systems[tunnel_name]["MEPs"][direction]

            if mep is not None:

                for basin_id, color in zip(
                    mep["basin_ids"],
                    mep["basin_colors"]
                ):

                    if basin_id in basin_colors:
                        basin_colors[basin_id] = color

    merges = []

    for ts in project.TS_list:

        if (
            ts.B_start in basin_ids
            and ts.B_end in basin_ids
            and ts.B_end != -1
        ):

            merges.append(
                (
                    ts.B_start,
                    ts.B_end,
                    ts.E,
                    ts.color
                )
            )

    return births, basin_colors, merges

def build_tree(births, merges):

    nodes = {
        b: Node(b, births[b])
        for b in births
    }

    cluster_id = 0

    parent = None

    for b1, b2, E, color in merges:

        n1 = nodes[b1]
        n2 = nodes[b2]

        parent = Node(
            f"C{cluster_id}",
            E,
            n1,
            n2,
            color
        )

        cluster_id += 1

        nodes[parent.name] = parent

        for k in list(nodes.keys()):

            if nodes[k] in (n1, n2):
                nodes[k] = parent

    return parent

def assign_x(node, x=0):

    if node.left is None:

        node.x = x

        return x + 1

    x = assign_x(node.left, x)
    x = assign_x(node.right, x)

    node.x = (
        node.left.x +
        node.right.x
    ) / 2

    return x

def plot_tree(node, basin_colors):

    if node.left is None:

        color = basin_colors[node.name]

        plt.scatter(
            node.x,
            node.energy,
            s=40,
            color=color,
            zorder=3
        )

        plt.text(
            node.x,
            node.energy - 0.5,
            str(node.name),
            ha="center",
            va="top",
            fontsize=8
        )

        return

    plt.plot(
        [node.left.x, node.left.x],
        [node.left.energy, node.energy],
        color=node.merge_color,
        linewidth=1.5
    )

    plt.plot(
        [node.right.x, node.right.x],
        [node.right.energy, node.energy],
        color=node.merge_color,
        linewidth=1.5
    )

    plt.plot(
        [node.left.x, node.right.x],
        [node.energy, node.energy],
        color=node.merge_color,
        linewidth=1.5
    )

    plot_tree(node.left, basin_colors)
    plot_tree(node.right, basin_colors)


def create_merge_trees(
    project,
    selected_tunnels,
    include_isolated_clusters=False
):

    plt.figure(
        figsize=(8, 6),
        dpi=150
    )

    x = 0

    groups = []

    for tunnel_name in selected_tunnels:

        groups.append(
            (
                tunnel_name,
                project.tunnel_systems[tunnel_name]
            )
        )

    if include_isolated_clusters:

        groups.append(
            (
                "Isolated clusters",
                {
                    "basins": [
                        b
                        for cluster in project.isolated_clusters.values()
                        for b in cluster["basins"]
                    ]
                }
            )
        )

    for label, group in groups:

        basin_ids = set(group["basins"])

        births, basin_colors, merges = extract_tree_data(
            project,
            basin_ids,
            tunnel_name
        )


        if not births:
            continue

        merges.sort(
            key=lambda entry: entry[2]
        )

        merges = filter_redundant_merges(
            merges
        )

        if len(merges) > 0:

            root = build_tree(
                births,
                merges
            )

            x = assign_x(root, x)

            label_y = root.energy + 1.0

            plt.text(
                root.x,
                label_y,
                label,
                ha="center",
                va="bottom",
                fontsize=10,
                fontweight="bold"
            )

            plot_tree(
                root,
                basin_colors
            )

            x += 1

        else:

            center_x = x + (len(births) - 1) / 2
            
            plt.text(
                center_x,
                max(births.values()) + 1.0,
                label,
                ha="center",
                va="bottom",
                fontsize=10,
                fontweight="bold"
            )


            for basin_id, energy in births.items():

                plt.scatter(
                    x,
                    energy,
                    color=basin_colors[basin_id],
                    s=40
                )

                plt.text(
                    x,
                    energy - 0.5,
                    str(basin_id),
                    ha="center",
                    va="top",
                    fontsize=8
                )

                x += 1

    plt.ylabel("Energy (kJ/mol)")
    plt.xticks([])
    plt.xlabel("")

    plt.grid(
        axis="y",
        linestyle="--",
        alpha=0.3
    )

    ax = plt.gca()

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.margins(x=0.05)

    plt.tight_layout()
    plt.show(block=False)