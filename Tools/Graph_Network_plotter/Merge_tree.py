import matplotlib.pyplot as plt
import read_basin_data
import os
import graph_network


# ---------------------------
# NODE CLASS
# ---------------------------
class Node:
    def __init__(self, name, energy, left=None, right=None):
        self.name = name
        self.energy = energy
        self.left = left
        self.right = right
        self.x = None


# ---------------------------
# UNION-FIND FILTER
# ---------------------------
def filter_redundant_merges(merges):
    parent = {}

    def find(x):
        if parent[x] != x:
            parent[x] = find(parent[x])
        return parent[x]

    def union(x, y):
        parent[find(y)] = find(x)

    nodes = set()
    for b1, b2, _ in merges:
        nodes.update([b1, b2])

    for n in nodes:
        parent[n] = n

    filtered = []
    for b1, b2, E in merges:
        if find(b1) != find(b2) and b2 != 0:
            filtered.append((b1, b2, E))
            union(b1, b2)

    return filtered


# ---------------------------
# TREE BUILDING
# ---------------------------
def build_tree(births, merges):
    nodes = {b: Node(b, births[b]) for b in births}
    cluster_id = 0

    for b1, b2, E in merges:
        n1, n2 = nodes[b1], nodes[b2]

        parent = Node(f"C{cluster_id}", E, n1, n2)
        cluster_id += 1

        nodes[parent.name] = parent

        # compress references
        for k in list(nodes.keys()):
            if nodes[k] in (n1, n2):
                nodes[k] = parent

    return parent


# ---------------------------
# X ASSIGNMENT
# ---------------------------
def assign_x(node, x=0):
    if node.left is None:
        node.x = x
        return x + 1

    x = assign_x(node.left, x)
    x = assign_x(node.right, x)

    node.x = (node.left.x + node.right.x) / 2
    return x


# ---------------------------
# PLOTTING
# ---------------------------
def plot_tree(node, color):
    if node.left is None:
        plt.scatter(node.x, node.energy, s=30, color=color, zorder=3)
        plt.text(node.x, node.energy - 0.5, node.name,
                 ha='center', va='top', fontsize=8)
        return

    for c in [node.left, node.right]:
        plt.plot([c.x, c.x], [c.energy, node.energy],
                 color=color, linewidth=1.5)

    plt.plot([node.left.x, node.right.x],
             [node.energy, node.energy],
             color=color, linewidth=1.5)

    plot_tree(node.left, color)
    plot_tree(node.right, color)


# ---------------------------
# LOAD DATA
# ---------------------------
path = os.path.join("TuTraSt_data", "basin_data.dat")
all_basins = read_basin_data.read_basins(path)

path = os.path.join("TuTraSt_data", "TS_data.dat")
all_merges = read_basin_data.read_merges(path)

# normalize types once
all_merges = [(int(a), int(b), e) for a, b, e in all_merges]


# ---------------------------
# GROUP ISOLATED BASINS
# ---------------------------
iso_basins = {str(b.ID) for b in all_basins if b.ts is None}
groups = {b: None for b in iso_basins}

group_id = 0
group_IDs = []

for b in iso_basins:
    if groups[b] is not None:
        continue

    stack = {b}
    changed = True

    while changed:
        changed = False
        for a, c, _ in all_merges:
            if str(a) in stack and str(c) not in stack:
                stack.add(str(c))
                changed = True
            elif str(c) in stack and str(a) not in stack:
                stack.add(str(a))
                changed = True

    for s in stack:
        groups[s] = group_id

    group_IDs.append(group_id)
    group_id += 1


# ---------------------------
# PLOT
# ---------------------------
plt.figure(figsize=(8, 6), dpi=150)
x = 0

#tunnel_systems = set(b.ts for b in all_basins)
tunnel_systems = sorted(set(b.ts for b in all_basins) - {None})
ts_cmap = plt.get_cmap("tab20")

tunnel_system_map = {}
ts_to_color = {
    ts: ts_cmap(i / max(len(tunnel_systems) - 1, 1))
    for i, ts in enumerate(tunnel_systems)
}

# --- tunnel systems ---
for ts in tunnel_systems:
    color = ts_to_color[ts]
    basins_ts = [b for b in all_basins if b.ts == ts]
    tunnel_system_map[ts] = basins_ts

    births = {
        int(b.ID): b.E
        for b in all_basins
        if b.ts == ts
    }

    merges = [
        m for m in all_merges
        if m[0] in births and m[1] in births
    ]

    if len(merges) == 0:
        continue

    merges.sort(key=lambda x: x[2])
    merges = filter_redundant_merges(merges)

    root = build_tree(births, merges)
    x = assign_x(root, x) + 1
    plot_tree(root, color)


# --- isolated groups ---
iso_cmap = plt.get_cmap("Dark2")

iso_to_color = {
    gid: iso_cmap(i / max(len(group_IDs) - 1, 1))
    for i, gid in enumerate(group_IDs)
}
isolated_groups_map = {}

for gid in group_IDs:
    color = iso_to_color[gid]

    iso_group = [b for b in all_basins if str(b.ID) in iso_basins and groups[str(b.ID)] == gid]
    isolated_groups_map[gid] = iso_group

    births = {
        int(b.ID): b.E
        for b in all_basins
        if str(b.ID) in iso_basins and groups[str(b.ID)] == gid
    }

    merges = [
        m for m in all_merges
        if m[0] in births and m[1] in births
    ]

    if len(births) == 0:
        continue

    if len(merges) > 0:
        merges.sort(key=lambda x: x[2])
        merges = filter_redundant_merges(merges)

        root = build_tree(births, merges)
        x = assign_x(root, x) + 1
        plot_tree(root, color)

    else:
        for b, E in births.items():
            plt.scatter(x, E, color=color, s=30)
            plt.text(x, E - 0.5, b,
                 ha='center', va='top', fontsize=8)
        x += 1

# ---------------------------
# STYLE
# ---------------------------
plt.ylabel("Energy (kJ/mol)", fontsize=12)
plt.xticks([])
plt.xlabel("")
plt.grid(axis="y", linestyle="--", alpha=0.3)

ax = plt.gca()
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.margins(x=0.05)
plt.tight_layout()

from matplotlib.patches import Patch

legend_elements = [
    Patch(facecolor='blue', edgecolor='blue', label='Tunnel systems'),
    Patch(facecolor='red', edgecolor='red', label='Inaccessible clusters')
]

plt.legend(handles=legend_elements, loc='best')

os.makedirs("Graph_network_plots", exist_ok=True)

plt.savefig(os.path.join("Graph_network_plots", "Merge_tree.png"), dpi = 300)

graph_network.plot_graph_network(tunnel_system_map, ts_to_color, isolated_groups_map, iso_to_color)