import numpy as np
from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QGroupBox, QHBoxLayout,
    QPushButton, QCheckBox, QLabel, QLineEdit, QWidget, QScrollArea, QColorDialog
)
from plotting.matplotlib.plot_histogram import plot_energy_histogram

class BasinSelectionDialog(QDialog):

    def __init__(self, project):
        super().__init__()

        self.project = project

        self.setWindowTitle("Plot individual basins")
        self.resize(600, 700)

        self.basin_widgets = {}
        self.tunnel_widgets = {}
        self.group_basin_widgets = {}

        main_layout = QVBoxLayout(self)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)

        content = QWidget()
        self.layout = QVBoxLayout(content)

        scroll.setWidget(content)
        main_layout.addWidget(scroll)

        # ✅ Apply button
        self.apply_btn = QPushButton("Apply changes")
        self.apply_btn.setFixedHeight(35)
        main_layout.addWidget(self.apply_btn)

        self.apply_btn.clicked.connect(self.apply_changes)

        self.build_ui()


    # ======================
    # ✅ UI BUILDING
    # ======================
    def build_ui(self):

        tunnel_groups = {}
        isolated_basins = []

        for b in self.project.Basin_list:
            if b.tunnel:
                key = f"Tunnel system {b.tunnel}"
                tunnel_groups.setdefault(key, []).append(b)
            else:
                isolated_basins.append(b)

        def extract_number(name):
            return int(name.split()[-1])

        sorted_keys = sorted(tunnel_groups.keys(), key=extract_number)

        for name in sorted_keys:
            self._add_group(name, tunnel_groups[name])

        if isolated_basins:
            self._add_group("Isolated clusters", isolated_basins)

        self.layout.addStretch()


    def _add_group(self, name, basins):

        group = QGroupBox(name)
        group_layout = QVBoxLayout()
        self.group_basin_widgets[name] = []

        # ======================
        # ✅ HEADER
        # ======================
        header_row = QHBoxLayout()

        # get plotting info
        if name == "Isolated clusters":
            info = self.project.isolated_clusters_plotting
        else:
            info = self.project.tunnel_systems_plotting.get(name)

        visible = info["basins"].get("visible", False)
        color = info["basins"].get("color", (0.5, 0.5, 1.0))
        opacity_val = info["basins"].get("opacity", 0.5)

        # toggle
        toggle = QCheckBox(name)
        toggle.setChecked(visible)
        toggle.setStyleSheet("font-weight: bold;")

        # color
        color_btn = QPushButton()
        color_btn.setFixedSize(25, 15)
        self._set_button_color(color_btn, color)


        # opacity
        opacity = QLineEdit(f"{opacity_val:.2f}")
        opacity.setFixedWidth(50)
        hist_btn = QPushButton("Histogram")
        hist_btn.setFixedWidth(90)

        hist_btn.clicked.connect(
            lambda _, n=name: self._histogram_group(n)
        )

        # store header widgets
        self.tunnel_widgets[name] = {
            "toggle": toggle,
            "color": color_btn,
            "opacity": opacity,
            "histogram": hist_btn,
            "object": info
        }

        # ✅ Connect group actions
        toggle.toggled.connect(lambda val, n=name: self._apply_group_toggle(n, val))
        color_btn.clicked.connect(lambda _, n=name, b=color_btn: self._apply_group_color(n, b))
        opacity.textChanged.connect(lambda t, n=name: self._apply_group_opacity(n, t))

        # layout
        header_row.addWidget(toggle)
        header_row.addStretch()
        color_label = QLabel("Color:")
        color_label.setStyleSheet("font-weight: bold;")
        header_row.addWidget(color_label)
        header_row.addWidget(color_btn)
        header_row.addSpacing(10)
        opacity_label = QLabel("Opacity:")
        opacity_label.setStyleSheet("font-weight: bold;")
        header_row.addWidget(opacity_label)
        header_row.addWidget(opacity)
        header_row.addWidget(hist_btn)

        group_layout.addLayout(header_row)

        # ======================
        # ✅ BASINS
        # ======================
        for b in basins:

            row = QHBoxLayout()

            cb = QCheckBox(f"B{b.ID}")
            cb.setChecked(b.visible)

            color_btn = QPushButton()
            color_btn.setFixedSize(25, 15)
            self._set_button_color(color_btn, b.color)
            color_btn.clicked.connect(
            lambda _, btn=color_btn, obj=b: self._change_single_color(btn, obj)
        )

            opacity = QLineEdit(f"{b.opacity:.2f}")
            opacity.setFixedWidth(50)


            hist_btn = QPushButton("Histogram")
            hist_btn.setFixedWidth(90)

            hist_btn.clicked.connect(
                lambda _, basin_id=b.ID: self._histogram_basin(basin_id)
            )


            row.addWidget(cb)
            row.addStretch()
            row.addWidget(QLabel("Color:"))
            row.addWidget(color_btn)
            row.addSpacing(10)
            row.addWidget(QLabel("Opacity:"))
            row.addWidget(opacity)
            row.addWidget(hist_btn)

            group_layout.addLayout(row)

            obj = {
                "checkbox": cb,
                "color": color_btn,
                "opacity": opacity,
                "histogram": hist_btn,
                "object": b
            }

            self.basin_widgets[b.ID] = obj
            self.group_basin_widgets[name].append(obj)

        group.setLayout(group_layout)
        self.layout.addWidget(group)


    # ======================
    # ✅ GROUP ACTIONS
    # ======================
    def _apply_group_toggle(self, name, state):

        for item in self.group_basin_widgets.get(name, []):
            item["checkbox"].setChecked(state)


    def _apply_group_color(self, name, button):

        color = QColorDialog.getColor()
        if not color.isValid():
            return

        rgb = (color.red(), color.green(), color.blue())

        self._set_button_color(button, [c/255 for c in rgb])

        for item in self.group_basin_widgets.get(name, []):
            self._set_button_color(item["color"], [c/255 for c in rgb])


    def _apply_group_opacity(self, name, text):

        for item in self.group_basin_widgets.get(name, []):
            item["opacity"].setText(text)


    # ======================
    # ✅ HELPERS
    # ======================
    def _set_button_color(self, button, color):

        r, g, b = color
        button.setStyleSheet(
            f"border:1px solid #888; border-radius:3px; background-color: rgb({int(r*255)}, {int(g*255)}, {int(b*255)});"
        )


    def _get_button_color(self, button):

        import re
        style = button.styleSheet()
        m = re.search(r"rgb\((\d+),\s*(\d+),\s*(\d+)\)", style)

        if m:
            return (
                int(m.group(1))/255,
                int(m.group(2))/255,
                int(m.group(3))/255
            )
        return (0.5, 0.5, 1.0)
    

    def _change_single_color(self, button, basin_obj):

        color = QColorDialog.getColor()
        if not color.isValid():
            return

        # ✅ update button visually
        button.setStyleSheet(
            f"border:1px solid #888; border-radius:3px; "
            f"background-color: rgb({color.red()}, {color.green()}, {color.blue()});"
        )

        # ✅ update underlying object immediately (optional but useful)
        basin_obj.color = (
            color.red()/255.0,
            color.green()/255.0,
            color.blue()/255.0
        )



    # ======================
    # ✅ APPLY
    # ======================
    def apply_changes(self):

        # tunnels
        for name, data in self.tunnel_widgets.items():

            obj = data["object"]

            obj["basins"]["visible"] = data["toggle"].isChecked()

            try:
                obj["basins"]["opacity"] = float(data["opacity"].text())
            except:
                pass

            obj["basins"]["color"] = self._get_button_color(data["color"])

        # basins
        for _, data in self.basin_widgets.items():

            b = data["object"]

            b.visible = data["checkbox"].isChecked()

            try:
                b.opacity = float(data["opacity"].text())
            except:
                pass

            b.color = self._get_button_color(data["color"])
        
        self.project.viewer.update_TuTraSt_plots(self.project)

        print("States updated!")

    

    def _histogram_basin(self, basin_id):

        mask = self.project.Basin_matrix == basin_id
        color = "steelblue"
        for b in self.project.Basin_list:
            if b.ID == basin_id:
                color = b.color

        plot_energy_histogram(
            self.project.Level_matrix,
            mask,
            self.project.grid, 
            self.project.E_levels["E_step"],
            f"Energy Histogram of Basin {basin_id}", 
            color = color
        )


    def _histogram_group(self, group_name):

        # -------------------------
        # Isolated clusters
        # -------------------------
        if group_name == "Isolated clusters":

            basin_ids = []

            for cluster in self.project.isolated_clusters.values():
                basin_ids.extend(cluster["basins"])

            color = self.project.isolated_clusters_plotting["basins"]["color"]

        # -------------------------
        # Tunnel systems
        # -------------------------
        else:

            basin_ids = self.project.tunnel_systems[group_name]["basins"]

            color = (
                self.project
                .tunnel_systems_plotting[group_name]
                ["basins"]["color"]
            )

        # -------------------------
        # Build mask
        # -------------------------
        mask = np.isin(
            self.project.Basin_matrix,
            basin_ids
        )

        # -------------------------
        # Plot
        # -------------------------
        plot_energy_histogram(
            self.project.Level_matrix,
            mask,
            self.project.grid,
            self.project.E_levels["E_step"],
            f"Energy Histogram of {group_name}",
            color=color
        )