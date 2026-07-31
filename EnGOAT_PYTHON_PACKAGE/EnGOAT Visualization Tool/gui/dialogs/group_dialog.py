from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QCheckBox, QDialog, QDoubleSpinBox, QGridLayout, QGroupBox,
    QHBoxLayout, QLabel, QPushButton, QLayout, QScrollArea,
    QSizePolicy, QVBoxLayout, QWidget, QColorDialog
)

from plotting.matplotlib.plot_histogram import plot_energy_histogram


class GroupDialog(QDialog):
    COLUMN_WIDTHS = {0: 80, 1: 80, 2: 150, 3: 80, 4: 120, 5: 150, 6: 80, 7: 70, 8: 45}
    PANEL_MAX_HEIGHT = 300

    def __init__(self, project, group_type, group_id):
        super().__init__()
        self.project = project
        self.group_type = group_type
        self.group_id = str(group_id)

        title = "Tunnel system" if group_type == "tunnel_systems" else "Isolated group"
        self.setWindowTitle(f"{title} {group_id}")
        self.resize(1000, 700)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(0)

        self.build_group_header(layout, title)
        self.build_basin_panel(layout)
        self.build_ts_panel(layout)

        self.connect_master_signals(self.basin_widgets, "basins", group_type, group_id)
        for item_id, widgets in self.basin_widgets["individual_basins"].items():
            self.connect_item_signals(widgets, "basins", item_id)

        if self.TS_widgets["individual_TS"]:
            self.connect_master_signals(self.TS_widgets, "TS", group_type, group_id)
            for item_id, widgets in self.TS_widgets["individual_TS"].items():
                self.connect_item_signals(widgets, "TS", item_id)

        layout.addStretch(1)

    # ------------------------------------------------------------------
    # Data access
    # ------------------------------------------------------------------

    def get_group_data(self):
        key = "tunnel_systems" if self.group_type == "tunnel_systems" else "isolated_groups"
        return getattr(self.project, key)[self.group_id]

    def get_basin_ids(self):
        return self.get_group_data()["basin_list"]

    def get_ts_ids(self):
        return self.get_group_data()["TS_list"]

    def get_group_items(self, item_type):
        if item_type == "basins":
            return self.get_basin_ids(), self.basin_widgets["individual_basins"]
        return self.get_ts_ids(), self.TS_widgets["individual_TS"]

    # ------------------------------------------------------------------
    # Master controls
    # ------------------------------------------------------------------

    def on_master_visibility_changed(self, item_type, state):
        ids, widgets = self.get_group_items(item_type)
        key = "individual_basins" if item_type == "basins" else "individual_TS"
        visibility_dict = self.project.visibility[key]
        for ID in ids:
            visibility_dict[str(ID)] = state
        self.project.set_visibility(key, visibility_dict)

        for item_id in ids:
            widgets[str(item_id)]["visible"].setChecked(state)

    def on_master_color_clicked(self, item_type, group_type, group_ID):

        ids, widgets = self.get_group_items(item_type)
        if not ids:
            return

        current_color = self.project.plotting_data[group_type][group_ID]["color"]
        color = self.choose_color(current_color)

        if color == current_color:
            return

        plotting_data = self.project.plotting_data[item_type]
        
        for ID in ids:
            ID = str(ID)
            plotting_data[ID]
            plotting_data[ID]["color"] = color
            self.set_button_color(widgets[ID]["color"], color)

        self.project.set_plotting_data(item_type, plotting_data)

        master = self.basin_widgets["master_row"] if item_type == "basins" else self.TS_widgets["master_row"]
        self.set_button_color(master["color"], color)

    def on_master_opacity_changed(self, item_type, value):
        ids, widgets = self.get_group_items(item_type)
        plotting_data = self.project.plotting_data[item_type]

        for item_id in ids:
            item_id = str(item_id)
            spinbox = widgets[item_id]["opacity"]
            spinbox.blockSignals(True)
            spinbox.setValue(value)
            spinbox.blockSignals(False)

            plotting_data[str(item_id)]["opacity"] = value

        self.project.set_plotting_data(item_type, plotting_data)

    def on_master_label_changed(self, item_type, state):
        ids, widgets = self.get_group_items(item_type)
        key = "basin_labels" if item_type == "basins" else "TS_labels"
        visibility_dict = self.project.visibility[key]
        for ID in ids:
            visibility_dict[str(ID)] = state
        self.project.set_visibility(key, visibility_dict)

        for item_id in ids:
            widget = widgets[str(item_id)]["label"]
            widget.blockSignals(True)
            widget.setChecked(state)
            widget.blockSignals(False)
            

    def on_master_energy_changed(self, item_type, state):
        ids, widgets = self.get_group_items(item_type)
        key = "basin_energies" if item_type == "basins" else "TS_energies"
        visibility_dict = self.project.visibility[key]
        for ID in ids:
            visibility_dict[str(ID)] = state
        self.project.set_visibility(key, visibility_dict)

        for item_id in ids:
            widget = widgets[str(item_id)]["energy_toggle"]
            widget.blockSignals(True)
            widget.setChecked(state)
            widget.blockSignals(False)

    # ------------------------------------------------------------------
    # Signal connections
    # ------------------------------------------------------------------

    def connect_master_signals(self, widgets, item_type, group_type, group_id):
        master = widgets["master_row"]
        master["visible"].toggled.connect(lambda state: self.on_master_visibility_changed(item_type, state))
        master["color"].clicked.connect(lambda _, type = item_type, group = group_type, g_ID = group_id: self.on_master_color_clicked(type, group, g_ID))
        master["opacity"].valueChanged.connect(lambda value, type = item_type: self.on_master_opacity_changed(type, value))
        master["label"].toggled.connect(lambda state: self.on_master_label_changed(item_type, state))
        master["energy"].toggled.connect(lambda state: self.on_master_energy_changed(item_type, state))

    def connect_item_signals(self, widgets, item_type, item_id):
        widgets["visible"].toggled.connect(
            lambda state: self.on_visibility_changed(item_type, item_id, state)
        )
        widgets["color"].clicked.connect(
            lambda: self.on_color_changed(
                item_type,
                item_id,
                self.choose_color(self.project.plotting_data[item_type][item_id]["color"])
            )
        )
        widgets["opacity"].valueChanged.connect(
            lambda value: self.on_opacity_changed(item_type, item_id, value)
        )
        widgets["label"].toggled.connect(
            lambda state: self.on_label_changed(item_type, item_id, state)
        )
        widgets["energy_toggle"].toggled.connect(
            lambda state: self.on_energy_changed(item_type, item_id, state)
        )
        widgets["histogram"].clicked.connect(
            lambda: self.on_histogram_clicked(item_type, item_id)
        )

    # ------------------------------------------------------------------
    # Individual controls
    # ------------------------------------------------------------------

    def on_visibility_changed(self, item_type, item_id, state):
        key = "individual_basins" if item_type == "basins" else "individual_TS"
        visibility = self.project.visibility[key]
        visibility[str(item_id)] = state
        self.project.set_visibility(key, visibility)

    def on_color_changed(self, item_type, item_id, color):
        plotting_data = self.project.plotting_data[item_type]
        plotting_data[item_id]["color"] = color
        self.project.set_plotting_data(item_type, plotting_data)

        widgets = self.basin_widgets if item_type == "basins" else self.TS_widgets
        self.set_button_color(widgets[f"individual_{item_type}"][str(item_id)]["color"], color)

    def on_opacity_changed(self, item_type, item_id, value):
        plotting_data = self.project.plotting_data[item_type]
        plotting_data[item_id]["opacity"] = value
        self.project.set_plotting_data(item_type, plotting_data)

    def on_label_changed(self, item_type, item_id, state):
        key = "basin_labels" if item_type == "basins" else "TS_labels"
        visibility = self.project.visibility[key]
        visibility[item_id] = state
        self.project.set_visibility(key, visibility)

    def on_energy_changed(self, item_type, item_id, state):
        key = "basin_energies" if item_type == "basins" else "TS_energies"
        visibility = self.project.visibility[key]
        visibility[item_id] = state
        self.project.set_visibility(key, visibility)

    def on_histogram_clicked(self, item_type, item_id):
        if item_type == "basins":
            data = self.project.basin_data[item_id]
            title = f"Energy-Volume Histogram of Basin {item_id}"
        else:
            data = self.project.TS_data[item_id]
            title = f"Energy-Volume Histogram of Transition State {item_id}"

        color = self.project.plotting_data[item_type][item_id]["color"]
        plot_energy_histogram(data["histogram"], title, color)

    # ------------------------------------------------------------------
    # Group header
    # ------------------------------------------------------------------

    def build_group_header(self, parent_layout, title_word):
        group_data = self.get_group_data()

        header = QWidget()
        header.setSizePolicy(QSizePolicy.Maximum, QSizePolicy.Maximum)
        header_layout = QVBoxLayout(header)
        header_layout.setContentsMargins(4, 4, 4, 8)
        header_layout.setSpacing(6)

        title_widget = QWidget()
        title_layout = QHBoxLayout(title_widget)
        title_layout.setContentsMargins(0, 0, 0, 0)
        title_layout.setSpacing(4)

        title = QLabel(f"{title_word} {self.group_id}")
        font = title.font()
        font.setBold(True)
        font.setPointSize(font.pointSize() + 2)
        title.setFont(font)

        histogram_button = self.create_histogram_button()
        group_type = self.group_type
        histogram = group_data["histogram"]
        hist_title = f"Energy-Volume Histogram of {title_word} {self.group_id}"
        color = self.project.plotting_data[group_type][self.group_id]["color"]
        histogram_button.clicked.connect(
            lambda: plot_energy_histogram(histogram, hist_title, color)
        )

        title_layout.addWidget(title)
        title_layout.addWidget(histogram_button)

        values = [
            ("Volume", f"{group_data['V']:.2f} Å³ ({group_data['V_rel'] * 100:.2f} %)"),
            ("Minimum E", f"{group_data['E_min']:.2f} kJ/mol"),
            ("Area", f"{group_data['A']:.2f} Å²"),
        ]

        if self.group_type == "tunnel_systems":
            values.append(("Dimensionality", str(group_data["dimensionality"])))

        values_layout = QGridLayout()
        values_layout.setContentsMargins(0, 0, 0, 0)
        values_layout.setHorizontalSpacing(18)
        values_layout.setVerticalSpacing(4)

        for index, (label_text, value_text) in enumerate(values):
            row, column = divmod(index, 2)
            label = QLabel(f"{label_text}:")
            self.make_bold(label)
            values_layout.addWidget(label, row, column * 2)
            values_layout.addWidget(QLabel(value_text), row, column * 2 + 1)

        header_layout.addWidget(title_widget, 0, Qt.AlignCenter)
        header_layout.addLayout(values_layout)
        parent_layout.addWidget(header, 0, Qt.AlignCenter)

    # ------------------------------------------------------------------
    # Panel construction
    # ------------------------------------------------------------------

    def build_basin_panel(self, parent_layout):
        self.basin_widgets = self.build_panel(
            parent_layout,
            "basins",
            "Basins",
            self.get_basin_ids(),
            self.project.basin_data,
            "individual_basins",
            "basins",
            "Frac. Coordinates",
            self.build_basin_row_data,
        )

    def build_ts_panel(self, parent_layout):
        self.TS_widgets = self.build_panel(
            parent_layout,
            "TS",
            "Transition States",
            self.get_ts_ids(),
            self.project.TS_data,
            "individual_TS",
            "TS",
            "Connecting basins",
            self.build_ts_row_data,
        )

    def build_panel(self, parent_layout, panel_type, title, ids, data, visibility_key, plotting_key, master_header, row_data_builder):
        widgets = {"master_row": {}, f"individual_{panel_type}": {}}

        group_box = QGroupBox()
        group_layout = QVBoxLayout(group_box)
        group_layout.setContentsMargins(6, 6, 6, 6)
        group_layout.setSpacing(4)
        group_layout.setSizeConstraint(QLayout.SetMinimumSize)
        group_box.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Maximum)

        title_label = QLabel(title)
        title_label.setStyleSheet("font-weight: bold;")
        group_layout.addWidget(title_label)

        if ids or panel_type == "basins":
            group_layout.addLayout(
                self.build_master_row(master_header, widgets, panel_type)
            )

        rows_layout = QGridLayout()
        rows_layout.setContentsMargins(8, 0, 8, 4)
        rows_layout.setHorizontalSpacing(2)
        rows_layout.setVerticalSpacing(3)

        for row, item_id in enumerate(ids):
            item_id = str(item_id)
            row_widgets = row_data_builder(
                item_id,
                data[item_id],
                visibility_key,
                plotting_key
            )
            self.add_row_to_layout(rows_layout, row, row_widgets)
            widgets[f"individual_{panel_type}"][item_id] = row_widgets

        self.set_column_widths(rows_layout, self.COLUMN_WIDTHS)

        content_widget = QWidget()
        content_widget.setLayout(rows_layout)
        group_layout.addWidget(self.create_scroll_area(content_widget, self.PANEL_MAX_HEIGHT))

        parent_layout.addWidget(group_box)
        return widgets

    # ------------------------------------------------------------------
    # Master row
    # ------------------------------------------------------------------

    def build_master_row(self, master_header, widgets, panel_type):
        layout = QGridLayout()
        layout.setContentsMargins(8, 2, 8, 2)
        layout.setHorizontalSpacing(2)
        layout.setVerticalSpacing(0)

        show_all = QCheckBox("Show:")
        show_all_labels = QCheckBox("Label")
        show_all_energy = QCheckBox("E")

        for checkbox in (show_all, show_all_labels, show_all_energy):
            self.make_bold(checkbox)

        if panel_type == "basins":
            ids = self.get_basin_ids()
            label_key, energy_key = "basin_labels", "basin_energies"
        else:
            ids = self.get_ts_ids()
            label_key, energy_key = "TS_labels", "TS_energies"

        ids = [str(item_id) for item_id in ids]

        show_all.setChecked(self.project.visibility[self.group_type][self.group_id])
        show_all_labels.setChecked(any(self.project.visibility[label_key][item_id] for item_id in ids))
        show_all_energy.setChecked(any(self.project.visibility[energy_key][item_id] for item_id in ids))

        plotting_data = self.project.plotting_data[self.group_type][self.group_id]
        color_button = self.create_color_button(plotting_data["color"])
        opacity_box = self.create_opacity_box(plotting_data["opacity"])

        layout.addWidget(show_all, 0, 0, Qt.AlignLeft)
        layout.addWidget(self.create_labeled_widget("<b>Color:</b>", color_button), 0, 1, Qt.AlignCenter)
        layout.addWidget(self.create_labeled_widget("<b>Opacity:</b>", opacity_box), 0, 2, Qt.AlignCenter)
        layout.addWidget(show_all_labels, 0, 3, Qt.AlignCenter)

        energy_header = QWidget()
        energy_layout = QHBoxLayout(energy_header)
        energy_layout.setContentsMargins(0, 0, 0, 0)
        energy_layout.setSpacing(2)
        energy_layout.addWidget(show_all_energy)

        energy_label = QLabel("E [kJ/mol]")
        energy_label.setAlignment(Qt.AlignCenter)
        energy_label.setStyleSheet("font-weight: bold;")
        energy_layout.addWidget(energy_label)

        layout.addWidget(energy_header, 0, 4, Qt.AlignCenter)

        for column, text in [(5, master_header), (6, "V [Å³]"), (7, "A [Å²]")]:
            label = QLabel(text)
            label.setAlignment(Qt.AlignCenter)
            label.setStyleSheet("font-weight: bold;")
            layout.addWidget(label, 0, column, Qt.AlignCenter)

        layout.addWidget(QWidget(), 0, 8)
        self.set_column_widths(layout, self.COLUMN_WIDTHS)

        widgets["master_row"] = {
            "visible": show_all,
            "color": color_button,
            "opacity": opacity_box,
            "label": show_all_labels,
            "energy": show_all_energy,
        }

        return layout

    # ------------------------------------------------------------------
    # Row data
    # ------------------------------------------------------------------

    def build_basin_row_data(self, item_id, data, visibility_key, plotting_key):
        plotting_data = self.project.plotting_data[plotting_key][item_id]

        visible_cb = QCheckBox(f"B{item_id}")
        visible_cb.setChecked(self.project.visibility[visibility_key][item_id])
        visible_cb.setStyleSheet("font-weight: bold;")

        center = data["center"]
        grid_shape = self.project.metadata["grid_shape"]
        center_label = QLabel(
            f"({center[0] / grid_shape[0]:.2f}, "
            f"{center[1] / grid_shape[1]:.2f}, "
            f"{center[2] / grid_shape[2]:.2f})"
        )

        return self.create_row_widgets(
            visible_cb,
            plotting_data,
            self.project.visibility["basin_labels"][item_id],
            self.project.visibility["basin_energies"][item_id],
            center_label,
            data,
        )

    def build_ts_row_data(self, item_id, data, visibility_key, plotting_key):
        plotting_data = self.project.plotting_data[plotting_key][item_id]

        visible_cb = QCheckBox(f"TS{item_id}")
        visible_cb.setChecked(self.project.visibility[visibility_key][item_id])
        visible_cb.setStyleSheet("font-weight: bold;")

        center_label = QLabel(f"{data['basins'][0]}, {data['basins'][1]}")

        return self.create_row_widgets(
            visible_cb,
            plotting_data,
            False,
            False,
            center_label,
            data,
        )

    def create_row_widgets(self, visible_cb, plotting_data, label_visible, energy_visible, center_label, data):
        label_cb = QCheckBox()
        label_cb.setChecked(label_visible)

        energy_cb = QCheckBox()
        energy_cb.setChecked(energy_visible)

        values = self.create_value_labels(
            data["E_min"],
            data["V"],
            data["V_rel"],
            data["A"],
        )

        return {
            "visible": visible_cb,
            "color": self.create_color_button(plotting_data["color"]),
            "opacity": self.create_opacity_box(plotting_data["opacity"]),
            "label": label_cb,
            "energy_toggle": energy_cb,
            "center": center_label,
            "energy": values["energy"],
            "volume": values["volume"],
            "area": values["area"],
            "histogram": self.create_histogram_button(),
        }

    # ------------------------------------------------------------------
    # Row layout
    # ------------------------------------------------------------------

    def add_row_to_layout(self, layout, row, widgets):
        widgets["center"].setAlignment(Qt.AlignCenter)

        for key in ("energy", "volume", "area"):
            widgets[key].setAlignment(Qt.AlignCenter)

        energy_widget = QWidget()
        energy_layout = QHBoxLayout(energy_widget)
        energy_layout.setContentsMargins(0, 0, 0, 0)
        energy_layout.setSpacing(2)
        energy_layout.addWidget(widgets["energy_toggle"])
        energy_layout.addWidget(widgets["energy"])

        row_widgets = [
            (widgets["visible"], 0, Qt.AlignLeft),
            (widgets["color"], 1, Qt.AlignCenter),
            (widgets["opacity"], 2, Qt.AlignCenter),
            (widgets["label"], 3, Qt.AlignCenter),
            (energy_widget, 4, Qt.AlignCenter),
            (widgets["center"], 5, Qt.AlignCenter),
            (widgets["volume"], 6, Qt.AlignCenter),
            (widgets["area"], 7, Qt.AlignCenter),
            (widgets["histogram"], 8, Qt.AlignCenter),
        ]

        for widget, column, alignment in row_widgets:
            layout.addWidget(widget, row, column, alignment)

    # ------------------------------------------------------------------
    # Widget helpers
    # ------------------------------------------------------------------

    @staticmethod
    def choose_color(current_color):
        color = QColorDialog.getColor()
        return color.getRgbF()[:3] if color.isValid() else current_color

    @staticmethod
    def create_labeled_widget(label, widget):
        container = QWidget()
        layout = QHBoxLayout(container)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(2)
        layout.setSizeConstraint(QLayout.SetFixedSize)
        layout.addWidget(QLabel(label))
        layout.addWidget(widget)
        return container

    def create_color_button(self, color):
        button = QPushButton()
        button.setFixedSize(25, 15)
        self.set_button_color(button, color)
        return button

    @staticmethod
    def create_opacity_box(opacity):
        box = QDoubleSpinBox()
        box.setRange(0.0, 1.0)
        box.setSingleStep(0.1)
        box.setDecimals(2)
        box.setValue(opacity)
        box.setFixedWidth(90)
        return box

    @staticmethod
    def create_value_labels(energy, volume, relative_volume, area):
        return {
            "energy": QLabel(f"{energy:.2f}"),
            "volume": QLabel(f"{volume:.2f} ({100 * relative_volume:.2f}%)"),
            "area": QLabel(f"{area:.2f}"),
        }

    @staticmethod
    def create_histogram_button():
        button = QPushButton("▃▇▂")
        button.setFixedWidth(45)
        return button

    # ------------------------------------------------------------------
    # Layout helpers
    # ------------------------------------------------------------------

    @staticmethod
    def set_column_widths(layout, column_widths):
        for column, width in column_widths.items():
            layout.setColumnMinimumWidth(column, width)
            layout.setColumnStretch(column, 0)

    @staticmethod
    def create_scroll_area(widget, max_height=300):
        scroll_area = QScrollArea()
        scroll_area.setWidget(widget)
        scroll_area.setWidgetResizable(True)
        scroll_area.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        scroll_area.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        scroll_area.setMaximumHeight(max_height)
        return scroll_area

    @staticmethod
    def make_bold(widget):
        font = widget.font()
        font.setBold(True)
        widget.setFont(font)

    @staticmethod
    def set_button_color(button, color):
        r, g, b = (int(channel * 255) for channel in color)
        button.setStyleSheet(
            f"border: 1px solid #888; border-radius: 3px; "
            f"background-color: rgb({r}, {g}, {b});"
        )