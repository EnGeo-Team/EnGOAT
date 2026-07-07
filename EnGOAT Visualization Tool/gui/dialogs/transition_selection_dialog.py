from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QGroupBox, QHBoxLayout,
    QPushButton, QCheckBox, QLabel, QLineEdit, QWidget, QScrollArea, QColorDialog
)


class TransitionSelectionDialog(QDialog):

    def __init__(self, project):
        super().__init__()

        self.project = project

        self.setWindowTitle("Plot individual transitions")
        self.resize(600, 700)

        self.ts_widgets = {}
        self.tunnel_widgets = {}
        self.group_ts_widgets = {}

        main_layout = QVBoxLayout(self)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)

        content = QWidget()
        self.layout = QVBoxLayout(content)

        scroll.setWidget(content)
        main_layout.addWidget(scroll)

        self.apply_btn = QPushButton("Apply changes")
        self.apply_btn.setFixedHeight(35)
        main_layout.addWidget(self.apply_btn)

        self.apply_btn.clicked.connect(self.apply_changes)

        self.build_ui()


    # ======================
    # ✅ BUILD UI
    # ======================
    def build_ui(self):

        tunnel_groups = {}
        isolated_ts = []

        for ts in self.project.TS_list:
            if ts.tunnel:
                key = f"Tunnel system {ts.tunnel}"
                tunnel_groups.setdefault(key, []).append(ts)
            else:
                isolated_ts.append(ts)

        def extract_number(name):
            return int(name.split()[-1])

        sorted_keys = sorted(tunnel_groups.keys(), key=extract_number)

        for name in sorted_keys:
            self._add_group(name, tunnel_groups[name])

        if isolated_ts:
            self._add_group("Isolated clusters", isolated_ts)

        self.layout.addStretch()


    def _add_group(self, name, ts_list):

        group = QGroupBox(name)
        group_layout = QVBoxLayout()

        self.group_ts_widgets[name] = []

        # ======================
        # ✅ HEADER
        # ======================
        header_row = QHBoxLayout()

        if name == "Isolated clusters":
            info = self.project.isolated_clusters_plotting
        else:
            info = self.project.tunnel_systems_plotting.get(name)

        plot = info["TS"]

        toggle = QCheckBox(name)
        toggle.setChecked(plot["visible"])
        toggle.setStyleSheet("font-weight: bold;")

        color_btn = QPushButton()
        color_btn.setFixedSize(25, 15)
        self._set_button_color(color_btn, plot["color"])

        opacity = QLineEdit(f"{plot['opacity']:.2f}")
        opacity.setFixedWidth(50)

        self.tunnel_widgets[name] = {
            "toggle": toggle,
            "color": color_btn,
            "opacity": opacity,
            "object": plot
        }

        # ✅ group actions
        toggle.toggled.connect(lambda val, n=name: self._apply_group_toggle(n, val))
        color_btn.clicked.connect(lambda _, n=name, b=color_btn: self._apply_group_color(n, b))
        opacity.textChanged.connect(lambda t, n=name: self._apply_group_opacity(n, t))

        header_row.addWidget(toggle)
        header_row.addStretch()
        header_row.addWidget(QLabel("Color:"))
        header_row.addWidget(color_btn)
        header_row.addSpacing(10)
        header_row.addWidget(QLabel("Opacity:"))
        header_row.addWidget(opacity)

        group_layout.addLayout(header_row)

        # ======================
        # ✅ TRANSITIONS
        # ======================
        for ts in ts_list:

            row = QHBoxLayout()

            cb = QCheckBox(f"TS{ts.ID}")
            cb.setChecked(ts.visible)

            color_btn = QPushButton()
            color_btn.setFixedSize(25, 15)
            self._set_button_color(color_btn, ts.color)

            # ✅ single TS color
            color_btn.clicked.connect(
                lambda _, btn=color_btn, obj=ts: self._change_single_color(btn, obj)
            )

            opacity = QLineEdit(f"{ts.opacity:.2f}")
            opacity.setFixedWidth(50)

            row.addWidget(cb)
            row.addStretch()
            row.addWidget(QLabel("Color:"))
            row.addWidget(color_btn)
            row.addSpacing(10)
            row.addWidget(QLabel("Opacity:"))
            row.addWidget(opacity)

            group_layout.addLayout(row)

            obj = {
                "checkbox": cb,
                "color": color_btn,
                "opacity": opacity,
                "object": ts
            }

            self.ts_widgets[ts.ID] = obj
            self.group_ts_widgets[name].append(obj)

        group.setLayout(group_layout)
        self.layout.addWidget(group)


    # ======================
    # ✅ GROUP ACTIONS
    # ======================
    def _apply_group_toggle(self, name, state):
        for item in self.group_ts_widgets.get(name, []):
            item["checkbox"].setChecked(state)


    def _apply_group_color(self, name, button):

        color = QColorDialog.getColor()
        if not color.isValid():
            return

        rgb = (color.red(), color.green(), color.blue())

        self._set_button_color(button, [c/255 for c in rgb])

        for item in self.group_ts_widgets.get(name, []):
            self._set_button_color(item["color"], [c/255 for c in rgb])


    def _apply_group_opacity(self, name, text):
        for item in self.group_ts_widgets.get(name, []):
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


    def _change_single_color(self, button, obj):

        color = QColorDialog.getColor()
        if not color.isValid():
            return

        button.setStyleSheet(
            f"border:1px solid #888; border-radius:3px; "
            f"background-color: rgb({color.red()}, {color.green()}, {color.blue()});"
        )

        obj.color = (
            color.red()/255,
            color.green()/255,
            color.blue()/255
        )


    # ======================
    # ✅ APPLY
    # ======================
    def apply_changes(self):

        for _, data in self.ts_widgets.items():

            ts = data["object"]

            ts.visible = data["checkbox"].isChecked()

            try:
                ts.opacity = float(data["opacity"].text())
            except:
                pass

            ts.color = self._get_button_color(data["color"])
        
        self.project.viewer.update_TuTraSt_plots(self.project)

        print("Transition states updated!")