from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QLabel, QTableWidget, QTableWidgetItem, QScrollArea, QWidget
)


class IsolatedClustersDialog(QDialog):

    def __init__(self, project, clusters):
        super().__init__()

        self.setWindowTitle("Isolated clusters")
        self.resize(700, 600)

        main_layout = QVBoxLayout(self)

        # ✅ Scroll area (important if many clusters)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)

        content = QWidget()
        layout = QVBoxLayout(content)

        scroll.setWidget(content)
        main_layout.addWidget(scroll)

        # ======================
        # ✅ Populate clusters
        # ======================
        for name, data in clusters.items():
            
            label = QLabel(name)
            label.setStyleSheet("font-weight: bold;")
            layout.addWidget(label)

            basin_ids = set(data["basins"])
            ts_ids = set(data["transitions"])

            # ✅ Get actual data from project
            basins = [b for b in project.Basin_list if b.ID in basin_ids]
            transitions = [ts for ts in project.TS_list if ts.ID in ts_ids]
            grid = project.grid
            Na, Nb, Nc = grid["grid_points"]

            # -------------------------
            # ✅ BASINS
            # -------------------------
            layout.addWidget(QLabel("Basins:"))

            basin_table = QTableWidget()
            basin_table.setColumnCount(5)
            basin_table.setHorizontalHeaderLabels(
                ["Label", "Center", "E [kJ/mol]", "V [Å³]", "A [Å²]"]
            )
            basin_table.setRowCount(len(basins))

            for row, b in enumerate(basins):
                basin_table.setItem(row, 0, QTableWidgetItem(f"B{b.ID}"))
                basin_table.setItem(row, 1, QTableWidgetItem(f"({b.center[0]/Na:.2f}, {b.center[1]/Nb:.2f}, {b.center[2]/Nc:.2f})"))
                basin_table.setItem(row, 2, QTableWidgetItem(f"{b.E:.3f}"))
                basin_table.setItem(row, 3, QTableWidgetItem(f"{b.V:.3f}"))
                basin_table.setItem(row, 4, QTableWidgetItem(f"{b.A:.3f}"))

            basin_table.setEditTriggers(QTableWidget.NoEditTriggers)
            layout.addWidget(basin_table)

            # -------------------------
            # ✅ TRANSITIONS (only if present)
            # -------------------------
            if len(transitions) > 0:
            
                layout.addWidget(QLabel("Transitions:"))

                ts_table = QTableWidget()
                ts_table.setColumnCount(5)
                ts_table.setHorizontalHeaderLabels(
                    ["Label", "Basin 1", "Basin 2", "E [kJ/mol]", "Cross vector"]
                )
                ts_table.setRowCount(len(transitions))

                for row, ts in enumerate(transitions):
                    ts_table.setItem(row, 0, QTableWidgetItem(f"TS{ts.ID}"))
                    ts_table.setItem(row, 1, QTableWidgetItem(f"B{str(ts.B_start)}"))
                    ts_table.setItem(row, 2, QTableWidgetItem(f"B{str(ts.B_end)}"))
                    ts_table.setItem(row, 3, QTableWidgetItem(f"{ts.E:.3f}"))
                    ts_table.setItem(row, 4, QTableWidgetItem(str(ts.cross_vector)))

                ts_table.setEditTriggers(QTableWidget.NoEditTriggers)
                layout.addWidget(ts_table)

            layout.addSpacing(15)