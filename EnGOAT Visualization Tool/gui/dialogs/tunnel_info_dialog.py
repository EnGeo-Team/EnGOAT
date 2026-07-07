from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QLabel, QTableWidget, QTableWidgetItem
)
from PySide6.QtWidgets import QGridLayout

class TunnelInfoDialog(QDialog):

    def __init__(self, tunnel_name, basins, transitions, grid, tunnel_info):
        super().__init__()

        self.setWindowTitle(tunnel_name)
        self.resize(600, 500)

        layout = QVBoxLayout(self)

        Na, Nb, Nc = grid["grid_points"]

        info = tunnel_info[tunnel_name]

        # ======================
        # ✅ HEADER INFO
        # ======================

        header_layout = QGridLayout()


        title = QLabel(tunnel_name)
        title.setStyleSheet("font-size: 14px; font-weight: bold;")
        layout.insertWidget(0, title)


        label = QLabel("Volume:")
        label.setStyleSheet("font-weight: bold;")
        header_layout.addWidget(label, 0, 0)
        
        header_layout.addWidget(QLabel(f"{info['V']:.2f} Å³"), 0, 1)
        
        
        label = QLabel("Surface area:")
        label.setStyleSheet("font-weight: bold;")
        header_layout.addWidget(label, 1, 0)
        
        header_layout.addWidget(QLabel(f"{info['A']:.2f} Å²"), 1, 1)
        
        
        label = QLabel("Dimensionality:")
        label.setStyleSheet("font-weight: bold;")
        header_layout.addWidget(label, 2, 0)
        
        header_layout.addWidget(QLabel(str(info["Dim"])), 2, 1)
        
        
        label = QLabel("Minimum E point:")
        label.setStyleSheet("font-weight: bold;")
        header_layout.addWidget(label, 3, 0)
        
        header_layout.addWidget(QLabel(f"{info['E_min']:.3f} kJ/mol"), 3, 1)

        layout.addLayout(header_layout)
        layout.addSpacing(10)

        # ======================
        # ✅ BASINS TABLE
        # ======================
        layout.addWidget(QLabel("Basins"))

        basin_table = QTableWidget()
        basin_table.setColumnCount(5)
        basin_table.setHorizontalHeaderLabels(
            ["Label", "Center", "E [kJ/mol]", "V [Å^3]", "A [Å²]"]
        )

        basin_table.setRowCount(len(basins))

        for row, b in enumerate(basins):
            basin_table.setItem(row, 0, QTableWidgetItem(f"B{str(b.ID)}"))
            basin_table.setItem(row, 1, QTableWidgetItem(f"({b.center[0]/Na:.2f}, {b.center[1]/Nb:.2f}, {b.center[2]/Nc:.2f})"))
            basin_table.setItem(row, 2, QTableWidgetItem(f"{b.E:.3f}"))
            basin_table.setItem(row, 3, QTableWidgetItem(f"{b.V:.3f}"))
            basin_table.setItem(row, 4, QTableWidgetItem(f"{b.A:.3f}"))

        basin_table.setEditTriggers(QTableWidget.NoEditTriggers)
        layout.addWidget(basin_table)

        # ======================
        # ✅ TRANSITIONS TABLE
        # ======================
        layout.addWidget(QLabel("Transition states"))

        ts_table = QTableWidget()
        ts_table.setColumnCount(5)
        ts_table.setHorizontalHeaderLabels(
            ["Label", "Basin 1", "Basin 2", "E [kJ/mol]", "Cross vector"]
        )

        ts_table.setRowCount(len(transitions))

        for row, ts in enumerate(transitions):
            ts_table.setItem(row, 0, QTableWidgetItem(f"TS{str(ts.ID)}"))
            ts_table.setItem(row, 1, QTableWidgetItem(str(ts.B_start)))
            ts_table.setItem(row, 2, QTableWidgetItem(str(ts.B_end)))
            ts_table.setItem(row, 3, QTableWidgetItem(f"{ts.E:.3f}"))
            ts_table.setItem(row, 4, QTableWidgetItem(str(ts.cross_vector)))

        ts_table.setEditTriggers(QTableWidget.NoEditTriggers)
        layout.addWidget(ts_table)
