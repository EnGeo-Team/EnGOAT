from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout,
    QCheckBox, QPushButton, QLabel, QScrollArea, QWidget
)
from PySide6.QtCore import Qt


class BondSelectionDialog(QDialog):

    def __init__(self, bonds, selected_bonds):
        super().__init__()

        self.setWindowTitle("Choose bonds")

        # ✅ fixed size
        self.setFixedSize(250, 350)

        main_layout = QVBoxLayout(self)

        # ✅ title label
        title = QLabel("Choose bonds:")
        title.setAlignment(Qt.AlignLeft)
        title.setStyleSheet("font-weight: bold;")
        main_layout.addWidget(title)

        # ✅ scroll area (important if many elements)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)

        container = QWidget()
        scroll_layout = QVBoxLayout(container)

        self.checkboxes = {}

        for bl in sorted(bonds):
            label = f"{bl[0]} – {bl[1]}"
            cb = QCheckBox(label)
            cb.setChecked(bl in selected_bonds)
            scroll_layout.addWidget(cb)
            self.checkboxes[bl] = cb

        scroll_layout.addStretch()
        scroll.setWidget(container)

        main_layout.addWidget(scroll)

        # ✅ buttons (aligned horizontally)
        button_layout = QHBoxLayout()

        btn_apply = QPushButton("Apply")
        btn_cancel = QPushButton("Cancel")

        button_layout.addStretch()
        button_layout.addWidget(btn_apply)
        button_layout.addWidget(btn_cancel)

        main_layout.addLayout(button_layout)

        # ✅ connections
        btn_apply.clicked.connect(self.accept)
        btn_cancel.clicked.connect(self.reject)

    def get_selection(self):
        return [
            bl for bl, cb in self.checkboxes.items() if cb.isChecked()
        ]