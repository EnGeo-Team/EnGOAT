from PySide6.QtWidgets import (
    QDialog,
    QVBoxLayout,
    QLabel,
    QCheckBox,
    QPushButton
)
from plotting.matplotlib.plot_merge_tree import create_merge_trees


class MergeTreeSelectionDialog(QDialog):

    def __init__(self, project, parent=None):
        super().__init__(parent)

        self.project = project

        self.setWindowTitle("Merge tree selection")
        self.resize(350, 300)

        layout = QVBoxLayout(self)

        layout.addWidget(
            QLabel("Select tunnel systems to include:")
        )

        self.checkboxes = {}

        for name in project.tunnel_systems.keys():

            cb = QCheckBox(name)

            layout.addWidget(cb)

            self.checkboxes[name] = cb

        if project.isolated_clusters:

            self.iso_checkbox = QCheckBox(
                "Isolated clusters"
            )

            layout.addWidget(self.iso_checkbox)

        else:
            self.iso_checkbox = None

        layout.addStretch()

        self.plot_button = QPushButton(
            "Plot merge tree"
        )

        layout.addWidget(self.plot_button)


        self.plot_button.clicked.connect(
            self.plot_merge_tree
        )

    def plot_merge_tree(self):

        selected_tunnels = [

            name

            for name, cb in self.checkboxes.items()

            if cb.isChecked()
        ]

        include_iso = (
            self.iso_checkbox is not None
            and self.iso_checkbox.isChecked()
        )

        self.accept()

        create_merge_trees(
            self.project,
            selected_tunnels,
            include_iso
        )


        
