import importlib
import subprocess
import sys


# Mapping:
# Python import name -> pip package name

REQUIRED_PACKAGES = {
    "numpy": "numpy",
    "scipy": "scipy",
    "ase": "ase",
    "PySide6": "PySide6",
    "pyvista": "pyvista",
    "pyvistaqt": "pyvistaqt",
    "matplotlib": "matplotlib",
    "mendeleev": "mendeleev",
}


def install_if_missing():

    for import_name, package_name in REQUIRED_PACKAGES.items():

        try:

            importlib.import_module(import_name)

            print(
                f"{package_name} is already installed."
            )

        except ImportError:

            print(
                f"{package_name} is missing. Installing..."
            )

            subprocess.check_call(
                [
                    sys.executable,
                    "-m",
                    "pip",
                    "install",
                    package_name,
                ]
            )

            print(
                f"{package_name} installed successfully."
            )


if __name__ == "__main__":

    install_if_missing()

    print(
        "\nAll required packages are installed."
    )