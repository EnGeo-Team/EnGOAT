from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parent

# List of Fortran modules to compile
modules = [
    ("find_neighbours", "find_neighbours.f90"),
    ("initiate_cluster", "initiate_cluster.f90"),
    ("get_process_cross_vector", "get_process_cross_vector.f90"),
]

for name, file in modules:
    subprocess.run(
        ["python3", "-m", "numpy.f2py", "-c", "-m", name, file],
        cwd=ROOT,
        check=True
    )
