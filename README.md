# EnGOAT
<p align="center">

  <img src="images/logo.png" width="250">
</p>
<p align="center">
  <img src="https://img.shields.io/badge/license-MIT-green">
  <img src="https://img.shields.io/badge/python-3.10-blue">
  <img src="https://img.shields.io/badge/fortran-90%2B-orange">
</p>

<p align="left">
  <strong>EnGOAT</strong> (Energy-GeOmetry Analysis Toolkit) is an open-source Python–Fortran package for analyzing potential energy surfaces of small species in crystalline solids. It is based on the TuTraSt algorithm, but extends its capabilities with additional functionality, improved performance, and a more user-friendly interface. The toolkit is designed to enable the analysis of structural and energetic properties, minimum energy pathways, and diffusion processes in solid-state systems.
</p>

---

## Overview

The structure of this repository is shown in the tree diagram below:

./
├── EnGOAT/                              #Directory containing the EnGOAT software
│   ├── EnGOAT.py
│   ├── Organize_TuTraSt.py
│   ├── PBC_minimax.py
│   ├── Read_cube_data.py
│   ├── data.py
│   ├── find_neighbours.f90
│   ├── find_new_clusters.py
│   ├── get_process_cross_vector.f90
│   ├── get_topological_descriptors.py
│   ├── grow_clusters.py
│   ├── initiate_cluster.f90
│   ├── install.py
│   └── kMC.py
├── Examples/                            #Directory 
│   └── Li_ions_in_LGPS
├── EnGOAT_Instruction_Manual.pdf
├── Tools/
│   ├── MSD_plotter.py
│   └── Minimum_E_Pathway_plotter.py
└── images/
    └── logo.png
├── LICENSE
├── README.md

---
