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

The structure of this repository is shown below:

```
.
├── EnGOAT/                              # Directory containing the EnGOAT software package
│   ├── EnGOAT.py
│   ├── kMC.py
│   ├── data.py
│   ├── install.py
│   ├── get_topological_descriptors.py
│   ├── find_new_clusters.py
│   ├── grow_clusters.py
│   ├── Organize_TuTraSt.py
│   ├── PBC_minimax.py
│   ├── Read_cube_data.py
│   ├── find_neighbours.f90
│   ├── get_process_cross_vector.f90
│   └── initiate_cluster.f90
│
├── Examples/                            # Examples 
│   └── Li_ions_in_LGPS/
│       ├── Input_files
│       ├── Output_files
│       └──Li_ions_in_LGPS_tutorial.pdf
│
├── Tools/                               # Data visualization scripts
│   ├── MSD_plotter.py
│   └── Minimum_E_Pathway_plotter.py
│
├── images/                              # Images for documentation
│   └── logo.png
│
├── EnGOAT_Instruction_Manual.pdf        # Instruction manual of the EnGOAT software
├── README.md                            # Project overview (this file)
└── LICENSE                              # License information
```

