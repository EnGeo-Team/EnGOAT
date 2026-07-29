# EnGOAT
<p align="center">

  <img src="logo.png" width="250">
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

## Features

- **PES-based** calculation of geometric descriptors, such as accessible volume and surface area  
- Identification and characterization of individual binding sites and transition states  
- Detection of diffusive channel networks with associated energy barriers  
- Efficient calculation of self-diffusion coefficients along crystallographic and Cartesian directions  
- Identification and visualization of minimum energy diffusion pathways

## Using EnGOAT

To use the software, download the EnGOAT package and run the installation script install.py to pip install all required python libraries. Then, simply run EnGOAT/EnGOAT.py to run the EnGOAT analysis software, or EnGOAT_visualization_tool/main.py to run the visualization application.

<!-- This is a comment 
This repository includes the **EnGOAT Instruction Manual**, which provides a comprehensive guide to using the software:

- **System Requirements & Installation** – Instructions to set up EnGOAT.
- **Input Files** – Detailed description of the files required to run the software.
- **Output Files** – Explanation of the generated results and metrics.
- **Visualization Tools** – Guide to using the tools in the `Tools` directory.
- **Appendix** – Detailed explanation of the underlying algorithm and calculation of metrics.

For hands-on guidance, the **Examples** directory contains tutorials showcasing potential use cases of EnGOAT. Each example includes:

- Required input files
- Generated output files
- A `.pdf` guide detailing:
  - How to run the software and tools
  - How to interpret the results

## Repository Overview

```
./
├── EnGOAT/                              # Directory containing the EnGOAT software package
│   └── ...
│
├── Examples/                            # Examples 
│   └── Li_ions_in_LGPS/
│       ├── Input_files
│       ├── Output_files
│       └── Li_ions_in_LGPS_tutorial.pdf
│
├── Tools/                               # Data visualization scripts
│   ├── MSD_plotter.py
│   └── Minimum_E_Pathway_plotter.py
│
├── EnGOAT_Instruction_Manual.pdf        # Instruction manual of the EnGOAT software
├── README.md                            # Project overview (this file)
├── LICENSE                              # License information
└── logo.png
```
-->
## Citation

If you use EnGOAT, please cite:
```
@article{TuTraSt,
  title={Automated multiscale approach to predict self-diffusion from a potential energy field},
  author={Mace, Amber and Barthel, Senja and Smit, Berend},
  journal={Journal of chemical theory and computation},
  volume={15},
  number={4},
  pages={2127--2141},
  year={2019},
  publisher={ACS Publications}
}
```
## License

This project is licensed under the MIT License.

## Authors

Matevž Turk,
Department of Chemistry - Ångström Laboratory; Structural Chemistry

Contact: matevz.turk@kemi.uu.se
