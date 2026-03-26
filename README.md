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

## Getting Started

For detailed information on system requirements and installation, see Section 2 of the _EnGOAT Instruction Manual_.

For a step-by-step guide on running the software, along with practical considerations and guidance on result analysis, refer to the tutorials provided in the _Examples_ directory.

## Citation

If you use EnGOAT, please cite:

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

## License

This project is licensed under the MIT License.

## Authors

Matevž Turk
Department of Chemistry - Ångström Laboratory; Structural Chemistry

Contact: matevz.turk@kemi.uu.se
