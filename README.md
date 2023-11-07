# MORIS

Multi-physics Optimization Research and Innovation System (2023)

## Overview

MORIS is a research code developed at the University of Colorado in the Aerospace Mechanics Research Center (AMRec)* by Kurt Maute's research group. AMRec was started in the 1990's under the name "Center for Aerospace Structures". It is housed in the Ann and H.J. Smead Department of Aerospace Engineering Sciences.

This code provides a research platform for PDE constrained optimization with focus on shape and topology optimization using an isogeometric formulation of the extended finite element method.

## Documentation

### Installation 
A script to install instructions can be found in [InstallScript](https://github.com/kkmaute/moris/tree/main/share/install/IntallScript.sh). Edit the script and verify/update a small number of installation parameters. Should the script fail, detailed install instructions can be found at [Install_Step_By_Step.txt](https://github.com/kkmaute/moris/tree/main/share/install/Install_Step_By_Step.txt)

### Mesh Generation for Interpolation Based Finite Element Analysis
Documentation for generating and outputting meshes and extraction operators as defined in this [paper](https://doi.org/10.1016/j.cma.2023.115890) can be found [here](https://github.com/kkmaute/moris/blob/main/share/doc/mesh_generation/main.pdf). GitHub's integrated PDF viewer doesn't support links, so please download the pdf to use them. Example input files can be found in [share/doc/mesh_generation/examples](https://github.com/kkmaute/moris/tree/main/share/doc/mesh_generation/examples).
