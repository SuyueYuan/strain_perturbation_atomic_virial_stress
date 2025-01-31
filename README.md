# strain_perturbation_atomic_virial_stress
This repository contains example scripts for the work "A Strain Perturbation Method for Atomic Stress Calculation with Machine-Learning Potentials."

### Workflow:
The strain perturbation (SP) method follows these steps: 
1. Compute the instantaneous atomic energy in the undeformed state.
2.  Apply predefined perturbation strains to the simulation box and remap atomic positions accordingly (affine transformation).
3.  Recompute atomic energy for the deformed state.
4.  Use the atomic energies from steps (1) and (3) to calculate atomic stress/virial using Eq. (7) in the revised manuscript. <br>

Steps (1)-(3) can be performed using most MD simulation packages (here, we provide LAMMPS scripts). Step (4) is a simple post-processing step, for which we provide a Python script.


### Files:
- lammps.in: Main LAMMPS script for performing strain perturbation and dumping per-atom energy and Voronoi volume. The example script uses the SNAP potential for Si simulation results, but the application can be extend to any classical/ML many-body potentials.
- purterb_strain.in: Subriutines read by lammps.in.
- diff_virial.py: Python script for post-processing LAMMPS dump files for finite difference calculations.
- parity_plots.ipynb: Jupyter notebook for generating parity plots comparing SP method results with those from the rigid pairwise virial formulation.


### Notes:
In principle, this workflow can be implemented as a single function within an MD simulation package. We are actively exploring this integration, and updates will be available in this repository.



### ACKNOWLEDGEMENT
This work was sponsored by the Office of Energy Efficiency and Renewable Energy, Vehicle Technologies Office and was performed under the auspices of the U.S. Department of Energy by Lawrence Livermore National Laboratory under Contract DE-AC52-07NA27344. Computing support for this work came from the Department of Energy's Office of Energy Efficiency and Renewable Energy located at the National Renewable Energy Laboratory and the Computing Grand Challenge program from Lawrence Livermore National Laboratory.
