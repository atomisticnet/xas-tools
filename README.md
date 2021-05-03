# xas-tools

**Tools related to X-ray absorption spectroscopy (XAS) calculations**.

Part of DOE VTO Project: Haoyue Guo, Nong Artrith, Deyu Lu, Alex Urban

Contact: hg2568@columbia.edu, nartrith@atomistic.net, and a.urban@columbia.edu

## Installation

Install using `pip`

    pip install --user .

Or install in editable (developer) mode with

    pip install --user -e .

## Core-hole potential calculations with VASP

The following list describes the “workflow” that we need to implement:

1. Take PBE structure as input
2. Determine symmetrically distinct S atoms and their weights
3. (Potentially) create supercell
4. Generate VASP input files for
   * Single-point LDA calculation of the PBE structure (for band alignment); original cell
   * Core-hole potential calculations for all symmetrically distinct S atoms; supercell
5. After the calculations are done (XAS database):
   * Apply peak alignment to distinct S atoms: `Escf` = `Etotal_energy_of_excited_state` – `Etotal_energy_of_ground_state`
   * Average aligned outputs with the correct weights to compute the XAS spectrum

Most of this will be implemented in [`xas_tools.vasp`](./xas_tools/vasp.py).

## Acknowledgments

We acknowledge financial support by the U.S. Department of Energy (DOE) Office of Energy Efficiency and Renewable Energy, Vehicle Technologies Office, Contract No. DE-SC0012704. DFT calculations and machine-learning model construction made use of the Extreme Science and Engineering Discovery Environment (XSEDE), which is supported by National Science Foundation grant number ACI-1053575 (allocation no. DMR14005). Calculations were also performed on the computational resources of the Center for Functional Nanomaterials, which is a U.S. DOE Office of Science Facility, at Brookhaven National Laboratory under Contract No. DE-SC0012704. We also acknowledge computing resources from Columbia University’s Shared Research Computing Facility project, which is supported by NIH Research Facility Improvement Grant 1G20RR030893-01, and associated funds from the New York State Empire State Development, Division of Science Technology and Innovation (NYSTAR) Contract C090171, both awarded April 15, 2010.
