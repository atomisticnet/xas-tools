# xas-tools

**Tools related to X-ray absorption spectroscopy (XAS) calculations**.

Contact: Haoyue Guo, Nong Artrith, Deyu Lu, Alex Urban

Emails: hg2568@columbia.edu, and a.urban@columbia.edu

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
5. After the calculations are done:
   * Apply peak alignment to distinct S atoms: `Efermi + Escf`, `Escf` = `Etotal_energy_of_excited_state` – `Etotal_energy_of_ground_state`
   * Average aligned outputs with the correct weights to compute the XAS spectrum

Most of this will be implemented in [`xas_tools.vasp`](./xas_tools/vasp.py).
