# xas-tools

**Tools related to X-ray absorption spectroscopy (XAS) calculations**.

This package is developed as part of a project funded by the
U.S. Department of Energy.  See acknowledgment below.

Haoyue Guo (<hguo1@bnl.gov>), Matthew Carbone, Nong Artrith
(<nartrith@atomistic.net>), Deyu Lu, Alex Urban (<a.urban@columbia.edu>)

## Installation

Install using `pip`

    pip install --user .

Or install in editable (developer) mode with

    pip install --user -e .

## XAS calculations with VASP

The package generates input files for XAS simulations with the supercell
core-hole approach [1].

The following list describes the generak “workflow”:

1. Given an input atomic structure
2. All distinct absorbers are determined based on the symmetry of the
   structure
3. Multiplicities are determined based on the number of symmetrically
   equivalent absorbers
4. If requested, a supercell is created
5. VASP input files are generated for
   * A ground-state calculation without core hole
   * Core-hole potential calculations for all symmetrically distinct
     absorbers
6. VASP calculations are performed as usual
7. The VASP output is post-processed by
   * Alignment of the XAS lines for individual absorbers
   * Combining and broadening the XAS lines of all individual absorbers

## Acknowledgments

We acknowledge financial support by the U.S. Department of Energy (DOE)
Office of Energy Efficiency and Renewable Energy, Vehicle Technologies
Office, Contract No. DE-SC0012704. DFT calculations and machine-learning
model construction made use of the Extreme Science and Engineering
Discovery Environment (XSEDE), which is supported by National Science
Foundation grant number ACI-1053575 (allocation
no. DMR14005). Calculations were also performed on the computational
resources of the Center for Functional Nanomaterials, which is a
U.S. DOE Office of Science Facility, at Brookhaven National Laboratory
under Contract No. DE-SC0012704. We also acknowledge computing resources
from Columbia University’s Shared Research Computing Facility project,
which is supported by NIH Research Facility Improvement Grant
1G20RR030893-01, and associated funds from the New York State Empire
State Development, Division of Science Technology and Innovation
(NYSTAR) Contract C090171, both awarded April 15, 2010.

## Reference

[1] F. Karsai, M. Humer, E. Flage-Larsen, P. Blaha and G. Kresse, *Phys. Rev. B* **98** (2018) 235205 .
