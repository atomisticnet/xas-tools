# xas-tools

**Tools related to X-ray absorption spectroscopy (XAS) calculations**.

Haoyue Guo (<haoyue1619@gmail.com>), Nong Artrith
(<n.artrith@uu.nl>), Alex Urban (<a.urban@columbia.edu>)

See also: <br/>
Haoyue Guo, Matthew R. Carbone, Chuntian Cao, Jianzhou Qu, Yonghua Du, Seong-Min Bak, Conan Weiland, Feng Wang, Shinjae Yoo, Nongnuch Artrith, Alexander Urban, Deyu Lu, ["Simulated sulfur K-edge X-ray absorption spectroscopy database of lithium thiophosphate solid electrolytes", **2023**, https://doi.org/10.48550/arXiv.2302.00126](https://doi.org/10.48550/arXiv.2302.00126) 

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

## Reference

[1] F. Karsai, M. Humer, E. Flage-Larsen, P. Blaha and G. Kresse, *Phys. Rev. B* **98** (2018) 235205 .
