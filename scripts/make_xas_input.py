#!/usr/bin/env python

"""
Generate VASP input files for XAS core-hole potential calculations
for a given input file in POSCAR format.

"""

import argparse
import os

import pymatgen as mg
from xas_tools.vasp import CHPCalculation

__author__ = "Alexander Urban"
__email__ = "aurban@atomistic.net"
__date__ = "2021-03-27"
__version__ = "0.11"


def make_input(poscar_files, supercell, band_multiple, element, level,
               core_hole, grid_points, ncore, vasp_set):

    for i, f in enumerate(poscar_files):
        struc = mg.core.Structure.from_file(f)
        chp = CHPCalculation(struc, element=element, n=level[0],
                             ell=level[1], z=core_hole)
        path = os.path.basename("XAS_" + f + "_{}".format(i))
        if vasp_set is not None:
            chp.write_vasp_input(supercell=supercell,
                                 path=path,
                                 band_multiple=band_multiple,
                                 grid_points=grid_points,
                                 ncore=ncore,
                                 vasp_set=vasp_set)
        else:
            chp.write_vasp_input(supercell=supercell,
                                 path=path,
                                 band_multiple=band_multiple,
                                 grid_points=grid_points,
                                 ncore=ncore)


def main():

    parser = argparse.ArgumentParser(
        description=__doc__+"\n{} {}".format(__date__, __author__),
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument(
        "poscar_files",
        help="Path to one or more atomic structure file(s) in POSCAR format.",
        nargs="+")

    parser.add_argument(
        "--supercell", "-s",
        help="Supercell multiples (default: None).",
        type=int,
        default=None,
        nargs=3)

    parser.add_argument(
        "--band-multiple", "-n",
        help="Number of bands per valence electron (default: 1).",
        type=int,
        default=1)

    parser.add_argument(
        "--element", "-e",
        help="Chemical symbol of the element that is probed (default: S)",
        default="S")

    parser.add_argument(
        "--level", "-l",
        help="Principal and angular momentum quantum numbers of the level "
             "that is probed.  For example, the K edge (excitation from "
             "the 1s orbitals) would be (1, 0).  Default: (1, 0)",
        type=int,
        default=[1, 0],
        nargs=2)

    parser.add_argument(
        "--core-hole", "-c",
        help="Fraction of the core hole in electrons (default: 1.0).",
        type=float,
        default=1.0)

    parser.add_argument(
        "--grid-points", "-g",
        help="Number of grid points for the representation of "
             "spectral lines (default: 40000).",
        type=int,
        default=40000)

    parser.add_argument(
        "--ncore",
        help="VASP's NCORE parameter, which is usually square root of "
             "the number of cores (default: 1).",
        type=int,
        default=1)

    parser.add_argument(
        "--vasp-set",
        help="Path to a YAML file with VASP input file parameters.",
        default=None)

    args = parser.parse_args()

    make_input(args.poscar_files, args.supercell, args.band_multiple,
               args.element, args.level, args.core_hole, args.grid_points,
               args.ncore, args.vasp_set)


if (__name__ == "__main__"):
    main()
