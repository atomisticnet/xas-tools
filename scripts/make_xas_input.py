#!/usr/bin/env python

"""
Generate VASP input files for XAS core-hole potential calculations
for a given input file in POSCAR format.

"""

import argparse

import pymatgen as mg
from xas_tools.vasp import CHPCalculation

__author__ = "Alexander Urban"
__email__ = "aurban@atomistic.net"
__date__ = "2021-03-27"
__version__ = "0.1"


def make_input(poscar_file, supercell, band_multiple):
    struc = mg.Structure.from_file(poscar_file)
    chp = CHPCalculation(struc)
    chp.write_vasp_input(supercell=supercell, band_multiple=band_multiple)


if (__name__ == "__main__"):

    parser = argparse.ArgumentParser(
        description=__doc__+"\n{} {}".format(__date__, __author__),
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument(
        "poscar_file",
        help="Path to the poscar file.")

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

    args = parser.parse_args()

    make_input(args.poscar_file, args.supercell, args.band_multiple)
