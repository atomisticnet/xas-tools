#!/usr/bin/env python

"""
Extract XAS spectra from VASP core-hole potential calculations.

"""

import argparse

from xas_tools.vasp import parse_vasp_chp_output

__author__ = "Alexander Urban"
__email__ = "aurban@atomistic.net"
__date__ = "2021-06-09"
__version__ = "0.11"


def analyze_output(paths):
    for p in paths:
        parse_vasp_chp_output(p)


def main():

    parser = argparse.ArgumentParser(
        description=__doc__+"\n{} {}".format(__date__, __author__),
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument(
        "paths",
        help="Paths to the base directories of VASP CHP calculations.",
        nargs="+")

    args = parser.parse_args()

    analyze_output(args.paths)


if (__name__ == "__main__"):
    main()
