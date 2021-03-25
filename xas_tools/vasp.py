# coding: utf-8
# Copyright (c) The Atomistic.net Team
# Distributed under the terms of the Mozilla
# Mozilla Public License, version 2.0
# (https://www.mozilla.org/en-US/MPL/2.0/)
"""
Set up and analyze XAS calculations performed with VASP.

"""

import numpy as np

import pymatgen as mg
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

__author__ = "Alex Urban, Haoyue Guo, Jianzhou Qu, Qian Wang, Nong Artrith"
__email__ = "a.urban@columbia.edu, hg2568@columbia.edu"
__maintainer__ = "Haoyue Guo, Alexander Urban"
__maintainer_email__ = "hg2568@columbia.edu, a.urban@columbia.edu"
__date__ = "2021-03-24"
__version__ = "0.1"


class CHPCalculation(object):

    def __init__(self, structure: mg.core.Structure, element: str = "S"):
        """
        Args:
          structure: Atomic structure
          species: Chemical symbol of the element that is probed

        """
        self.structure = structure
        self.xas_element = mg.Element(element)
        self.atoms = [(i, s) for i, s in enumerate(structure)
                      if s.specie == self.xas_element]
        self.equivalent_atoms = self._check_equivalent_atoms()
        self.weights = [len(group) for group in self.equivalent_atoms]

        # mark active atoms with the group that they belong to
        select_for_core_hole = np.array([0 for i in self.structure])
        for i, equiv_set in enumerate(self.equivalent_atoms):
            select_for_core_hole[equiv_set] = i + 1
        self.structure.add_site_property("ch_select", select_for_core_hole)

    def __str__(self):
        return

    def _check_equivalent_atoms(self):
        """
        Analyze symmetry to determine which atoms of the XAS species are
        equivalent.

        """
        sga = SpacegroupAnalyzer(self.structure)
        symmops = sga.get_symmetry_operations()
        equivalent_atoms = []
        assigned = []
        all_coords = [s.frac_coords for _, s in self.atoms]
        for i, atom in self.atoms:
            if i in assigned:
                continue
            coo1 = atom.frac_coords
            equivalent_atoms.append([i])
            assigned.append(i)
            for so in symmops[1:]:
                symm_coo1 = so.operate(coo1)
                distances = self.structure.lattice.get_all_distances(
                    [symm_coo1], all_coords)[0]
                equiv = [j for j in range(len(self.atoms))
                         if distances[j] <= 1.0e-3]
                for atom2 in equiv:
                    if self.atoms[atom2][0] not in assigned:
                        equivalent_atoms[-1].append(self.atoms[atom2][0])
                        assigned.append(self.atoms[atom2][0])
        return equivalent_atoms
