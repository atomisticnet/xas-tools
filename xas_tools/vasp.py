# coding: utf-8
# Copyright (c) The Atomistic.net Team
# Distributed under the terms of the Mozilla
# Mozilla Public License, version 2.0
# (https://www.mozilla.org/en-US/MPL/2.0/)
"""
Set up and analyze XAS calculations performed with VASP.

"""

import os

import yaml
import numpy as np

import pymatgen as mg
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp.sets import DictSet

__author__ = "Alex Urban, Haoyue Guo, Jianzhou Qu, Qian Wang, Nong Artrith"
__email__ = "a.urban@columbia.edu, hg2568@columbia.edu"
__maintainer__ = "Haoyue Guo, Alexander Urban"
__maintainer_email__ = "hg2568@columbia.edu, a.urban@columbia.edu"
__date__ = "2021-03-24"
__version__ = "0.1"

HERE = os.path.dirname(os.path.abspath(__file__))


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

    @property
    def active_atom_types(self):
        """
        Returns the different types of atoms that are relevant for the XAS
        calculation.

        """
        return set([s.properties["ch_select"] for s in self.structure]) - {0}

    @property
    def active_mult(self):
        """
        Returns the multiplicity of each type of active atoms, as determined
        from symmetry analysis.  (e.g., how many atoms of sulfur type 1
        are present?)

        """
        atom_types = [s.properties["ch_select"] for s in self.structure]
        return {t: atom_types.count(t) for t in self.active_atom_types}

    def write_vasp_input(self, supercell=None, band_multiple=1, path=None):
        """
        Write VASP input files to the selected path.

        Args:
          supercell: Tuple with multiples of the three lattice vectors
          band_multiple: Set NBANDS to band_multiple times the number of
            valence electrons
          path: base path name for input file directories

        """
        struc_super = self.structure.copy()
        if supercell is not None:
            struc_super.make_supercell(supercell)

        with open(os.path.join(HERE, "LDA-XAS.yaml")) as fp:
            config = yaml.load(fp, Loader=yaml.FullLoader)

        for t in self.active_atom_types:
            # get sites of all active atoms of the same type
            active = [i for i, s in enumerate(struc_super)
                      if s.properties["ch_select"] == t]

            # get a list of all sites in the structure with one of the
            # active atoms placed at the beginning and create a new
            # structure
            site_list = struc_super.sites.copy()
            site_list = [site_list[active[0]]] + site_list
            del site_list[active[0]+1]
            xas_struc = mg.Structure.from_sites(site_list)

            vaspset = DictSet(xas_struc, config,
                              sort_structure=False,
                              force_gamma=True)

            # number of valence electrons for each atomic species,
            # identified by the atomic number
            num_elec = {
                p.atomic_no: np.sum(
                    [n for _, _, n in p.electron_configuration])
                for p in vaspset.potcar
            }

            # total number of valence electrons
            num_valence = np.sum(
                [num_elec[s.specie.number] for s in vaspset.structure])

            # set the number of bands according to the user-defined
            # multiple
            vaspset.user_incar_settings['NBANDS'] = num_valence*band_multiple

            if path is None:
                path = "XAS_input"
            vaspset.write_input(
                "{}_{}_{}".format(path, t, self.active_mult[t]),
                make_dir_if_not_present=True)
