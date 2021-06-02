# coding: utf-8
# Copyright (c) The Atomistic.net Team
# Distributed under the terms of the Mozilla
# Mozilla Public License, version 2.0
# (https://www.mozilla.org/en-US/MPL/2.0/)
"""
Set up and analyze XAS calculations performed with VASP.

"""

import glob
import os
import re

import json
import yaml
import numpy as np

import pandas as pd

import pymatgen as mg
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp.sets import DictSet
from pymatgen.io.vasp.outputs import Outcar, Oszicar

import monty

__author__ = "Alex Urban, Haoyue Guo, Jianzhou Qu, Qian Wang, Nong Artrith"
__email__ = "a.urban@columbia.edu, hg2568@columbia.edu"
__maintainer__ = "Haoyue Guo, Nong Artrith, Alexander Urban"
__maintainer_email__ = ("hg2568@columbia.edu, nartrith@atomistic.net"
                        + ", a.urban@columbia.edu")
__date__ = "2021-03-24"
__version__ = "0.1"

HERE = os.path.dirname(os.path.abspath(__file__))


class CHPCalculation(object):

    def __init__(self, structure: mg.core.Structure, element: str = "S",
                 n: int = 1, ell: int = 0, z: float = 1.0):
        """
        Args:
          structure: Atomic structure
          element: Chemical symbol of the element that is probed
          n: principal quantum number of the excited electron
          ell: angular momentum quantum number of the excited electron
          z: electron count, i.e., for fractional core holes is shielding
            is considered

          Example: the K edge would correspond to (n, ell) = (1, 0), i.e.,
                   excitation from the 1s orbitals.

        """
        self.structure = structure
        self.xas_element = mg.Element(element)
        self.ch_n = n
        self.ch_l = ell
        self.ch_z = z
        self.atoms = [(i, s) for i, s in enumerate(structure)
                      if s.specie == self.xas_element]
        self.equivalent_atoms = self._check_equivalent_atoms()
        self.weights = [len(group) for group in self.equivalent_atoms]

        # mark active atoms with the symmetry group that they belong to
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
        cnt = {t: atom_types.count(t) for t in self.active_atom_types}
        n = np.gcd.reduce(list(cnt.values()))
        return {t: cnt[t]//n for t in cnt}

    def write_vasp_input(self, supercell=None, band_multiple=1,
                         path=".", vasp_set="LDA-XAS", write_no_ch=True):
        """
        Write VASP input files to the selected path.

        Args:
          supercell: Tuple with multiples of the three lattice vectors
          band_multiple: Set NBANDS to band_multiple times the number of
            valence electrons
          path: base path name for input file directories
          vasp_set: name of the VASP input parameter set or path to a
            YAML file; currently only one named set is available (LDA-XAS)
          write_no_ch: If True, will also generate input files for a
            regular SCF calculations without core hole using the same
            input file parameters.

        """
        struc_super = self.structure.copy()
        if supercell is not None:
            struc_super.make_supercell(supercell)

        try:
            with open(os.path.join(HERE, "{}.yaml".format(vasp_set))) as fp:
                config = yaml.load(fp, Loader=yaml.FullLoader)
        except FileNotFoundError:
            with open(vasp_set) as fp:
                config = yaml.load(fp, Loader=yaml.FullLoader)

        vaspset = DictSet(struc_super, config,
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

        if write_no_ch:
            # set the number of bands according to the user-defined
            # multiple:
            vaspset.user_incar_settings.update({
                'NBANDS': num_valence*band_multiple
            })
            vaspset.write_input(
                os.path.join(path, "XAS_input_SCF"),
                make_dir_if_not_present=True)

        w1 = len(str(max(self.active_atom_types)))
        w2 = len(str(max(self.active_mult.values())))
        dir_frmt = ("XAS_input_{:0" + "{}".format(w1)
                    + "d}_{:0" + "{}".format(w2) + "d}")

        for t in self.active_atom_types:
            # Get sites of all active atoms of the same type
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

            # set the number of bands according to the user-defined
            # multiple and set core-hole potential parameters:
            vaspset.user_incar_settings.update({
                'NBANDS': num_valence*band_multiple,
                'ICORELEVEL': 2,
                'CLNT': 1,
                'CLN': self.ch_n,
                'CLL': self.ch_l,
                'CLZ': self.ch_z
            })

            vaspset.write_input(
                os.path.join(path, dir_frmt.format(t, self.active_mult[t])),
                make_dir_if_not_present=True)


def parse_vasp_chp_output(base_path, output_path='XAS_output'):
    """
    Read and process output from core-hole potential calculations.

    Args:
      base_path (str): Path to a directory with VASP output subdirectories.
        The directory needs to contain an scf calculation in a directory
        named "XAS_input_SCF" and the core-hole potential calculations
        in directories named "XAS_input_*_*".
      output_path (str): Directory in which the processed XAS output data
        should be collected.

    Returns:
      An AbsorptionSpectrum object

    """

    no_ch_path = glob.glob(os.path.join(base_path, 'XAS_input_SCF'))[0]
    ch_paths = glob.glob(os.path.join(base_path, 'XAS_input_*_*'))

    if not os.path.exists(output_path):
        os.makedirs(output_path)
        outdir = os.path.join(
            output_path, 'structure_{}'.format(1))
    else:
        num_old = len(glob.glob(os.path.join(output_path, 'structure_*')))
        outdir = os.path.join(
            output_path, 'structure_{}'.format(num_old + 1))
    os.mkdir(outdir)

    # read the total energy of the calculation without core hole
    oszicar_path_no_ch = glob.glob(os.path.join(no_ch_path, "OSZICAR"))[0]
    oszicar_no_ch = Oszicar(oszicar_path_no_ch)
    etot_no_ch = oszicar_no_ch.ionic_steps[-1]['E0']

    # read also the structures to determine the size of the supercell
    structure_no_ch = mg.Structure.from_file(
        os.path.join(no_ch_path, 'CONTCAR'))
    structure_ch = mg.Structure.from_file(os.path.join(ch_paths[0], 'CONTCAR'))
    N_no_ch = len(structure_no_ch.sites)
    N_ch = len(structure_ch.sites)
    supercell = N_ch/N_no_ch

    multiplicity = [int(os.path.basename(p).split('_')[-1]) for p in ch_paths]
    metadata = {
        'origin': os.path.relpath(base_path, outdir),
        'supercell_size': supercell,
        'supercell_composition': structure_ch.composition.formula,
        'unique_atoms': len(ch_paths),
        'multiplicity': multiplicity,
        'total_energy_scf': etot_no_ch*supercell,
        'total_energy_ch': [],
        'fermi_energy_ch': [],
        'scf_path': os.path.relpath(no_ch_path, outdir),
        'ch_path': [],
        'spectrum_path': []
    }

    w1 = len(str(len(ch_paths)))
    w2 = len(str(max(multiplicity)))
    csv_path_frmt = ("atom_{:0" + "{}".format(w1)
                     + "d}_{:0" + "{}".format(w2) + "d}.csv")

    for i, path in enumerate(ch_paths):
        outcar_path = glob.glob(os.path.join(path, "OUTCAR*"))[0]
        oszicar_path = glob.glob(os.path.join(path, "OSZICAR*"))[0]

        # read total energy
        oszicar = Oszicar(oszicar_path)
        etot_ch = oszicar.ionic_steps[-1]['E0']
        dE = etot_ch - etot_no_ch*supercell

        # read Fermi energy
        oc = Outcar(outcar_path)

        # extract the raw XAS spectrum
        xas = []
        with monty.io.zopen(outcar_path) as fp:
            line = fp.readline()
            while not re.search('IMAGINARY DIELECTRIC FUNCTION',
                                line.decode('utf-8')):
                line = fp.readline()
            fp.readline()
            fp.readline()
            line = fp.readline().decode('utf-8')
            while len(line.strip()) > 0:
                dielectric = [float(a) for a in line.split()]
                xas.append([dielectric[0], dielectric[0] + dE,
                            np.sum(dielectric[1:])])
                line = fp.readline().decode('utf-8')
        xas = np.array(xas)

        df = pd.DataFrame(data=xas,
                          columns=['Raw Energy (eV)',
                                   'Aligned Energy (eV)',
                                   'Intensity'])
        csv_path = csv_path_frmt.format(i+1, multiplicity[i])
        df.to_csv(os.path.join(outdir, csv_path), index=False)

        # keep track of meta data
        metadata['ch_path'].append(os.path.relpath(path, outdir))
        metadata['total_energy_ch'].append(etot_ch)
        metadata['fermi_energy_ch'].append(oc.efermi)
        metadata['spectrum_path'].append(csv_path)

    with open(os.path.join(outdir, 'metadata.json'), 'w') as fp:
        json.dump(metadata, fp)
