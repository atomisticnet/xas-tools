# coding: utf-8
# Copyright (c) The Atomistic.net Team
# Distributed under the terms of the Mozilla
# Mozilla Public License, version 2.0
# (https://www.mozilla.org/en-US/MPL/2.0/)
"""
Object representations of spectroscopy data.

"""

import json
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from .util import broaden

__author__ = "Alex Urban, Nong Artrith"
__email__ = "a.urban@columbia.edu"
__maintainer__ = "Nong Artrith, Alexander Urban"
__maintainer_email__ = ("nartrith@atomistic.net"
                        + ", a.urban@columbia.edu")
__date__ = "2021-05-20"
__version__ = "0.1"


class AbsorptionSpectrum(object):

    def __init__(self, path):
        self.path = path
        with open(os.path.join(path, 'metadata.json')) as fp:
            self.metadata = json.load(fp)
        self.multiplicity = self.metadata['multiplicity']
        self.raw_data = []
        for csv_path in self.metadata['spectrum_path']:
            self.raw_data.append(pd.read_csv(os.path.join(path, csv_path)))
        self.E_min, self.E_max = self._get_edge_energies()
        self.broadened = None

    def __str__(self):
        out = 'Absorption Spectrum\n'
        out += '  data path: {}\n'.format(self.path)
        return out

    def _get_edge_energies(self):
        E_min = self.raw_data[0]['Aligned Energy (eV)'].values[-1]
        E_max = self.raw_data[0]['Aligned Energy (eV)'].values[0]
        for line in self.raw_data:
            idx = line['Intensity'].values > 0
            E = line['Aligned Energy (eV)'].values[idx]
            E_min = min(E_min, np.min(E))
            E_max = max(E_max, np.max(E))
        return E_min, E_max

    @property
    def spectrum(self):
        if self.broadened is None:
            return None
        else:
            X, Y = self.broadened[0]
            for b in self.broadened[1:]:
                Y += b[1]
            return X, Y

    def calculate_broadened(self, gauss_fwhm=None, lorentz_fwhm1=None,
                            lorentz_fwhm2=None, n=100, energy_range=None,
                            dE=1.0):
        """
        Calculate broadened atomic spectra and store them in the attribute
        `self.broadened`.

        """

        self.broadened = []
        E_range = [self.E_min - dE, self.E_max + dE]
        if energy_range is not None:
            if energy_range[0] is not None:
                # E_range[0] = max(E_range[0], energy_range[0])
                E_range[0] = energy_range[0]
            if energy_range[1] is not None:
                # E_range[1] = min(E_range[1], energy_range[1])
                E_range[1] = energy_range[1]
        for i, line in enumerate(self.raw_data):
            m = self.multiplicity[i]
            X = line['Aligned Energy (eV)'].values
            Y = line['Intensity'].values*m
            if gauss_fwhm is not None:
                X2, Y2 = broaden(X, Y, gauss_fwhm, xlim=E_range, n=n)
            else:
                X2, Y2 = X, Y
            if lorentz_fwhm1 is not None:
                if lorentz_fwhm2 is None:
                    raise ValueError("Both Lorentz parametera required.")
                else:
                    X, Y = broaden(X2, Y2, lorentz_fwhm1,
                                   fwhm2=lorentz_fwhm2, xlim=E_range,
                                   n=n, lorentz=True)
            else:
                X, Y = X2, Y2
            self.broadened.append((X, Y))

    def plot_atomic_lines(self, dE=1.0, **kwargs):
        """
        Plot the XAS lines for all individual atoms, weighted with their
        multiplicity.

        Args:
          dE (float): energy margin below and above lowest and highest
            energy, respectively
          **kwargs: Additional keyword arguments are used to update
            matplotlib's rcparams.

        Returns:
          fig, ax - matplotlib figure and axis

        """
        rcparams = {"font.size": 14,
                    "legend.frameon": False,
                    "xtick.top": True,
                    "xtick.direction": "in",
                    "xtick.minor.visible": True,
                    "xtick.major.size": 8,
                    "xtick.minor.size": 4,
                    "ytick.right": True,
                    "ytick.direction": "in",
                    "ytick.minor.visible": True,
                    "ytick.major.size": 8,
                    "ytick.minor.size": 4}
        if kwargs:
            rcparams.update(kwargs)
        plt.rcParams.update(rcparams)
        fig, ax = plt.subplots()
        for i, line in enumerate(self.raw_data):
            m = self.multiplicity[i]
            X = line['Aligned Energy (eV)'].values
            Y = line['Intensity']*m
            ax.plot(X, Y)
        ax.set_xlim([self.E_min - dE, self.E_max + dE])
        ax.set_yticks([])
        ax.set_xlabel("Energy (eV)")
        ax.set_ylabel("Intensity")
        return fig, ax

    def plot_broadened(self, dE=1.0, **kwargs):
        """
        Plot the XAS line from the combination of all atoms after
        broadening.  Broadened lines must first be calculated with the
        method `self.calculated_broadened`.

        Args:
          dE (float): energy margin below and above lowest and highest
            energy, respectively
          **kwargs: Additional keyword arguments are used to update
            matplotlib's rcparams.

        Returns:
          fig, ax - matplotlib figure and axis

        """
        if self.broadened is None:
            print("First call `calculate_broadened()`")
            return

        rcparams = {"font.size": 14,
                    "legend.frameon": False,
                    "xtick.top": True,
                    "xtick.direction": "in",
                    "xtick.minor.visible": True,
                    "xtick.major.size": 8,
                    "xtick.minor.size": 4,
                    "ytick.right": True,
                    "ytick.direction": "in",
                    "ytick.minor.visible": True,
                    "ytick.major.size": 8,
                    "ytick.minor.size": 4}
        if kwargs:
            rcparams.update(kwargs)
        plt.rcParams.update(rcparams)
        fig, ax = plt.subplots()
        X, Y = self.spectrum
        ax.plot(X, Y)
        ax.set_yticks([])
        ax.set_xlabel("Energy (eV)")
        ax.set_ylabel("Intensity")
        return fig, ax
