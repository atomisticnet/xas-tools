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

from .util import broaden, variable_convolution, gauss_function, lorentz_function

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
        self.has_xyz_components = False
        self.xyz_components = None

    def __str__(self):
        out = 'Absorption Spectrum\n'
        out += '  data path: {}\n'.format(self.path)
        return out

    def _get_edge_energies(self):
        E_min = self.raw_data[0]['Energy (eV)'].values[-1]
        E_max = self.raw_data[0]['Energy (eV)'].values[0]
        for line in self.raw_data:
            idx = line['Intensity'].values > 0
            E = line['Energy (eV)'].values[idx]
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

    @property
    def xyz_broadened(self):
        if self.xyz_components is None:
            return None
        else:
            X, Yx, Yy, Yz = self.xyz_components[0]
            for b in self.xyz_components[1:]:
                Yx += b[1]
                Yy += b[2]
                Yz += b[3]
            return X, Yx, Yy, Yz

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
                E_range[0] = energy_range[0]
            if energy_range[1] is not None:
                E_range[1] = energy_range[1]
        for i, line in enumerate(self.raw_data):
            m = self.multiplicity[i]
            X = line['Energy (eV)'].values
            Y = line['Intensity'].values*m
            try:
                omega_x = line['omega_x'].values*m
                omega_y = line['omega_y'].values*m
                omega_z = line['omega_z'].values*m
                self.xyz_components = []
                self.has_xyz_components = True
            except KeyError:
                self.has_xyz_components = False
            if gauss_fwhm is not None:
                sigma = gauss_fwhm/(2.0*np.sqrt(2.0*np.log(2.0)))
                X2, Y2 = variable_convolution(X, Y, gauss_function,
                                              sigma, num_points=n,
                                              x_range=E_range)
                if self.has_xyz_components:
                    _, Y2x = variable_convolution(
                        X, omega_x, gauss_function, sigma, num_points=n,
                        x_range=E_range)
                    _, Y2y = variable_convolution(
                        X, omega_y, gauss_function, sigma, num_points=n,
                        x_range=E_range)
                    _, Y2z = variable_convolution(
                        X, omega_z, gauss_function, sigma, num_points=n,
                        x_range=E_range)
            else:
                X2, Y2 = X, Y
                if self.has_xyz_components:
                    Y2x = omega_x
                    Y2y = omega_y
                    Y2z = omega_z
            if lorentz_fwhm1 is not None:
                if lorentz_fwhm2 is None:
                    raise ValueError("Both Lorentz parametera required.")
                else:
                    X, Y = variable_convolution(
                        X2, Y2, lorentz_function, lorentz_fwhm1,
                        sigma_end=lorentz_fwhm2, num_points=n,
                        x_range=E_range)
                    if self.has_xyz_components:
                        _, Yx = variable_convolution(
                            X2, Y2x, lorentz_function, lorentz_fwhm1,
                            sigma_end=lorentz_fwhm2, num_points=n,
                            x_range=E_range)
                        _, Yy = variable_convolution(
                            X2, Y2y, lorentz_function, lorentz_fwhm1,
                            sigma_end=lorentz_fwhm2, num_points=n,
                            x_range=E_range)
                        _, Yz = variable_convolution(
                            X2, Y2z, lorentz_function, lorentz_fwhm1,
                            sigma_end=lorentz_fwhm2, num_points=n,
                            x_range=E_range)
            else:
                X, Y = X2, Y2
                if self.has_xyz_components:
                    Yx = Y2x
                    Yy = Y2y
                    Yz = Y2z
            self.broadened.append((X, Y))
            if self.has_xyz_components:
                self.xyz_components.append((X, Yx, Yy, Yz))

    def calculate_broadened_old(self, gauss_fwhm=None, lorentz_fwhm1=None,
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
            X = line['Energy (eV)'].values
            Y = line['Intensity'].values*m
            try:
                omega_x = line['omega_x'].values*m
                omega_y = line['omega_y'].values*m
                omega_z = line['omega_z'].values*m
                self.xyz_components = []
                self.has_xyz_components = True
            except KeyError:
                self.has_xyz_components = False
            if gauss_fwhm is not None:
                X2, Y2 = broaden(X, Y, gauss_fwhm, xlim=E_range, n=n)
                if self.has_xyz_components:
                    _, Y2x = broaden(X, omega_x, gauss_fwhm, xlim=E_range, n=n)
                    _, Y2y = broaden(X, omega_y, gauss_fwhm, xlim=E_range, n=n)
                    _, Y2z = broaden(X, omega_z, gauss_fwhm, xlim=E_range, n=n)
            else:
                X2, Y2 = X, Y
                if self.has_xyz_components:
                    Y2x = omega_x
                    Y2y = omega_y
                    Y2z = omega_z
            if lorentz_fwhm1 is not None:
                if lorentz_fwhm2 is None:
                    raise ValueError("Both Lorentz parametera required.")
                else:
                    X, Y = broaden(X2, Y2, lorentz_fwhm1,
                                   fwhm2=lorentz_fwhm2, xlim=E_range,
                                   n=n, lorentz=True)
                    if self.has_xyz_components:
                        _, Yx = broaden(X2, Y2x, lorentz_fwhm1,
                                        fwhm2=lorentz_fwhm2, xlim=E_range,
                                        n=n, lorentz=True)
                        _, Yy = broaden(X2, Y2y, lorentz_fwhm1,
                                        fwhm2=lorentz_fwhm2, xlim=E_range,
                                        n=n, lorentz=True)
                        _, Yz = broaden(X2, Y2z, lorentz_fwhm1,
                                        fwhm2=lorentz_fwhm2, xlim=E_range,
                                        n=n, lorentz=True)
            else:
                X, Y = X2, Y2
                if self.has_xyz_components:
                    Yx = Y2x
                    Yy = Y2y
                    Yz = Y2z
            self.broadened.append((X, Y))
            if self.has_xyz_components:
                self.xyz_components.append((X, Yx, Yy, Yz))

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
            X = line['Energy (eV)'].values
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
