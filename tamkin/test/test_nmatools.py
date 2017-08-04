# -*- coding: utf-8 -*-
# TAMkin is a post-processing toolkit for normal mode analysis, thermochemistry
# and reaction kinetics.
# Copyright (C) 2008-2012 Toon Verstraelen <Toon.Verstraelen@UGent.be>, An Ghysels
# <An.Ghysels@UGent.be> and Matthias Vandichel <Matthias.Vandichel@UGent.be>
# Center for Molecular Modeling (CMM), Ghent University, Ghent, Belgium; all
# rights reserved unless otherwise stated.
#
# This file is part of TAMkin.
#
# TAMkin is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# In addition to the regulations of the GNU General Public License,
# publications and communications based in parts on this program or on
# parts of this program are required to cite the following article:
#
# "TAMkin: A Versatile Package for Vibrational Analysis and Chemical Kinetics",
# An Ghysels, Toon Verstraelen, Karen Hemelsoet, Michel Waroquier and Veronique
# Van Speybroeck, Journal of Chemical Information and Modeling, 2010, 50,
# 1736-1750W
# http://dx.doi.org/10.1021/ci100099g
#
# TAMkin is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
#--


import os
import numpy as np
import pkg_resources
import unittest

from molmod import lightspeed, centimeter
from molmod.test.common import tmpdir

from tamkin import *


__all__ = ["NMAToolsTestCase"]


class NMAToolsTestCase(unittest.TestCase):

    def test_load_coordinates_charmm(self):
        molecule = load_molecule_charmm(
            pkg_resources.resource_filename(__name__, "../data/test/an/ethanol.cor"),
            pkg_resources.resource_filename(__name__, "../data/test/an/ethanol.hess.full"))
        coordinates, masses, symbols = load_coordinates_charmm(
            pkg_resources.resource_filename(__name__, "../data/test/an/ethanol.cor"))
        for at in range(9):
            self.assertAlmostEqual(molecule.masses[at], masses[at], 3)
            for j in range(3):
                self.assertAlmostEqual(molecule.coordinates[at,j], coordinates[at,j], 3)

    def test_load_modes_charmm(self):
        molecule = load_molecule_charmm(
            pkg_resources.resource_filename(__name__, "../data/test/an/ethanol.cor"),
            pkg_resources.resource_filename(__name__, "../data/test/an/ethanol.hess.full"))

        # full Hessian
        nma = NMA(molecule)
        modes2, freqs2, masses2 = load_modes_charmm(
            pkg_resources.resource_filename(__name__, "../data/test/an/ethanol.modes.full"))
        for index in range(6,27):
            for j in range(27):
                self.assertAlmostEqual(abs(nma.modes[j,index]), abs(modes2[j,index]), 7)
        for at in range(27):
            self.assertAlmostEqual(nma.freqs[at], freqs2[at], 7)

        # MBH
        nma = NMA(molecule, MBH([[5,6,7,8]]))
        modes2, freqs2, masses2 = load_modes_charmm(
            pkg_resources.resource_filename(__name__, "../data/test/an/ethanol.modes.mbh"))
        for index in range(6,21):
            for j in range(27):
                self.assertAlmostEqual(abs(nma.modes[j,index]), abs(modes2[j,index]), 7)
        for at in range(21):
            self.assertAlmostEqual(nma.freqs[at], freqs2[at], 7)

    def test_overlap(self):
        molecule = load_molecule_charmm(
            pkg_resources.resource_filename(__name__, "../data/test/an/ethanol.cor"),
            pkg_resources.resource_filename(__name__, "../data/test/an/ethanol.hess.full"))
        nma1 = NMA(molecule)
        fixed = load_indices(
            pkg_resources.resource_filename(__name__, "../data/test/an/fixed.06.txt"))
        nma2 = NMA(molecule, PHVA(fixed))
        overlap = compute_overlap(nma1, nma2)
        overlap = compute_overlap((nma1.modes, nma1.freqs), (nma2.modes, nma2.freqs))
        overlap = compute_overlap(nma1.modes, nma2.modes)
        overlap = compute_overlap(nma1.modes[:,0], nma2.modes[:,0])
        # TODO
        #self.assertAlmostEqual()

    def test_delta_vector(self):
        # from charmmcor
        coor1,masses1,symb1 = load_coordinates_charmm(
            pkg_resources.resource_filename(__name__, "../data/test/an/ethanol.cor"))
        coor2,masses2,symb2 = load_coordinates_charmm(
            pkg_resources.resource_filename(__name__, "../data/test/an/ethanol.2.cor"))
        delta = compute_delta(coor1, coor2)
        # TODO
        #self.assertAlmostEqual()
        delta = compute_delta(coor1, coor2, masses=masses1, normalize=True)
        #self.assertAlmostEqual()

    def test_eigenvalue_sensitivity(self):
        molecule = load_molecule_charmm(
            pkg_resources.resource_filename(__name__, "../data/test/an/ethanol.cor"),
            pkg_resources.resource_filename(__name__, "../data/test/an/ethanol.hess.full"))
        nma = NMA(molecule)
        for i in range(7,27):
            sensit = compute_sensitivity_freq(nma, i)
            self.assertAlmostEqual(np.sum((np.dot(sensit,nma.modes)-nma.modes)**2,0)[i],0.0,9)

    def test_create_blocks_peptide_charmm(self):
        blocks1 = create_blocks_peptide_charmm(
            pkg_resources.resource_filename(__name__, "../data/test/charmm/crambin.crd"),
            "RTB", blocksize=1)
        blocks2 = create_blocks_peptide_charmm(
            pkg_resources.resource_filename(__name__, "../data/test/charmm/crambin.crd"),
            "RTB", blocksize=2)
        self.assertEqual(len(blocks1)/2+1, len(blocks2))
        subs1 = create_subs_peptide_charmm(
            pkg_resources.resource_filename(__name__, "../data/test/charmm/crambin.crd"),
            frequency=1)
        subs2 = create_subs_peptide_charmm(
            pkg_resources.resource_filename(__name__, "../data/test/charmm/crambin.crd"),
            frequency=2)
        self.assertEqual(len(subs1)/2, len(subs2))

        blocks = create_blocks_peptide_charmm(
            pkg_resources.resource_filename(__name__, "../data/test/charmm/crambin.crd"))
        self.assertEqual(len(blocks), 91)
        blocks = create_blocks_peptide_charmm(
            pkg_resources.resource_filename(__name__, "../data/test/charmm/crambin.crd"),
            "dihedral")
        self.assertEqual(len(blocks), 91)
        blocks = create_blocks_peptide_charmm(
            pkg_resources.resource_filename(__name__, "../data/test/charmm/crambin.crd"),
            "RHbending")
        self.assertEqual(len(blocks),136)
        blocks = create_blocks_peptide_charmm(
            pkg_resources.resource_filename(__name__, "../data/test/charmm/crambin.crd"),
            "normal")
        self.assertEqual(len(blocks), 91)

    def test_plot_spectrum(self):
        molecule = load_molecule_charmm(
            pkg_resources.resource_filename(__name__, "../data/test/an/ethanol.cor"),
            pkg_resources.resource_filename(__name__, "../data/test/an/ethanol.hess.full"))
        nma = NMA(molecule)
        invcm = lightspeed/centimeter
        with tmpdir(__name__, 'test_plot_spectrum') as dn:
            plot_spectrum_lines(
                os.path.join(dn, "spectrum-lines.1.png"),
                [nma.freqs, nma.freqs], title="standard settings")
            plot_spectrum_lines(
                os.path.join(dn, "spectrum-lines.2.png"),
                [nma.freqs, nma.freqs], low=-10.0*invcm, high=500.0*invcm, title="zoom")

            plot_spectrum_dos(
                os.path.join(dn, "spectrum-dos.1.png"),
                [nma.freqs], title="standard settings")
            plot_spectrum_dos(
                os.path.join(dn, "spectrum-dos.2.png"),
                [nma.freqs], low=-10.0*invcm, high=1500.0*invcm, title="zoom")
            plot_spectrum_dos(
                os.path.join(dn, "spectrum-dos.3.png"),
                [nma.freqs], low=-10.0*invcm, high=1500.0*invcm, width=50.0*invcm, title="width")
            plot_spectrum_dos(
                os.path.join(dn, "spectrum-dos.4.png"),
                [nma.freqs], low=-10.0*invcm, high=1500.0*invcm, width=50.0*invcm, step=20.0*invcm, title="step size")
            plot_spectrum_dos(
                os.path.join(dn, "spectrum-dos.5.png"),
                [nma.freqs, nma.freqs*1.1], title="two spectra")
            plot_spectrum_dos(
                os.path.join(dn, "spectrum-dos.6.png"),
                [nma.freqs, nma.freqs*1.1], all_amps=[1.0,2.0], title="different amplitude")

    def test_create_enm_molecule(self):
        molecule = load_molecule_charmm(
            pkg_resources.resource_filename(__name__, "../data/test/an/ethanol.cor"),
            pkg_resources.resource_filename(__name__, "../data/test/an/ethanol.hess.full"))
        selected = range(5)
        mol = create_enm_molecule(molecule, selected)
        nma = NMA(mol)

        mol = create_enm_molecule(molecule.coordinates, rcut=5)
        nma = NMA(mol)

        mol = create_enm_molecule(molecule.coordinates, selected, masses=np.ones(molecule.size)*2.0, rcut=5)
        nma = NMA(mol)
