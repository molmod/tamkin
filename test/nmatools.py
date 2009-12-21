# TAMkin is a post-processing toolkit for thermochemistry and kinetics analysis.
# Copyright (C) 2008-2009 Toon Verstraelen <Toon.Verstraelen@UGent.be>,
# Matthias Vandichel <Matthias.Vandichel@UGent.be> and
# An Ghysels <An.Ghysels@UGent.be>, Center for Molecular Modeling (CMM), Ghent
# University, Ghent, Belgium; all rights reserved unless otherwise stated.
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
# parts of this program are required to cite the following five articles:
#
# "Vibrational Modes in partially optimized molecular systems.", An Ghysels,
# Dimitri Van Neck, Veronique Van Speybroeck, Toon Verstraelen and Michel
# Waroquier, Journal of Chemical Physics, Vol. 126 (22): Art. No. 224102, 2007
# DOI:10.1063/1.2737444
#
# "Cartesian formulation of the Mobile Block Hesian Approach to vibrational
# analysis in partially optimized systems", An Ghysels, Dimitri Van Neck and
# Michel Waroquier, Journal of Chemical Physics, Vol. 127 (16), Art. No. 164108,
# 2007
# DOI:10.1063/1.2789429
#
# "Calculating reaction rates with partial Hessians: validation of the MBH
# approach", An Ghysels, Veronique Van Speybroeck, Toon Verstraelen, Dimitri Van
# Neck and Michel Waroquier, Journal of Chemical Theory and Computation, Vol. 4
# (4), 614-625, 2008
# DOI:10.1021/ct7002836
#
# "Mobile Block Hessian approach with linked blocks: an efficient approach for
# the calculation of frequencies in macromolecules", An Ghysels, Veronique Van
# Speybroeck, Ewald Pauwels, Dimitri Van Neck, Bernard R. Brooks and Michel
# Waroquier, Journal of Chemical Theory and Computation, Vol. 5 (5), 1203-1215,
# 2009
# DOI:10.1021/ct800489r
#
# "Normal modes for large molecules with arbitrary link constraints in the
# mobile block Hessian approach", An Ghysels, Dimitri Van Neck, Bernard R.
# Brooks, Veronique Van Speybroeck and Michel Waroquier, Journal of Chemical
# Physics, Vol. 130 (18), Art. No. 084107, 2009
# DOI:10.1063/1.3071261
#
# TAMkin is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
# --


from tamkin import *

import unittest, numpy
import pylab

__all__ = ["NMAToolsTestCase"]


class NMAToolsTestCase(unittest.TestCase):

    def test_load_coordinates_charmm(self):
        molecule = load_molecule_charmm("input/an/ethanol.cor", "input/an/ethanol.hess.full")
        symbols, coordinates, masses = load_coordinates_charmm("input/an/ethanol.cor")
        for at in range(9):
            self.assertAlmostEqual(molecule.masses[at], masses[at], 3)
            for j in range(3):
                self.assertAlmostEqual(molecule.coordinates[at,j], coordinates[at,j], 3)

    def test_load_modes_charmm(self):
        molecule = load_molecule_charmm("input/an/ethanol.cor", "input/an/ethanol.hess.full")

        # full Hessian
        nma = NMA(molecule)
        freqs2, modes2 = load_modes_charmm("input/an/ethanol.modes.full")
        for index in range(6,27):
            for j in range(27):
                self.assertAlmostEqual(abs(nma.modes[j,index]), abs(modes2[j,index]), 7)
        for at in range(27):
            self.assertAlmostEqual(nma.freqs[at], freqs2[at], 7)

        # MBH
        nma = NMA(molecule, MBH([[5,6,7,8]]))
        freqs2, modes2 = load_modes_charmm("input/an/ethanol.modes.mbh")
        for index in range(6,21):
            for j in range(27):
                self.assertAlmostEqual(abs(nma.modes[j,index]), abs(modes2[j,index]), 7)
        for at in range(21):
            self.assertAlmostEqual(nma.freqs[at], freqs2[at], 7)


    def test_overlap(self):
        molecule = load_molecule_charmm("input/an/ethanol.cor","input/an/ethanol.hess.full")
        nma1 = NMA(molecule)
        fixed = load_fixed_txt("input/an/fixed.06.txt")
        nma2 = NMA(molecule, PHVA(fixed))
        overlap = calculate_overlap_nma(nma1, nma2)
        # TODO
        #self.assertAlmostEqual()


    def test_Delta_vector(self):
        # from charmmcor
        symb1,coor1,masses1 = load_coordinates_charmm("input/an/ethanol.cor")
        symb2,coor2,masses2 = load_coordinates_charmm("input/an/ethanol.2.cor")
        Delta = get_Delta_vector(coor1, coor2)
        # TODO
        #self.assertAlmostEqual()
        Delta = get_Delta_vector(coor1, coor2, masses = masses1, normalize = True)
        #self.assertAlmostEqual()

    def test_eigenvalue_sensitivity(self):
        molecule = load_molecule_charmm("input/an/ethanol.cor","input/an/ethanol.hess.full")
        nma = NMA(molecule)
        for i in range(7,27):
            sensit = calculate_sensitivity_freq(nma, i)
            self.assertAlmostEqual(numpy.sum((numpy.dot(sensit,nma.modes)-nma.modes)**2,0)[i],0.0,9)




