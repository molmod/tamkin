# TAMkin is a post-processing toolkit for thermochemistry and kinetics analysis.
# Copyright (C) 2008-2010 Toon Verstraelen <Toon.Verstraelen@UGent.be>,
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

from molmod.periodic import periodic
from molmod.units import angstrom, amu, calorie, avogadro, electronvolt
from molmod.constants import lightspeed

import unittest, numpy


__all__ = ["DataTestCase"]


class DataTestCase(unittest.TestCase):

    def test_get_submolecule(self):
        molecule = load_molecule_charmm("input/an/ethanol.cor","input/an/ethanol.hess.full")
        select = [3,4,2]
        molecule2 = molecule.get_submolecule(select)
        for i,at in enumerate(select):
            self.assertAlmostEqual(molecule.numbers[at], molecule2.numbers[i])
            self.assertAlmostEqual(molecule.masses[at], molecule2.masses[i])
            self.assertEqual(molecule.symbols[at], molecule2.symbols[i])
            for mu in range(3):
                self.assertAlmostEqual(molecule.coordinates[at,mu],molecule2.coordinates[i,mu])
                self.assertAlmostEqual(molecule.gradient[at,mu],molecule2.gradient[i,mu])
                for j,at2 in enumerate(select):
                    for nu in range(3):
                        self.assertAlmostEqual(molecule.hessian[3*at+mu,3*at2+nu],molecule2.hessian[3*i+mu,3*j+nu])
        self.assertEqual(molecule.multiplicity, molecule2.multiplicity)
        self.assertEqual(molecule.symmetry_number, molecule2.symmetry_number)
        self.assertEqual(molecule.periodic, molecule2.periodic)
        self.assertEqual(molecule.energy, molecule2.energy)

    def test_get_submolecule_cp2k(self):
        molecule = load_molecule_cp2k("input/cp2k/pentane/opt.xyz", "input/cp2k/pentane/sp.out", "input/cp2k/pentane/freq.out")
        select = range(5)+[9,11,14]
        molecule2 = molecule.get_submolecule(select, title="this is submol", energy=5., periodic=False, symmetry_number=6)  # just trying out something
        for i,at in enumerate(select):
            self.assertAlmostEqual(molecule.numbers[at], molecule2.numbers[i])
            self.assertAlmostEqual(molecule.masses[at], molecule2.masses[i])
            for mu in range(3):
                self.assertAlmostEqual(molecule.coordinates[at,mu],molecule2.coordinates[i,mu])
                self.assertAlmostEqual(molecule.gradient[at,mu],molecule2.gradient[i,mu])
                for j,at2 in enumerate(select):
                    for nu in range(3):
                        self.assertAlmostEqual(molecule.hessian[3*at+mu,3*at2+nu],molecule2.hessian[3*i+mu,3*j+nu])
        self.assertEqual(molecule.multiplicity, molecule2.multiplicity)
        self.assertEqual(molecule2.symmetry_number, 6)
        self.assertEqual(molecule2.periodic, False)
        self.assertEqual(molecule2.energy, 5.)

    def test_translate_pbc(self):
        molecule = load_molecule_cp2k("input/cp2k/pentane/opt.xyz", "input/cp2k/pentane/sp.out", "input/cp2k/pentane/freq.out")
        self.assertAlmostEqual(molecule.unit_cell.matrix[1,1]/angstrom, 30.000,3)
        self.assertAlmostEqual(molecule.unit_cell.matrix[0,1]/angstrom, 0.000,3)
        self.assertAlmostEqual(molecule.coordinates[5,1]/angstrom, 13.9457396458)
        selected = range(6)+[11,14]
        molecule2 = translate_pbc(molecule, selected, [1,-1,0])
        self.assertAlmostEqual(molecule2.coordinates[5,1]/angstrom, 13.9457396458-30.)
        for i in range(molecule.size):
            self.assertAlmostEqual(molecule.numbers[i], molecule2.numbers[i])
            self.assertAlmostEqual(molecule.masses[i], molecule2.masses[i])
            for mu in range(3):
                if i not in selected:
                    self.assertAlmostEqual(molecule.coordinates[i,mu],molecule2.coordinates[i,mu])
                else:
                    if mu is 0:
                        self.assertAlmostEqual(molecule.coordinates[i,mu]+30.*angstrom,molecule2.coordinates[i,mu])
                    if mu is 1:
                        self.assertAlmostEqual(molecule.coordinates[i,mu]-30.*angstrom,molecule2.coordinates[i,mu])
                    if mu is 2:
                        self.assertAlmostEqual(molecule.coordinates[i,mu],molecule2.coordinates[i,mu])
                self.assertAlmostEqual(molecule.gradient[i,mu],molecule2.gradient[i,mu])
                for j in selected:
                    for nu in range(3):
                        self.assertAlmostEqual(molecule.hessian[3*i+mu,3*j+nu],molecule2.hessian[3*i+mu,3*j+nu])
        self.assertEqual(molecule.multiplicity, molecule2.multiplicity)
        self.assertEqual(molecule.symmetry_number, molecule2.symmetry_number)
        self.assertAlmostEqual(molecule2.unit_cell.matrix[0,0]/angstrom, 30.000,3)
        self.assertAlmostEqual(molecule2.unit_cell.matrix[1,2]/angstrom, 0.000,3)

