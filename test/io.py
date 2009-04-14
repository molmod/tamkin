# TAMkin is a post-processing toolkit for thermochemistry and kinetics analysis.
# Copyright (C) 2008-2009 Toon Verstraelen <Toon.Verstraelen@UGent.be>,
# Matthias Vandichel <Matthias.Vandichel@UGent.be> and
# An Ghysels <An.Ghysels@UGent.be>
#
# This file is part of TAMkin.
#
# TAMkin is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
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

from molmod.data.periodic import periodic
from molmod.units import angstrom, cm, amu
from molmod.constants import lightspeed

import unittest, numpy


__all__ = ["IOTestCase"]


class IOTestCase(unittest.TestCase):
    def test_load_fixed_g03com(self):
        fixed_atoms = load_fixed_g03com("input/mat/Zp_p_prod.18aug.com")
        self.assertEqual(len(fixed_atoms), 48)
        self.assertEqual(fixed_atoms, range(114,114+48))

    def test_load_molecule_g03fchk(self):
        atoms = 181
        molecule = load_molecule_g03fchk("input/mat/Zp_p_react.28aug.fchk","input/mat/Zp_p_react.14mei.fchk")

        self.assertEqual(molecule.hessian.shape,(atoms*3,atoms*3))
        self.assertAlmostEqual(molecule.energy, -18613.135744186180, 7)

    def test_load_molecule_cp2k(self):
        molecule = load_molecule_cp2k("input/cp2k/pentane/opt.xyz", "input/cp2k/pentane/sp.out", "input/cp2k/pentane/freq.out")
        self.assertAlmostEqual(molecule.energy, 0.012255059530862)
        self.assertEqual(molecule.multiplicity, 1)
        self.assertEqual(molecule.numbers[0], 6)
        self.assertEqual(molecule.numbers[4], 1)
        self.assertAlmostEqual(molecule.masses[0], periodic[6].mass)
        self.assertAlmostEqual(molecule.masses[4], periodic[1].mass)
        self.assertAlmostEqual(molecule.coordinates[5,1]/angstrom, 13.9457396458)
        self.assertAlmostEqual(molecule.gradient[0,2], 0.0000000038, 9)
        self.assertAlmostEqual(molecule.gradient[11,0], 0.0000000177, 9)
        self.assertAlmostEqual(molecule.hessian[0,0], 1.08660340, 6)
        self.assertAlmostEqual(molecule.hessian[-1,-1], 0.528947590, 6)

    def test_load_molecule_cpmd(self):
        molecule = load_molecule_cpmd("input/cpmd/damp.out", "input/cpmd/GEOMETRY.xyz", "input/cpmd/MOLVIB")
        self.assertAlmostEqual(molecule.energy, -17.14142079)
        self.assertEqual(molecule.multiplicity, 1)
        self.assertEqual(molecule.numbers[0], 8)
        self.assertEqual(molecule.numbers[2], 1)
        self.assertAlmostEqual(molecule.masses[0]/amu, 15.999400)
        self.assertAlmostEqual(molecule.masses[2]/amu, 1.007970)
        self.assertAlmostEqual(molecule.coordinates[2,0]/angstrom, .907207019799)
        self.assertAlmostEqual(molecule.gradient[0,2], 0.0)
        self.assertAlmostEqual(molecule.gradient[2,0], 0.0)
        self.assertAlmostEqual(molecule.hessian[0,0], 0.530679165332463953, 6)
        self.assertAlmostEqual(molecule.hessian[-1,-1], 0.921045226686428159E-01, 6)

    def test_checkpoint(self):
        molecule = load_molecule_cp2k("input/cp2k/pentane/opt.xyz", "input/cp2k/pentane/sp.out", "input/cp2k/pentane/freq.out")
        nma1 = NMA(molecule)
        nma1.write_to_file("output/test.chk")
        nma2 = NMA.read_from_file("output/test.chk")

        self.assertEqual(nma1.freqs.shape, nma2.freqs.shape)
        self.assertEqual(nma1.modes.shape, nma2.modes.shape)
        self.assertEqual(nma1.masses.shape, nma2.masses.shape)
        self.assertEqual(nma1.numbers.shape, nma2.numbers.shape)
        self.assertEqual(nma1.coordinates.shape, nma2.coordinates.shape)
        self.assertEqual(nma1.inertia_tensor.shape, nma2.inertia_tensor.shape)

        self.assert_(abs(nma1.freqs - nma2.freqs).max()/abs(nma1.freqs).max() < 1e-15)
        self.assert_(abs(nma1.modes - nma2.modes).max()/abs(nma1.modes).max() < 1e-15)
        self.assert_(abs(nma1.masses - nma2.masses).max()/abs(nma1.masses).max() < 1e-15)
        self.assert_(abs(nma1.coordinates - nma2.coordinates).max()/abs(nma1.coordinates).max() < 1e-15)
        self.assert_(abs(nma1.inertia_tensor - nma2.inertia_tensor).max()/abs(nma1.inertia_tensor).max() < 1e-15)
        self.assert_((nma1.numbers==nma2.numbers).all())

        self.assertAlmostEqual(nma1.mass, nma2.mass)
        self.assertAlmostEqual(nma1.energy, nma2.energy)
        self.assertEqual(nma1.multiplicity, nma2.multiplicity)
        self.assertEqual(nma1.symmetry_number, nma2.symmetry_number)



