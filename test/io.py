# TAMkin is a post-processing toolkit for thermochemistry and kinetics analysis.
# Copyright (C) 2008 Toon Verstraelen <Toon.Verstraelen@UGent.be>,
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


from tamkin.partf import PartFun, ExternalTranslation, ExternalRotation, Electronic
from tamkin.io import load_fixed_g03com, load_molecule_g03fchk, load_molecule_cp2k

from molmod.data.periodic import periodic
from molmod.units import angstrom

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

        pf = PartFun(molecule, [])
        expected_freqs = numpy.array([ # taken from cp2k output file
            125.923216, 132.737491, 251.881675, 258.656045, 260.028863,
            489.690362, 522.430525, 836.589205, 924.800235, 1069.391906,
            1152.109253, 1182.947666, 1226.649618, 1279.047068, 1494.291889,
            1505.963333, 1526.247839, 1526.770449, 1527.220869, 1554.895642,
            1574.963368, 1664.772539, 1726.540705, 1730.035399, 1730.370082,
            1740.966064, 1740.971286, 1758.656380, 1767.704604, 1790.726516,
            1825.670669, 1980.257356, 2107.134204, 4515.156994, 4526.218064,
            4540.443901, 4578.619346, 4578.624300, 4581.926079, 4582.125771,
            4646.530637, 4657.432936, 4671.093295, 4751.240132, 4751.291450
        ])
        for i in xrange(len(expected_freqs)):
            self.assertAlmostEqual(expected_freqs[i], pf.vibrational.freqs[i+6], 1)

