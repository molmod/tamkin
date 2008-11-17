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


import unittest
from tamkin.io import load_fixed_g03com, load_kin_g03fchk

__all__ = ["IOTestCase"]


class IOTestCase(unittest.TestCase):
    def test_load_fixed_g03com(self):
        fixed_atoms = load_fixed_g03com("input/Zp_p_prod.18aug.com")
        self.assertEqual(len(fixed_atoms), 48)
        self.assertEqual(fixed_atoms, range(114,114+48))

    def test_load_kin_g03fchk(self):
        atoms = 181
        data = load_kin_g03fchk("input/Zp_p_react.28aug.fchk","input/Zp_p_react.14mei.fchk")

        self.assertEqual(data["hessian"].shape,(atoms*3,atoms*3))
        self.assertAlmostEqual(data["energy"], -18613.135744186180, 7)



