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


from tamkin.tools import ReactionAnalysis
from tamkin.partf import PartFun, ExternalTranslation, ExternalRotation, Electronic
from tamkin.io import load_molecule_g03fchk

from molmod.units import kjmol, atm
from molmod.constants import boltzmann

import unittest
import numpy

__all__ = ["ToolsTestCase"]

class ToolsTestCase(unittest.TestCase):
    def test_reaction_analysis(self):
        pf_react1 = PartFun(load_molecule_g03fchk("input/sterck/aa.fchk"), [ExternalTranslation(), ExternalRotation(1), Electronic()])
        pf_react2 = PartFun(load_molecule_g03fchk("input/sterck/aarad.fchk"), [ExternalTranslation(), ExternalRotation(1), Electronic()])
        pf_trans = PartFun(load_molecule_g03fchk("input/sterck/paats.fchk"), [ExternalTranslation(), ExternalRotation(1), Electronic()])

        ra = ReactionAnalysis([pf_react1, pf_react2], pf_trans, 280, 360)
        # not a very accurate check because the fit is carried out differently
        # in the fancy excel file where these numbers come from.
        self.assertAlmostEqual(ra.Ea/kjmol, 25.96, 1)
        self.assertAlmostEqual(ra.A/ra.unit, 2.29E+02, -1)

        ra.plot("output/arrhenius_aa.png")
        ra.write_to_file("output/reaction_aa.txt")



