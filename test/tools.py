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


from tamkin import *

from molmod.units import kjmol, atm
from molmod.constants import boltzmann

import unittest
import numpy


__all__ = ["ToolsTestCase"]


class ToolsTestCase(unittest.TestCase):
    def test_reaction_analysis_sterck(self):
        pf_react1 = PartFun(NMA(load_molecule_g03fchk("input/sterck/aa.fchk")), [ExternalTranslation(), ExternalRotation(1)])
        pf_react2 = PartFun(NMA(load_molecule_g03fchk("input/sterck/aarad.fchk")), [ExternalTranslation(), ExternalRotation(1)])
        pf_trans = PartFun(NMA(load_molecule_g03fchk("input/sterck/paats.fchk")), [ExternalTranslation(), ExternalRotation(1)])

        ra = ReactionAnalysis([pf_react1, pf_react2], pf_trans, 280, 360)
        # not a very accurate check because the fit is carried out differently
        # in the fancy excel file where these numbers come from.
        self.assertAlmostEqual(ra.Ea/kjmol, 25.96, 1)
        self.assertAlmostEqual(numpy.log(ra.A/ra.unit), numpy.log(2.29E+02), 1)

        ra.plot("output/arrhenius_aa.png")
        ra.monte_carlo()
        ra.write_to_file("output/reaction_aa.txt")
        ra.plot_parameters("output/parameters_aa.png")

    def test_reaction_analysis_mat(self):
        pf_react = PartFun(NMA(load_molecule_g03fchk("input/mat/5Te_etheen_react_deel2.fchk")), [ExternalTranslation(), ExternalRotation(1)])
        pf_trans = PartFun(NMA(load_molecule_g03fchk("input/mat/5Te_etheen_ts_deel2_punt108_freq.fchk")), [ExternalTranslation(), ExternalRotation(1)])

        ra = ReactionAnalysis([pf_react], pf_trans, 100, 1200)
        # not a very accurate check because the fit is carried out differently
        # in the fancy excel file where these numbers come from.
        self.assertAlmostEqual(ra.Ea/kjmol, 160.6, 0)
        self.assertAlmostEqual(numpy.log(ra.A/ra.unit), numpy.log(3.33e10), 0)
        ra.plot("output/arrhenius_mat1.png")
        ra.monte_carlo()
        ra.write_to_file("output/reaction_mat1.txt")
        ra.plot_parameters("output/parameters_mat1.png")

        ra = ReactionAnalysis([pf_react], pf_trans, 670, 770)
        # not a very accurate check because the fit is carried out differently
        # in the fancy excel file where these numbers come from.
        self.assertAlmostEqual(ra.Ea/kjmol, 161.9, 1)
        self.assertAlmostEqual(numpy.log(ra.A/ra.unit), numpy.log(4.08e10), 0)
        ra.plot("output/arrhenius_mat2.png")
        ra.monte_carlo()
        ra.write_to_file("output/reaction_mat2.txt")
        ra.plot_parameters("output/parameters_mat2.png")

    def test_thermo_analysis_mat(self):
        # just a blind test to see test whether the code does not crash.
        pf = PartFun(NMA(load_molecule_g03fchk("input/mat/5Te_etheen_react_deel2.fchk")), [ExternalTranslation(), ExternalRotation(1)])
        ta = ThermoAnalysis(pf, [200,300,400,500,600,700,800,900])
        ta.write_to_file("output/thermo_mat2.csv")

