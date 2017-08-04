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
import pkg_resources
import numpy as np
import unittest

from molmod.units import kjmol, atm, meter, mol, second
from molmod.constants import boltzmann
from molmod.test.common import tmpdir

from tamkin import *


__all__ = ["PFToolsTestCase"]


class PFToolsTestCase(unittest.TestCase):
    def test_reaction_analysis_sterck(self):
        pf_react1 = PartFun(NMA(load_molecule_g03fchk(
            pkg_resources.resource_filename(__name__, "../data/test/sterck/aa.fchk"))),
            [ExtTrans(cp=False), ExtRot(1)])
        pf_react2 = PartFun(NMA(load_molecule_g03fchk(
            pkg_resources.resource_filename(__name__, "../data/test/sterck/aarad.fchk"))),
            [ExtTrans(cp=False), ExtRot(1)])
        pf_ts = PartFun(NMA(load_molecule_g03fchk(
            pkg_resources.resource_filename(__name__, "../data/test/sterck/paats.fchk"))),
            [ExtTrans(cp=False), ExtRot(1)])

        km = KineticModel([pf_react1, pf_react2], pf_ts)
        ra = ReactionAnalysis(km, 280, 360)
        # not a very accurate check because the fit is carried out differently
        # in the fancy excel file where these numbers come from.
        self.assertAlmostEqual(ra.Ea/kjmol, 25.96, 1)
        self.assertAlmostEqual(km.unit, meter**3/mol/second)
        self.assertEqual(km.unit_name, "m^3 mol^-1 s^-1")
        self.assertAlmostEqual(np.log(ra.A/km.unit), np.log(2.29E+02), 1)

        with tmpdir(__name__, 'test_reaction_analysis_sterck') as dn:
            ra.plot_arrhenius(os.path.join(dn, "arrhenius_aa.png"))
            ra.monte_carlo()
            ra.write_to_file(os.path.join(dn, "reaction_aa.txt"))
            ra.plot_parameters(os.path.join(dn, "parameters_aa.png"))

    def test_reaction_analysis_mat(self):
        pf_react = PartFun(NMA(load_molecule_g03fchk(
            pkg_resources.resource_filename(__name__, "../data/test/mat5T/react.fchk"))), [])
        pf_ts = PartFun(NMA(load_molecule_g03fchk(
            pkg_resources.resource_filename(__name__, "../data/test/mat5T/ts.fchk"))), [])

        km = KineticModel([pf_react], pf_ts)
        ra = ReactionAnalysis(km, 100, 1200, temp_step=50)
        # not a very accurate check because the fit is carried out differently
        # in the fancy excel file where these numbers come from.
        self.assertAlmostEqual(ra.Ea/kjmol, 160.6, 0)
        self.assertAlmostEqual(km.unit, 1.0/second)
        self.assertEqual(km.unit_name, "s^-1")
        self.assertAlmostEqual(np.log(ra.A/km.unit), np.log(3.33e10), 0)
        with tmpdir(__name__, 'test_reaction_analysis_mat1') as dn:
            ra.plot_arrhenius(os.path.join(dn, "arrhenius_mat1.png"))
            ra.monte_carlo()
            ra.write_to_file(os.path.join(dn, "reaction_mat1.txt"))
            ra.plot_parameters(os.path.join(dn, "parameters_mat1.png"))

        wigner = Wigner(pf_ts) # Blind test of the wigner correction and
        # the corrected reaction analysis.
        km = KineticModel([pf_react], pf_ts, tunneling=wigner)
        ra = ReactionAnalysis(km, 100, 1200, temp_step=50)
        with tmpdir(__name__, 'test_reaction_analysis_mat2') as dn:
            ra.plot_arrhenius(os.path.join(dn, "arrhenius_mat1w.png"))
            ra.monte_carlo()
            ra.write_to_file(os.path.join(dn, "reaction_mat1w.txt"))
            ra.plot_parameters(os.path.join(dn, "parameters_mat1w.png"))

        km = KineticModel([pf_react], pf_ts)
        ra = ReactionAnalysis(km, 670, 770)
        # not a very accurate check because the fit is carried out differently
        # in the fancy excel file where these numbers come from.
        self.assertAlmostEqual(ra.Ea/kjmol, 161.9, 1)
        self.assertAlmostEqual(np.log(ra.A/km.unit), np.log(4.08e10), 0)
        with tmpdir(__name__, 'test_reaction_analysis_mat3') as dn:
            ra.plot_arrhenius(os.path.join(dn, "arrhenius_mat2.png"))
            ra.monte_carlo()
            ra.write_to_file(os.path.join(dn, "reaction_mat2.txt"))
            ra.plot_parameters(os.path.join(dn, "parameters_mat2.png"))

    def test_thermo_analysis_mat(self):
        # just a blind test to see test whether the code does not crash.
        pf = PartFun(NMA(load_molecule_g03fchk(
            pkg_resources.resource_filename(__name__, "../data/test/mat5T/react.fchk"))),
            [ExtTrans(), ExtRot(1)])
        ta = ThermoAnalysis(pf, [200,300,400,500,600,700,800,900])
        with tmpdir(__name__, 'test_thermo_analysis_mat') as dn:
            ta.write_to_file(os.path.join(dn, "thermo_mat2.csv"))
