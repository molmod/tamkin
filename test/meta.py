# -*- coding: utf-8 -*-
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
# "TAMkin: A Versatile Package for Vibrational Analysis and Chemical Kinetics",
# An Ghysels, Toon Verstraelen, Karen Hemelsoet, Michel Waroquier and Veronique
# Van Speybroeck, Journal of Chemical Information and Modeling, Articles ASAP
# (As Soon As Publishable)
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
# --


from tamkin import *

import unittest, os, glob


__all__ = ["MetaTestCase"]


class MetaTestCase(unittest.TestCase):
    def check_example(self, dirname, fn_py):
        root = "../examples"
        self.assert_(os.path.isdir(root))
        retcode = os.system("cd %s; cd %s; ./%s 1> /dev/null 2> /dev/null" % (root, dirname, fn_py))
        self.assertEqual(retcode, 0)

    def test_example_001(self):
        self.check_example("001_ethane", "./thermo.py")

    def test_example_002(self):
        self.check_example("002_linear_co2", "./thermo.py")

    def test_example_003(self):
        self.check_example("003_pentane", "./thermo.py")

    def test_example_005(self):
        self.check_example("005_acrylamide_reaction", "./reaction.py")

    def test_example_006(self):
        self.check_example("006_5T_ethene_reaction", "./reaction.py")

    def test_example_007(self):
        self.check_example("007_mfi_propene_reaction", "./reaction.py")

    def test_example_008(self):
        self.check_example("008_ethane_rotor", "./thermo.py")

    def test_example_009(self):
        self.check_example("009_ethyl_ethene", "./reaction.py")

    def test_example_012(self):
        self.check_example("012_ethyl_ethene_scaling", "./reaction.py")

    def test_example_013(self):
        self.check_example("013_butane", "./thermo.py")

    def test_example_014(self):
        self.check_example("014_pentane_mbh", "./thermo.py")

    def test_example_015(self):
        self.check_example("015_kie", "./reaction.py")

    def test_example_016(self):
        self.check_example("016_modes", "./modes.py")

    def test_example_017(self):
        self.check_example("017_activationkineticmodel", "./reaction.py")

    def test_example_018(self):
        self.check_example("018_physisorption", "./adsorption.py")

    def test_example_019(self):
        self.check_example("019_ethyl_ethene_simple", "./kinetic.py")

    def test_code_quality(self):
        root = "../lib"
        self.assert_(os.path.isdir(root))
        white = (" ", "\t")
        for fn in glob.glob("%s/*.py") + glob.glob("%s/io/*.py"):
            f = file(fn)
            for line in f:
                if line[-2] in white:
                    self.fail("Trailing whitespace in %s." % fn)
