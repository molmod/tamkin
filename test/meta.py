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

    def test_code_quality(self):
        root = "../lib"
        self.assert_(os.path.isdir(root))
        white = (" ", "\t")
        for fn in glob.glob("%s/*.py") + glob.glob("%s/io/*.py"):
            f = file(fn)
            for line in f:
                if line[-2] in white:
                    self.fail("Trailing whitespace in %s." % fn)
