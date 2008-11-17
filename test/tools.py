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


from tamkin.tools import fit_kin

import unittest
import numpy

__all__ = ["FitKinCase"]

class FitKinCase(unittest.TestCase):
    def test_fit_kin(self):
        k = numpy.array([7.9473102E+05, 9.8300444E+05, 1.2085262E+06, 1.4771808E+06, 1.7955340E+06, 2.1708793E+06, 2.6112829E+06, 3.1256298E+06, 3.7236678E+06, 4.4160510E+06, 5.2143822E+06])
        temps = numpy.array([670,680,690,700,710,720,730,740,750,760,770])
        fit_kin(temps,k)



