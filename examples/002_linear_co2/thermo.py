#!/usr/bin/env python
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


# Import the tamkin libarary.
from tamkin import *
# Load the gaussian data.
molecule = load_molecule_g03fchk("gaussian.fchk")

# Perform the normal mode analysis
nma = NMA(molecule)
# Construct a partition function object with the typical gas phase contributions.
pf = PartFun(nma, [ExternalTranslation(), ExternalRotation(2)])

# Write some general information about the molecule and the partition function
# to a file.
pf.write_to_file("partfun.txt")

# Write an extensive overview of the thermodynamic properties to a file:
ta = ThermoAnalysis(pf, [300,400,500,600])
ta.write_to_file("thermo.csv")




