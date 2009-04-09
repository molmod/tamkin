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
# Import unit conversin factors
from molmod.units import *
# Import numpy stuff
from numpy import *


# Load the gaussian data.
molecule = load_molecule_g03fchk("gaussian.fchk")
print "Energy [kJ/mol]:", molecule.energy/kjmol


# A) Perform a standard normal mode analysis
nma1 = NMA(molecule)
print "The zero eigenmodes: %s" % nma1.zeros
# Construct a partition function with the typical gas phase contributions.
pf1 = PartFun(nma1, [ExternalTranslation(), ExternalRotation(2)])
# 2 is the rotational symmetry number.
print "Heat capacity at 300K [J/(mol*K)]:", pf1.heat_capacity(300*K)/(J/mol/K)
# Write some general information about the molecule and the partition function
# to a file.
pf1.write_to_file("partfun1.txt")
# Write an extensive overview of the thermodynamic properties to a file:
ta1 = ThermoAnalysis(pf1, arange(300,601,5))
ta1.write_to_file("thermo1.csv")


# B) Perform a normal mode analysis with constrained external degrees of freedom.
# This implies that the vibrational analysis is performed in 3N-6 dof.
nma2 = NMA(molecule, ConstrainExt())
print "The zero eigenmodes: %s" % nma2.zeros
# Construct a partition function with the typical gas phase contributions.
pf2 = PartFun(nma2, [ExternalTranslation(), ExternalRotation(2)])
# 2 is the rotational symmetry number.
print "Heat capacity at 300K [J/(mol*K)]:", pf2.heat_capacity(300*K)/(J/mol/K)
# Write some general information about the molecule and the partition function
# to a file.
pf2.write_to_file("partfun2.txt")
# Write an extensive overview of the thermodynamic properties to a file:
ta2 = ThermoAnalysis(pf2, arange(300,601,5))
ta2.write_to_file("thermo2.csv")




