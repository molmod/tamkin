#!/usr/bin/env python
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


# Import the tamkin libarary.
from tamkin import *
# Import unit conversin factors
from molmod import *
# Import numpy stuff
from numpy import *


# Load the gaussian data.
molecule = load_molecule_g03fchk("gaussian.fchk")
print "Energy [kJ/mol]:", molecule.energy/kjmol


# A) Perform a standard normal mode analysis
nma1 = NMA(molecule)
print "The zero eigenmodes: %s" % nma1.zeros
# Construct a partition function with the typical gas phase contributions.
pf1 = PartFun(nma1, [ExtTrans(), ExtRot()])
# 2 is the rotational symmetry number.
print "Heat capacity at 300K [J/(mol*K)]:", pf1.heat_capacity(300*kelvin)/(joule/mol/kelvin)
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
pf2 = PartFun(nma2, [ExtTrans(), ExtRot(2)])
# 2 is the rotational symmetry number.
print "Heat capacity at 300K [J/(mol*K)]:", pf2.heat_capacity(300*kelvin)/(joule/mol/kelvin)
# Write some general information about the molecule and the partition function
# to a file.
pf2.write_to_file("partfun2.txt")
# Write an extensive overview of the thermodynamic properties to a file:
ta2 = ThermoAnalysis(pf2, arange(300,601,5))
ta2.write_to_file("thermo2.csv")


