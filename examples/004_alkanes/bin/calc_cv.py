#!/usr/bin/python
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


# This mini program computes the heat capacity of a molecular system based on
# the variance of the kinetic energy from an NVE simulation.

# import conversion factors and physical constants
from molmod.units import *
from molmod.constants import *

# hardcoded parameters
N = 8 # number of atoms
DOF = 3*N-6 # degrees of freedom
temp = 299.337869691 # kelvin
var = 4.67336041784e-06 # atomic units

# Allen and Tildesley, page 53, formula 2.82
Cv = 0.5*DOF*boltzmann/(1-var/(0.5*DOF*boltzmann**2*temp**2))
print "Heat capacity [J/(mol*K)]: %.f " % (Cv/(J/mol/K))
