#!/usr/bin/env python
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


# Import the tamkin library.
from tamkin import *

# Load the gaussian data:
mol_react = load_molecule_g03fchk("react.fchk")
mol_ts = load_molecule_g03fchk("ts.fchk")
# Perform the normal mode analysis
nma_react = NMA(mol_react, ConstrainExt())
nma_ts = NMA(mol_ts, ConstrainExt(2e-4))
# Construct the two partition functions.
pf_react = PartFun(nma_react, [])
pf_ts = PartFun(nma_ts, [])

# Define a kinetic model for the chemical reaction. These are the mandatory arguments:
#  1) a list of reactant partition functions
#     (one for unimolecular, two for bimolecular, ...)
#  2) the transition state partition function
# There is one more optional argument:
#  3) tunneling: a model for the tunelling correction
km = KineticModel([pf_react], pf_ts, tunneling=None)

# Analyze the chemical reaction. These are the arguments:
#  1) A kinetic model
#  2) the starting temperature for the fit
#  3) the final temperature for the fit
# The following argument is optional:
#  4) temp_step: The interval on the temperature grid in Kelvin, 10 is default
ra = ReactionAnalysis(km, 670, 770)
# Make the Arrhenius plot
ra.plot_arrhenius("arrhenius.png")
# Write all results to a file.
ra.write_to_file("reaction.txt")
