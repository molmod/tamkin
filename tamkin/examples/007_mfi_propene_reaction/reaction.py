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

# load the fixed atoms
#fixed_atoms = load_fixed_g03com("Zp_p_react.14mei.com")
fixed_react = load_fixed_g03com("Zp_p_react.28aug.com")
fixed_trans = load_fixed_g03com("Zp_p_TS.28aug.com")
# load the gaussian data
mol_react = load_molecule_g03fchk("Zp_p_react.28aug.fchk", "Zp_p_react.14mei.fchk")
mol_trans = load_molecule_g03fchk("Zp_p_TS.28aug.fchk", "5Tp_p_TS.oniom21apr_HF.fchk")
# Construct the three partition functions.
pf_react = PartFun(NMA(mol_react, PHVA(fixed_react)))
pf_trans = PartFun(NMA(mol_trans, PHVA(fixed_trans)))


# Define a kinetic model for the chemical reaction. These are the mandatory arguments:
#  1) a list of reactant partition functions
#     (one for unimolecular, two for bimolecular, ...)
#  2) the transition state partition function
# There is one more optional argument:
#  3) tunneling: a model for the tunelling correction
km = KineticModel([pf_react], pf_trans, tunneling=None)

# Analyze the chemical reaction. These are the arguments:
#  1) A kinetic model
#  2) the starting temperature for the fit
#  3) the final temperature for the fit
# The following argument is optional:
#  4) temp_step: The interval on the temperature grid in Kelvin, 10 is default
ra = ReactionAnalysis(km, 670, 770)
ra.plot_arrhenius("arrhenius.png") # make the Arrhenius plot

# Estimate the error on the kinetic parameters due to level of theory artifacts
# with Monte Carlo sampling. The monte_carlo method takes three optional
# arguments:
#  1) freq_error: the absolute stochastic error on the frequencies (default=1*invcm)
#  2) energy_error: the absolute error on the energy (default=0.0)
#  3) num_iter: the number of monte carlo samples (default=100)
ra.monte_carlo()
# plot the parameters, this includes the monte carlo results
ra.plot_parameters("parameters.png")
# write all results to a file.
ra.write_to_file("reaction.txt") # summary of the analysis
