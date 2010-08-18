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
# There are two more optional arguments
#  3) cp: model at constant pressure, default=True
#  4) tunneling: a model for the tunelling correction
km = KineticModel([pf_react], pf_trans, cp=True, tunneling=None)

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
