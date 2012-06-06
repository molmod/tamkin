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


# Import the tamkin libarary.
from tamkin import *
# Import unit conversin factors
from molmod import *
# Import numpy stuff
from numpy import *

###
### reaction can be split up into two parts (HAA = acetic acid):
### activation part: VO(AA)(OH)--H2O + TBHP <=> VO(AA)(OOtBu) + 2H2O
### reaction part (oxygen transfer reaction): VO(AA)(OOtBu) + cyclohexene -> VO(ACAC)(OtBu) + cyclohexene-epoxide
###
# yclohexene  reaction.py  TBHP  TS  VO_AA_OOtBu  water

nma_H2O = NMA(load_molecule_g03fchk("water/gaussian.fchk"))
nma_VO_AA_OH_H2O = NMA(load_molecule_g03fchk("VO_AA_OH_H2O/gaussian.fchk"))
nma_VO_AA_OOtBu = NMA(load_molecule_g03fchk("VO_AA_OOtBu/gaussian.fchk"))
nma_TBHP = NMA(load_molecule_g03fchk("TBHP/gaussian.fchk"))
nma_TS = NMA(load_molecule_g03fchk("TS/TS.fchk"))
nma_cyclohexene = NMA(load_molecule_g03fchk("cyclohexene/gaussian.fchk"))

pf_H2O = PartFun(nma_H2O, [ExtTrans(), ExtRot(2),Vibrations(freq_scaling=1.0, zp_scaling=1.0)])
pf_VO_AA_OH_H2O = PartFun(nma_VO_AA_OH_H2O, [ExtTrans(), ExtRot(1),Vibrations(freq_scaling=1.0, zp_scaling=1.0)])
pf_VO_AA_OOtBu = PartFun(nma_VO_AA_OOtBu, [ExtTrans(), ExtRot(1),Vibrations(freq_scaling=1.0, zp_scaling=1.0)])
pf_TBHP = PartFun(nma_TBHP, [ExtTrans(), ExtRot(1),Vibrations(freq_scaling=1.0, zp_scaling=1.0)])
pf_TS =  PartFun(nma_TS, [ExtTrans(), ExtRot(1),Vibrations(freq_scaling=1.0, zp_scaling=1.0)])
pf_cyclohexene = PartFun(nma_cyclohexene, [ExtTrans(), ExtRot(1),Vibrations(freq_scaling=1.0, zp_scaling=1.0)])

# One can compute the equilibrium coefficient for the activation reaction and the forward coefficient for the reaction separately:

# Thermodynamic model
tm = ThermodynamicModel([pf_VO_AA_OH_H2O,pf_TBHP],[pf_VO_AA_OOtBu,pf_H2O,pf_H2O])
print "Thermodynamic model:"
print "equilibrium constant and gibbs free energy difference"
print tm.equilibrium_constant(323.0,do_log=False)/tm.unit, tm.unit_name,"    " ,tm.free_energy_change(323.0)/kjmol, "kJ/mol"
print
# Kinetic model
km = KineticModel([pf_VO_AA_OOtBu,pf_cyclohexene],pf_TS,tunneling=None)
print "Kinetic model:"
print "rate constant and gibbs free energy difference"
print km.rate_constant(323.0,do_log=False)/km.unit, "%s" % km.unit_name,"   " , km.free_energy_change(323.0)/kjmol, "kJ/mol"
print
ra = ReactionAnalysis(km, 273, 373.1, temp_step=10)
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


# Or one can determine a global k value and fit this within a temperature region by making use of the class ActivationKineticModel
akm = ActivationKineticModel(tm,km)
rakm = ReactionAnalysis(akm,273,373.1,temp_step=10)
rakm.plot_arrhenius("arrhenius_akm.png")
rakm.monte_carlo()
rakm.plot_parameters("parameters_akm.png")
print "Activation kinetic model:"
print "global rate constant and global gibbs free energy difference"
print akm.rate_constant(323.0)/akm.unit, "%s" % akm.unit_name,"    ", akm.free_energy_change(323.0)/kjmol, "kJ/mol"
akm.write_to_file("activation_model.txt")
