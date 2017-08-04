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
# import matplotlib.pyplot as pt for plotting
import matplotlib.pyplot as pt


# Define an auxiliary function
def load_rotor(mol, filename, rotsym, even, top_indexes=None):
    rot_scan = load_rotscan_g03log(filename, top_indexes)
    rotor = Rotor(rot_scan, mol, rotsym=rotsym, even=even)
    return rotor

mol_ethyl = load_molecule_g03fchk("ethyl/freq/gaussian.fchk")
mol_ethene = load_molecule_g03fchk("ethene/freq/gaussian.fchk")
mol_ts_gauche = load_molecule_g03fchk("ts_ad1/freq_gauche/gaussian.fchk")
mol_ts_trans = load_molecule_g03fchk("ts_ad1/freq_trans/gaussian.fchk")
# Perform normal mode analysis on the molecules
nma_ethyl = NMA(mol_ethyl, ConstrainExt())
nma_ethene = NMA(mol_ethene, ConstrainExt())
nma_ts_gauche = NMA(mol_ts_gauche, ConstrainExt())
nma_ts_trans = NMA(mol_ts_trans, ConstrainExt())
# Construct the rotors
rotor_ethyl = load_rotor(mol_ethyl, "ethyl/scan/gaussian.log", 6, True)
rotor1_ts_gauche = load_rotor(mol_ts_gauche, "ts_ad1/scan1/gaussian.log", 3, False, [4,5,6])
rotor2_ts_gauche = load_rotor(mol_ts_gauche, "ts_ad1/scan2/gaussian.log", 1, True)
rotor1_ts_trans = load_rotor(mol_ts_trans, "ts_ad1/scan1/gaussian.log", 3, False, [4,5,6])
rotor2_ts_trans = load_rotor(mol_ts_trans, "ts_ad1/scan2/gaussian.log", 1, True)
# Construct the partition functions.
pf_ethyl = PartFun(nma_ethyl, [
    ExtTrans(), ExtRot(), Electronic(2),
    Vibrations(freq_scaling=0.9614, zp_scaling=0.9806),
    rotor_ethyl,
])
pf_ethene = PartFun(nma_ethene, [
    ExtTrans(), ExtRot(),
    Vibrations(freq_scaling=0.9614, zp_scaling=0.9806),
])
pf_ts_gauche = PartFun(nma_ts_gauche, [
    ExtTrans(), ExtRot(), Electronic(2),
    Vibrations(freq_scaling=0.9614, zp_scaling=0.9806),
    rotor1_ts_gauche, rotor2_ts_gauche,
])
pf_ts_trans = PartFun(nma_ts_trans, [
    ExtTrans(), ExtRot(), Electronic(2),
    Vibrations(freq_scaling=0.9614, zp_scaling=0.9806),
    rotor1_ts_trans, rotor2_ts_trans,
])

# Analyse the partition functions in detail
ta_ethyl = ThermoAnalysis(pf_ethyl, [300,400,500,600])
ta_ethyl.write_to_file("thermo_ethyl.csv")
ta_ethene = ThermoAnalysis(pf_ethene, [300,400,500,600])
ta_ethene.write_to_file("thermo_ethene.csv")
ta_ts_gauche = ThermoAnalysis(pf_ts_gauche, [300,400,500,600])
ta_ts_gauche.write_to_file("thermo_ts_gauche.csv")
ta_ts_trans = ThermoAnalysis(pf_ts_trans, [300,400,500,600])
ta_ts_trans.write_to_file("thermo_ts_trans.csv")

# Define kinetic models for the chemical reaction. These are the mandatory arguments:
#  1) a list of reactant partition functions
#     (one for unimolecular, two for bimolecular, ...)
#  2) the transition state partition function
# There is one more optional argument. (not used here).
#  3) tunneling: a model for the tunelling correction
km_gauche = KineticModel([pf_ethyl, pf_ethene], pf_ts_gauche)
km_trans = KineticModel([pf_ethyl, pf_ethene], pf_ts_trans)

# Analyze the chemical reactions. These are the arguments:
#  1) A kinetic model
#  2) the starting temperature for the fit
#  3) the final temperature for the fit
# The following argument is optional:
#  4) temp_step: The interval on the temperature grid in Kelvin, 10 is default
ra_gauche = ReactionAnalysis(km_gauche, 300, 600)
ra_trans = ReactionAnalysis(km_trans, 300, 600)

# make the Arrhenius plots
pt.clf()
ra_gauche.plot_arrhenius(label="gauche", color="red")
ra_trans.plot_arrhenius(label="trans", color="blue")
pt.legend(loc=0)
pt.savefig("arrhenius.png")

# Estimate the error on the kinetic parameters due to level of theory artifacts
# with Monte Carlo sampling. The monte_carlo method takes three optional
# arguments:
#  1) freq_error: the absolute stochastic error on the frequencies (default=1*invcm)
#  2) energy_error: the absolute error on the energy (default=0.0)
#  3) num_iter: the number of monte carlo samples (default=100)
ra_gauche.monte_carlo(num_iter=10)
ra_trans.monte_carlo(num_iter=10)
# plot the parameters, this includes the monte carlo results
pt.clf()
ra_gauche.plot_parameters(label="gauche", color="red")
ra_trans.plot_parameters(label="trans", color="blue", marker="^")
pt.legend(loc=0, numpoints=1)
pt.xlabel("E$_a$ [kJ mol$^{-1}$]")
pt.ylabel("ln(A) [ln(m$^3$ mol$^{-1}$ s$^{-1}$)]")
pt.xlim(33.24,33.75)
pt.savefig("parameters.png")
# write all results to a file.
ra_gauche.write_to_file("reaction_gauche.txt")
ra_trans.write_to_file("reaction_trans.txt")
# Plot the energy levels and the potential of the hindered rotor. The
# temperature argument is used to indicate the population of each level in the
# plot.
rotor_ethyl.plot_levels("rotor_ethyl_energy_levels.png", 300)
rotor1_ts_gauche.plot_levels("rotor1_ts_gauche_energy_levels.png", 300)
rotor2_ts_gauche.plot_levels("rotor2_ts_gauche_energy_levels.png", 300)
rotor1_ts_trans.plot_levels("rotor1_ts_trans_energy_levels.png", 300)
rotor2_ts_trans.plot_levels("rotor2_ts_trans_energy_levels.png", 300)
