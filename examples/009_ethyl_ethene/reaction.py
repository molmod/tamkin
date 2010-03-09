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
# Import pylab for plotting
import pylab


# Define an auxiliary function
def load_rotor(mol, filename, rotsym, even):
    rot_scan = load_rotscan_g03log(filename)
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
rotor1_ts_gauche = load_rotor(mol_ts_gauche, "ts_ad1/scan1/gaussian.log", 3, False)
rotor2_ts_gauche = load_rotor(mol_ts_gauche, "ts_ad1/scan2/gaussian.log", 1, True)
rotor1_ts_trans = load_rotor(mol_ts_trans, "ts_ad1/scan1/gaussian.log", 3, False)
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

# Analyze the chemical reaction. These are the arguments:
#  1) a list of reactant partition functions
#     (one for unimolecular, two for bimolecular, ...)
#  2) the transition state partition function
#  3) the starting temperature for the fit
#  4) the final temperature for the fit
# The following are optional arguments:
#  6) temp_step: The interval on the temperature grid in Kelvin, 10 is default
#  5) mol_volume: a function that returns the molecular volume (inverse of the
#     concentration) for a given temperature. If not given, it assumes that
#     this is given: FixedVolume(temp=298.15*K, pressure=1*atm)
#  7) tunneling: a tunneling correction object
ra_gauche = ReactionAnalysis([pf_ethyl, pf_ethene], pf_ts_gauche, 300, 600, 10)
ra_trans = ReactionAnalysis([pf_ethyl, pf_ethene], pf_ts_trans, 300, 600, 10)

# make the Arrhenius plots
pylab.clf()
ra_gauche.plot_arrhenius(label="gauche", color="red")
ra_trans.plot_arrhenius(label="trans", color="blue")
pylab.legend(loc=0)
pylab.savefig("arrhenius.png")

# Estimate the error on the kinetic parameters due to level of theory artifacts
# with Monte Carlo sampling. The monte_carlo method takes three optional
# arguments:
#  1) freq_error: the absolute stochastic error on the frequencies (default=1*invcm)
#  2) energy_error: the absolute error on the energy (default=0.0)
#  3) num_iter: the number of monte carlo samples (default=100)
ra_gauche.monte_carlo(num_iter=1000)
ra_trans.monte_carlo(num_iter=1000)
# plot the parameters, this includes the monte carlo results
pylab.clf()
ra_gauche.plot_parameters(label="gauche", color="red")
ra_trans.plot_parameters(label="trans", color="blue", marker="^")
pylab.legend(loc=0, numpoints=1)
pylab.xlabel("E$_a$ [kJ mol$^{-1}$]")
pylab.ylabel("ln(A) [ln(m$^3$ mol$^{-1}$ s$^{-1}$)]")
pylab.xlim(33.24,33.75)
pylab.savefig("parameters.png")
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



