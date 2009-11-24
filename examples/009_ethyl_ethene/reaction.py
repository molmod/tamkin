#!/usr/bin/env python
# TAMkin is a post-processing toolkit for thermochemistry and kinetics analysis.
# Copyright (C) 2008-2009 Toon Verstraelen <Toon.Verstraelen@UGent.be>,
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


# Define an auxiliary function
def load_rotor(mol, filename, rotsym, even):
    dihedral, angles, energies, geometries, top_indexes = load_rotscan_g03(filename)
    cancel_freq = compute_cancel_frequency(mol, top_indexes)
    rotor = Rotor(
        top_indexes, cancel_freq, rotsym=rotsym, even=even,
        potential=(angles, energies, 5), num_levels=50
    )
    return rotor

mol_ethyl = load_molecule_g03fchk("ethyl/freq/gaussian.fchk")
mol_ethene = load_molecule_g03fchk("ethene/freq/gaussian.fchk")
mol_trans = load_molecule_g03fchk("trans/freq2/gaussian.fchk")
# Perform normal mode analysis on the three molecules
nma_ethyl = NMA(mol_ethyl, ConstrainExt())
nma_ethene = NMA(mol_ethene, ConstrainExt())
nma_trans = NMA(mol_trans, ConstrainExt())
# Construct the rotor about the forming bond in the transition state
rotor_ethyl = load_rotor(mol_ethyl, "ethyl/scan/gaussian.log", 6, True)
rotor_trans1 = load_rotor(mol_trans, "trans/scan1/gaussian.log", 3, False)
rotor_trans2 = load_rotor(mol_trans, "trans/scan2/gaussian.log", 1, True)
# Construct the three partition functions.
pf_ethyl = PartFun(nma_ethyl, [
    ExternalTranslation(), ExternalRotation(1), Electronic(2),
    Vibrations(freq_scaling=0.9614, zp_scaling=0.9806),
    rotor_ethyl
])
pf_ethene = PartFun(nma_ethene, [
    ExternalTranslation(), ExternalRotation(4),
    Vibrations(freq_scaling=0.9614, zp_scaling=0.9806),
])
pf_trans = PartFun(nma_trans, [
    ExternalTranslation(), ExternalRotation(1), Electronic(2),
    Vibrations(freq_scaling=0.9614, zp_scaling=0.9806),
    rotor_trans1, rotor_trans2
])

# Analyse the partition functions in detail
ta_ethyl = ThermoAnalysis(pf_ethyl, [300,400,500,600])
ta_ethyl.write_to_file("thermo_ethyl.csv")
ta_ethene = ThermoAnalysis(pf_ethene, [300,400,500,600])
ta_ethene.write_to_file("thermo_ethene.csv")
ta_trans = ThermoAnalysis(pf_trans, [300,400,500,600])
ta_trans.write_to_file("thermo_trans.csv")

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
ra = ReactionAnalysis([pf_ethyl, pf_ethene], pf_trans, 300, 600, 10)
ra.plot("arrhenius.png") # make the Arrhenius plot

# Estimate the error on the kinetic parameters due to level of theory artifacts
# with Monte Carlo sampling. The monte_carlo method takes three optional
# arguments:
#  1) freq_error: the relative systematic error on the frequencies
#  2) freq_energy: the relative error on the energy
#  4) num_iter: the number of monte carlo samples
ra.monte_carlo(0.05, 0.00, 100)
# plot the parameters, this includes the monte carlo results
ra.plot_parameters("parameters.png")
# write all results to a file.
ra.write_to_file("reaction.txt") # summary of the analysis
# Plot the energy levels and the potential of the hindered rotor. The
# temperature argument is used to indicate the population of each level in the
# plot.
rotor_ethyl.plot_levels("rotor_ethyl_energy_levels.png", 300)
rotor_trans1.plot_levels("rotor_trans1_energy_levels.png", 300)
rotor_trans2.plot_levels("rotor_trans2_energy_levels.png", 300)



