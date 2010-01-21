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
# Import units
from molmod.units import kjmol
# Import pylab for plotting
import pylab
# Import standard python libraries
import os, sys


def load_rotor(mol, filename, rotsym, even, expansion=5):
    dihedral, angles, energies, geometries, top_indexes = load_rotscan_g03(filename)
    cancel_freq = compute_cancel_frequency(mol, dihedral, top_indexes)
    rotor = Rotor(
        dihedral, top_indexes, cancel_freq, rotsym=rotsym, even=even,
        potential=(angles, energies, expansion), num_levels=50
    )
    return rotor


def run(do_rotor, do_counterpoise, load_sp):
    prefix = {True: "ir", False: "ho"}[do_rotor]
    prefix += {True: "_cps", False: "_bss"}[do_counterpoise]

    if load_sp:
        mol_ethyl = load_molecule_g03fchk("ethyl__freq/gaussian.fchk", "ethyl__sp/gaussian.fchk")
        mol_ethene = load_molecule_g03fchk("ethene__freq/gaussian.fchk", "ethene__sp/gaussian.fchk")
    else:
        mol_ethyl = load_molecule_g03fchk("ethyl__freq/gaussian.fchk")
        mol_ethene = load_molecule_g03fchk("ethene__freq/gaussian.fchk")
    if do_counterpoise:
        mol_ts_gauche = load_molecule_g03fchk("ts_ad1_gauche__freq/gaussian.fchk", "ts_ad1_gauche__bsse/gaussian.fchk")
        mol_ts_trans = load_molecule_g03fchk("ts_ad1_trans__freq/gaussian.fchk", "ts_ad1_trans__bsse/gaussian.fchk")
    else:
        if load_sp:
            mol_ts_gauche = load_molecule_g03fchk("ts_ad1_gauche__freq/gaussian.fchk", "ts_ad1_gauche__sp/gaussian.fchk")
            mol_ts_trans = load_molecule_g03fchk("ts_ad1_trans__freq/gaussian.fchk", "ts_ad1_gauche__sp/gaussian.fchk")
        else:
            mol_ts_gauche = load_molecule_g03fchk("ts_ad1_gauche__freq/gaussian.fchk")
            mol_ts_trans = load_molecule_g03fchk("ts_ad1_trans__freq/gaussian.fchk")
    # Perform normal mode analysis on the molecules
    nma_ethyl = NMA(mol_ethyl, ConstrainExt(gradient_threshold=1e-3))
    nma_ethene = NMA(mol_ethene, ConstrainExt(gradient_threshold=1e-3))
    nma_ts_gauche = NMA(mol_ts_gauche, ConstrainExt(gradient_threshold=1e-3))
    nma_ts_trans = NMA(mol_ts_trans, ConstrainExt(gradient_threshold=1e-3))
    if do_rotor:
        # Construct the rotors
        rotor_ethyl = load_rotor(mol_ethyl, "ethyl__scan_methyl/gaussian.log", 6, True, 1)
        rotor1_ts_gauche = load_rotor(mol_ts_gauche, "ts_ad1_gauche__scan_methyl/gaussian.log", 3, False)
        rotor2_ts_gauche = load_rotor(mol_ts_gauche, "ts_ad1_trans__scan_forming_bond/gaussian.log", 1, True)
        rotor1_ts_trans = load_rotor(mol_ts_trans, "ts_ad1_trans__scan_methyl/gaussian.log", 3, False)
        rotor2_ts_trans = load_rotor(mol_ts_trans, "ts_ad1_trans__scan_forming_bond/gaussian.log", 1, True)
        # Construct the partition functions.
        pf_ethyl = PartFun(nma_ethyl, [
            ExternalTranslation(), ExternalRotation(1), Electronic(2),
            #Vibrations(freq_scaling=0.9614, zp_scaling=0.9806),
            rotor_ethyl,
        ])
        pf_ethene = PartFun(nma_ethene, [
            ExternalTranslation(), ExternalRotation(4),
            #Vibrations(freq_scaling=0.9614, zp_scaling=0.9806),
        ])
        pf_ts_gauche = PartFun(nma_ts_gauche, [
            ExternalTranslation(), ExternalRotation(1), Electronic(2),
            #Vibrations(freq_scaling=0.9614, zp_scaling=0.9806),
            rotor1_ts_gauche, rotor2_ts_gauche,
        ])
        pf_ts_trans = PartFun(nma_ts_trans, [
            ExternalTranslation(), ExternalRotation(1), Electronic(2),
            #Vibrations(freq_scaling=0.9614, zp_scaling=0.9806),
            rotor1_ts_trans, rotor2_ts_trans,
        ])
    else:
        # Construct the partition functions.
        pf_ethyl = PartFun(nma_ethyl, [
            ExternalTranslation(), ExternalRotation(1), Electronic(2),
            #Vibrations(freq_scaling=0.9614, zp_scaling=0.9806),
        ])
        pf_ethene = PartFun(nma_ethene, [
            ExternalTranslation(), ExternalRotation(4),
            #Vibrations(freq_scaling=0.9614, zp_scaling=0.9806),
        ])
        pf_ts_gauche = PartFun(nma_ts_gauche, [
            ExternalTranslation(), ExternalRotation(1), Electronic(2),
            #Vibrations(freq_scaling=0.9614, zp_scaling=0.9806),
        ])
        pf_ts_trans = PartFun(nma_ts_trans, [
            ExternalTranslation(), ExternalRotation(1), Electronic(2),
            #Vibrations(freq_scaling=0.9614, zp_scaling=0.9806),
        ])

    if do_rotor:
        # Plot the energy levels and the potential of the hindered rotor. The
        # temperature argument is used to indicate the population of each level in the
        # plot.
        rotor_ethyl.plot_levels("rotor_ethyl_energy_levels.png", 300)
        rotor1_ts_gauche.plot_levels("rotor1_ts_gauche_energy_levels.png", 300)
        rotor2_ts_gauche.plot_levels("rotor2_ts_gauche_energy_levels.png", 300)
        rotor1_ts_trans.plot_levels("rotor1_ts_trans_energy_levels.png", 300)
        rotor2_ts_trans.plot_levels("rotor2_ts_trans_energy_levels.png", 300)


    # Analyse the partition functions in detail
    ta_ethyl = ThermoAnalysis(pf_ethyl, [300,400,500,600])
    ta_ethyl.write_to_file("%s_thermo_ethyl.csv" % prefix)
    ta_ethene = ThermoAnalysis(pf_ethene, [300,400,500,600])
    ta_ethene.write_to_file("%s_thermo_ethene.csv" % prefix)
    ta_ts_gauche = ThermoAnalysis(pf_ts_gauche, [300,400,500,600])
    ta_ts_gauche.write_to_file("%s_thermo_ts_gauche.csv" % prefix)
    ta_ts_trans = ThermoAnalysis(pf_ts_trans, [300,400,500,600])
    ta_ts_trans.write_to_file("%s_thermo_ts_trans.csv" % prefix)

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
    pylab.savefig("%s_arrhenius.png" % prefix)

    # Estimate the error on the kinetic parameters due to level of theory artifacts
    # with Monte Carlo sampling. The monte_carlo method takes three optional
    # arguments:
    #  1) freq_error: the relative systematic error on the frequencies
    #  2) freq_energy: the relative error on the energy
    #  4) num_iter: the number of monte carlo samples
    ra_gauche.monte_carlo(0.05, 0.05, 100)
    ra_trans.monte_carlo(0.05, 0.05, 100)
    # plot the parameters, this includes the monte carlo results
    pylab.clf()
    ra_gauche.plot_parameters(label="gauche", color="red")
    ra_trans.plot_parameters(label="trans", color="blue")
    pylab.legend(loc=0)
    pylab.savefig("%s_parameters.png" % prefix)
    # write all results to a file.
    ra_gauche.write_to_file("%s_reaction_gauche.txt" % prefix)
    ra_trans.write_to_file("%s_reaction_trans.txt" % prefix)

    def write_ra_summary(fn, ra):
        f = file(fn, "w")
        print >> f, "% 10.5e % 10.5e % 10.5e % 10.5e    %10.5e %10.2e    %10.2e %10.2e" % (
            ra.compute_rate_coeff(300)/ra.unit,
            ra.compute_rate_coeff(400)/ra.unit,
            ra.compute_rate_coeff(500)/ra.unit,
            ra.compute_rate_coeff(600)/ra.unit,
            ra.A/ra.unit,
            ra.Ea/kjmol,
            ra.compute_delta_G(0.0)/kjmol,
            ra.compute_delta_E()/kjmol,
        )
        f.close()

    write_ra_summary("%s_summary_gauche.txt" % prefix, ra_gauche)
    write_ra_summary("%s_summary_trans.txt" % prefix, ra_trans)


usage = """USAGE: ./reaction.py dirname

Analyzes the kinetics of the addition of ethene to ethyl with and without
hindered rotors. It looks at both the trans and gauce pathways.
"""


def main():
    from optparse import OptionParser
    parser = OptionParser(usage)
    options, args = parser.parse_args()

    if len(args) != 1:
        parser.error("Expecting exactly on argument")

    print "Processing %s" % args[0]
    os.chdir(args[0])
    load_sp = args[0].startswith("GEO")
    for do_rotor in True, False:
        for do_counterpoise in True, False:
            try:
                run(do_rotor, do_counterpoise, load_sp)
                print "      OK: do_rotor=%i, do_counterpoise=%i, load_sp=%i" % (do_rotor, do_counterpoise, load_sp)
            except (IOError, KeyError), e:
                print "  Failed: do_rotor=%i, do_counterpoise=%i, load_sp=%i" % (do_rotor, do_counterpoise, load_sp)
                print e


if __name__ == "__main__":
    main()

