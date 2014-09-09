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
# Import units
from molmod.units import kjmol
# import matplotlib.pyplot as pt for plotting
import matplotlib.pyplot as pt
# Import standard python libraries
import os, sys


def load_rotor(mol, filename, rotsym, even, expansion=5):
    rot_scan = load_rotscan_g03log(filename)
    rotor = Rotor(rot_scan, mol, rotsym=rotsym, even=even)
    return rotor


def load_cps_barrier(fn_template):
    """Loads the counterpoise corrected barrier from five sp calculations"""
    from molmod.io import FCHKFile
    def load_ener(fn_fchk):
        fchk = FCHKFile(fn_fchk, field_labels=["Total Energy"])
        return fchk.fields["Total Energy"]
    return (
        load_ener(fn_template % "full_0") +
        load_ener(fn_template % "sole_1") - load_ener(fn_template % "full_1") +
        load_ener(fn_template % "sole_2") - load_ener(fn_template % "full_2")
    )



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
        try:
            mol_ts_gauche = load_molecule_g03fchk("ts_ad1_gauche__freq/gaussian.fchk", "ts_ad1_gauche__bsse/gaussian.fchk")
        except (IOError, KeyError):
            cps_barrier_gauche = load_cps_barrier("ts_ad1_gauche__cps_%s/gaussian.fchk")
            mol_ts_gauche = load_molecule_g03fchk("ts_ad1_gauche__freq/gaussian.fchk", energy=cps_barrier_gauche)
        try:
            mol_ts_trans = load_molecule_g03fchk("ts_ad1_trans__freq/gaussian.fchk", "ts_ad1_trans__bsse/gaussian.fchk")
        except (IOError, KeyError):
            cps_barrier_trans = load_cps_barrier("ts_ad1_trans__cps_%s/gaussian.fchk")
            mol_ts_trans = load_molecule_g03fchk("ts_ad1_trans__freq/gaussian.fchk", energy=cps_barrier_trans)
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
            ExtTrans(), ExtRot(), Electronic(2),
            #Vibrations(freq_scaling=0.9614, zp_scaling=0.9806),
            rotor_ethyl,
        ])
        pf_ethene = PartFun(nma_ethene, [
            ExtTrans(), ExtRot(),
            #Vibrations(freq_scaling=0.9614, zp_scaling=0.9806),
        ])
        pf_ts_gauche = PartFun(nma_ts_gauche, [
            ExtTrans(), ExtRot(), Electronic(2),
            #Vibrations(freq_scaling=0.9614, zp_scaling=0.9806),
            rotor1_ts_gauche, rotor2_ts_gauche,
        ])
        pf_ts_trans = PartFun(nma_ts_trans, [
            ExtTrans(), ExtRot(), Electronic(2),
            #Vibrations(freq_scaling=0.9614, zp_scaling=0.9806),
            rotor1_ts_trans, rotor2_ts_trans,
        ])
    else:
        # Construct the partition functions.
        pf_ethyl = PartFun(nma_ethyl, [
            ExtTrans(), ExtRot(1), Electronic(2),
            #Vibrations(freq_scaling=0.9614, zp_scaling=0.9806),
        ])
        pf_ethene = PartFun(nma_ethene, [
            ExtTrans(), ExtRot(4),
            #Vibrations(freq_scaling=0.9614, zp_scaling=0.9806),
        ])
        pf_ts_gauche = PartFun(nma_ts_gauche, [
            ExtTrans(), ExtRot(1), Electronic(2),
            #Vibrations(freq_scaling=0.9614, zp_scaling=0.9806),
        ])
        pf_ts_trans = PartFun(nma_ts_trans, [
            ExtTrans(), ExtRot(1), Electronic(2),
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

    # Define kinetic models for the chemical reaction. These are the mandatory arguments:
    #  1) a list of reactant partition functions
    #     (one for unimolecular, two for bimolecular, ...)
    #  2) the transition state partition function
    # There are two more optional arguments
    #  3) cp: model at constant pressure, default=True
    #  4) tunneling: a model for the tunelling correction
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
    pt.savefig("%s_arrhenius.png" % prefix)

    # Estimate the error on the kinetic parameters due to level of theory artifacts
    # with Monte Carlo sampling. The monte_carlo method takes three optional
    # arguments:
    #  1) freq_error: the absolute stochastic error on the frequencies (default=1*invcm)
    #  2) energy_error: the absolute error on the energy (default=0.0)
    #  3) num_iter: the number of monte carlo samples (default=100)
    ra_gauche.monte_carlo(num_iter=100)
    ra_trans.monte_carlo(num_iter=100)
    # plot the parameters, this includes the monte carlo results
    pt.clf()
    ra_gauche.plot_parameters(label="gauche", color="red")
    ra_trans.plot_parameters(label="trans", color="blue")
    pt.legend(loc=0)
    pt.savefig("%s_parameters.png" % prefix)
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
            #try:
                run(do_rotor, do_counterpoise, load_sp)
                print "      OK: do_rotor=%i, do_counterpoise=%i, load_sp=%i" % (do_rotor, do_counterpoise, load_sp)
            #except (IOError, KeyError), e:
            #    print "  Failed: do_rotor=%i, do_counterpoise=%i, load_sp=%i" % (do_rotor, do_counterpoise, load_sp)
            #    print e


if __name__ == "__main__":
    main()
