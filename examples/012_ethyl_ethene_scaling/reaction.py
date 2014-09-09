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


# load the molecules
mol_ethyl = load_molecule_g03fchk("ethyl/freq/gaussian.fchk")
mol_ethene = load_molecule_g03fchk("ethene/freq/gaussian.fchk")
mol_ts_gauche = load_molecule_g03fchk("ts_ad1/freq_gauche/gaussian.fchk")

# Perform normal mode analysis on the molecules
nma_ethyl = NMA(mol_ethyl, ConstrainExt())
nma_ethene = NMA(mol_ethene, ConstrainExt())
nma_ts_gauche = NMA(mol_ts_gauche, ConstrainExt())

cases = [ #(freq_scaling, color)
   (0.95, "blue"),
   (1.0, "purple"),
   (1.05, "red"),
   (1.10, "orange"),
   (1.15, "green"),
]

# start with clear figure:
pt.clf()

for s, color in cases:
    # For each scaling factor, a new curve is plotted on the Arrhenius plot

    # Construct the partition functions.
    pf_ethyl = PartFun(nma_ethyl, [ExtTrans(), ExtRot(), Electronic(2), Vibrations(freq_scaling=s)])
    pf_ethene = PartFun(nma_ethene, [ExtTrans(), ExtRot(), Vibrations(freq_scaling=s)])
    pf_ts_gauche = PartFun(nma_ts_gauche, [ExtTrans(), ExtRot(), Electronic(2), Vibrations(freq_scaling=s)])

    # Define a kinetic model for the chemical reaction. These are the mandatory arguments:
    #  1) a list of reactant partition functions
    #     (one for unimolecular, two for bimolecular, ...)
    #  2) the transition state partition function
    # There is one more optional argument.
    #  3) tunneling: a model for the tunelling correction
    km_gauche = KineticModel([pf_ethyl, pf_ethene], pf_ts_gauche)

    # Analyze the chemical reaction. These are the arguments:
    #  1) A kinetic model
    #  2) the starting temperature for the fit
    #  3) the final temperature for the fit
    # The following argument is optional:
    #  4) temp_step: The interval on the temperature grid in Kelvin, 10 is default
    ra_gauche = ReactionAnalysis(km_gauche, 300, 600)

    # Add a line to the Arrhenius plots
    ra_gauche.plot_arrhenius(label="gauche %.2f" % s, color=color)

# Add a legend (loc=0 means that the legend is put on an empty place in the plot)
pt.legend(loc=0)
# Write the figure to file.
pt.savefig("arrhenius.png")
