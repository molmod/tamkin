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
pylab.clf()

for s, color in cases:
   # For each scaling factor, a new curve is plotted on the Arrhenius plot

   # Construct the partition functions.
   pf_ethyl = PartFun(nma_ethyl, [ExtTrans(), ExtRot(), Electronic(2), Vibrations(freq_scaling=s)])
   pf_ethene = PartFun(nma_ethene, [ExtTrans(), ExtRot(), Vibrations(freq_scaling=s)])
   pf_ts_gauche = PartFun(nma_ts_gauche, [ExtTrans(), ExtRot(), Electronic(2), Vibrations(freq_scaling=s)])

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

   # Add a line to the Arrhenius plots
   ra_gauche.plot_arrhenius(label="gauche %.2f" % s, color=color)

# Add a legend (loc=0 means that the legend is put on an empty place in the plot)
pylab.legend(loc=0)
# Write the figure to file.
pylab.savefig("arrhenius.png")


