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


# Import libraries.
from tamkin import *   # The TAMkin library

# Load the molecules (including the Hessian etc.)
mol_ethyl = load_molecule_g03fchk("ethyl.fchk")
mol_ethene = load_molecule_g03fchk("ethene.fchk")
mol_ts_trans = load_molecule_g03fchk("ts_trans.fchk")
# Perform normal mode analysis on the molecules
nma_ethyl = NMA(mol_ethyl, ConstrainExt())
nma_ethene = NMA(mol_ethene, ConstrainExt())
nma_ts_trans = NMA(mol_ts_trans, ConstrainExt())
# Construct the partition functions.
pf_ethyl = PartFun(nma_ethyl, [ExtTrans(), ExtRot()])
pf_ethene = PartFun(nma_ethene, [ExtTrans(), ExtRot()])
pf_ts_trans = PartFun(nma_ts_trans, [ExtTrans(), ExtRot()])
# Define a kinetic model for the chemical reaction.
km_trans = KineticModel([pf_ethyl, pf_ethene], pf_ts_trans)
# Write tables with the principal energies at 300K, 400K, 500K and 600K
km_trans.write_table(300, "kinetic300.csv")
km_trans.write_table(400, "kinetic400.csv")
km_trans.write_table(500, "kinetic500.csv")
km_trans.write_table(600, "kinetic600.csv")
# Analyze the chemical reactions.
ra_trans = ReactionAnalysis(km_trans, 300, 600)
# Make the Arrhenius plot
ra_trans.plot_arrhenius("arrhenius.png")
# Write the analysis to a file
ra_trans.write_to_file("kinetic.txt")
