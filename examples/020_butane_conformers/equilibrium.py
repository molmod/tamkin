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
mol_trans = load_molecule_g03fchk("trans.fchk")
mol_gauche = load_molecule_g03fchk("gauche.fchk")
# Perform normal mode analysis on the molecules
nma_trans = NMA(mol_trans, ConstrainExt())
nma_gauche = NMA(mol_gauche, ConstrainExt())
# Construct the partition functions.
pf_trans = PartFun(nma_trans, [ExtTrans(), ExtRot()])
pf_gauche = PartFun(nma_gauche, [ExtTrans(), ExtRot(), Electronic(multiplicity=2)])
# Define a kinetic model for the chemical reaction.
tm = ThermodynamicModel([pf_trans], [pf_gauche])
# Write tables with the principal energies at 300K, 400K, 500K and 600K
tm.write_table(300, "equilibrium300.csv")
tm.write_table(400, "equilibrium400.csv")
tm.write_table(500, "equilibrium500.csv")
tm.write_table(600, "equilibrium600.csv")
# Write an overview of the thermodynamic model to a file
tm.write_to_file("equilibrium.txt")
