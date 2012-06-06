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
mol_oxygen = load_molecule_g03fchk("oxygen.fchk")
mol_hydrogen = load_molecule_g03fchk("hydrogen.fchk")
mol_water = load_molecule_g03fchk("water.fchk")
# Perform normal mode analysis on the molecules
nma_oxygen = NMA(mol_oxygen, ConstrainExt())
nma_hydrogen = NMA(mol_hydrogen, ConstrainExt())
nma_water = NMA(mol_water, ConstrainExt())
# Construct the partition functions.
pf_oxygen = PartFun(nma_oxygen, [ExtTrans(), ExtRot()])
pf_hydrogen = PartFun(nma_hydrogen, [ExtTrans(), ExtRot()])
pf_water = PartFun(nma_water, [ExtTrans(), ExtRot()])
# Define a kinetic model for the chemical reaction.
tm = ThermodynamicModel([pf_oxygen, (pf_hydrogen, 2)], [(pf_water, 2)])
# Write tables with the principal energies at 298.15K
tm.write_table(298.15, "formation.csv")
# Write an overview of the thermodynamic model to a file
tm.write_to_file("formation.txt")
