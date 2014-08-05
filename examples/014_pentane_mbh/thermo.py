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


# Import the tamkin libarary.
from tamkin import *
# Import unit conversin factors
from molmod import *
# Import numpy stuff
from numpy import *


# Load the gaussian data.
molecule = load_molecule_g03fchk("gaussian.fchk")
print "Energy [kJ/mol]:", molecule.energy/kjmol


# A) Detect all blocks such that only the dihedral angles remain free, and apply
# MBH
# detect the graph based on inter-atomic distances
molecule.set_default_graph()
# define a block for each atom and its neighbors (if it has more than one neighbor)
blocks = []
for central, neighbors in molecule.graph.neighbors.iteritems():
    if len(neighbors) > 1:
        block = list(neighbors)
        block.append(central)
        blocks.append(block)
print "Blocks:", blocks
nma = NMA(molecule, MBH(blocks))
print "The zero eigenmodes: %s" % nma.zeros
# write the modes to file
dump_modes_gaussian("molden.vib", nma)

# B) Construct a partition function with the typical gas phase contributions.
pf = PartFun(nma, [ExtTrans(), ExtRot()])
print "Heat capacity at 300K, constant pressure [J/(mol*K)]:", pf.heat_capacity(300*kelvin)/(joule/mol/kelvin)
# Write some general information about the molecule and the partition function
# to a file.
pf.write_to_file("partfun.txt")
# Write an extensive overview of the thermodynamic properties to a file:
ta = ThermoAnalysis(pf, arange(300,601,5))
ta.write_to_file("thermo.csv")
