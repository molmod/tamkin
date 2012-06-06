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


# Import the tamkin library
from tamkin import *
# Import unit conversion factors and physical constants
from molmod.units import *
from molmod.constants import *

import sys, os, numpy # standard libraries


# Parse the command line arguments
args = sys.argv[1:]
if len(args) != 1:
    print "One arguments are required: alkane_n"
    sys.exit()

# Hardcoded parameters
temp = 298.15
pressure = 1*atm

# Load the cp2k data and construct a partition function object
molecule = load_molecule_cp2k(
    os.path.join(args[0], "opt.xyz"),
    os.path.join(args[0], "sp/sp.out"),
    os.path.join(args[0], "freq/freq.out"),
    is_periodic=False,
)
nma = NMA(molecule, ConstrainExt())
pf = PartFun(
    nma,
    [ExtTrans(gaslaw=IdealGasLaw(pressure)), ExtRot(), Vibrations(classical=False)],
)


# Write the frequencies to a csv file
f = file(os.path.join(args[0], "freqs.csv"), "w")
print >> f, '"Frequency","Wavenumber","Vibrational temperature"'
print >> f, '"Atomic units","1/cm","K"'
for i in xrange(len(pf.vibrational.freqs)):
    freq = pf.vibrational.freqs[i]
    print >> f, '%e,%f,%f' % (freq, freq/lightspeed*centimeter, 2*numpy.pi*freq/boltzmann)
f.close()


# B) Generate a thermo analysis report
ta = ThermoAnalysis(pf, [100, 200, 300, 400, 500, 600, 700, 800, 900, 100])
ta.write_to_file(os.path.join(args[0], "static_thermo.csv"))
