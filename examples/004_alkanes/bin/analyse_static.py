#!/usr/bin/python
# TAMkin is a post-processing toolkit for thermochemistry and kinetics analysis.
# Copyright (C) 2008-2009 Toon Verstraelen <Toon.Verstraelen@UGent.be>,
# Matthias Vandichel <Matthias.Vandichel@UGent.be> and
# An Ghysels <An.Ghysels@UGent.be>
#
# This file is part of TAMkin.
#
# TAMkin is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
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


# Import the tamkin library
from tamkin import *
# Import unit conversion factors and physical constants
from molmod.units import *
from molmod.constants import *

import sys, os, numpy # standard libraries


# Parse the command line arguments
args = sys.argv[1:]
if len(args) != 2:
    print "Two arguments are required: alkane_n and the symmetry number"
    sys.exit()
symmetry_number = int(args[1])

# Hardcoded parameters
temp = 298.15
pressure = 1*atm

# Load the cp2k data and construct a partition function object
molecule = load_molecule_cp2k(
    os.path.join(args[0], "opt.xyz"),
    os.path.join(args[0], "sp/sp.out"),
    os.path.join(args[0], "freq/freq.out"),
)
nma = NMA(molecule, ConstrainExt())
pf = PartFun(nma, [
    ExternalTranslation(FixedVolume(temp,pressure)),
    ExternalRotation(symmetry_number),
], Vibrations(classical=False)) # change this boolean to get classical vibrations


# Write the frequencies to a csv file
f = file(os.path.join(args[0], "freqs.csv"), "w")
print >> f, '"Frequency","Wavenumber","Vibrational temperature"'
print >> f, '"Atomic units","1/cm","K"'
for i in xrange(pf.vibrational.num_free):
    freq = pf.vibrational.freqs[i]
    print >> f, '%e,%f,%f' % (freq, freq/lightspeed*cm, 2*numpy.pi*freq/boltzmann)
f.close()


# B) Generate a thermo analysis report
ta = ThermoAnalysis(pf, [100, 200, 300, 400, 500, 600, 700, 800, 900, 100])
ta.write_to_file(os.path.join(args[0], "static_thermo.csv"))


