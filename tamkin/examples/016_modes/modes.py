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


# Import the tamkin libarary
from tamkin import *
# Import other libraries
import numpy as np # numerical and array routines
from molmod import * # for units

# Load the gamess data
molecule = load_molecule_pcgamess_punch("PUNCH")
# Perform the normal mode analysis
#nma = NMA(molecule, do_modes=True)
nma = NMA(molecule, ConstrainExt(), do_modes=True)

# Print modes and wavenumbers on screen
invcm = lightspeed/centimeter
for i in xrange(len(nma.freqs)):
    print "Mode:", i
    print "Wavenumber:", nma.freqs[i]/invcm
    print "Is zero:", i in nma.zeros
    print "Norm of eigenmode:", np.linalg.norm(nma.modes[:,i])
    print "Components:"
    print nma.modes[:,i]
    print
