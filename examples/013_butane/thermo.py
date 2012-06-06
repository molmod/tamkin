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
molecule = load_molecule_g03fchk("freq/gaussian.fchk")
print "Energy [kJ/mol]:", molecule.energy/kjmol
rot_scan = load_rotscan_g03log("scan/gaussian.log")


# A) Perform a standard normal mode analysis
nma = NMA(molecule, ConstrainExt())
rotor = Rotor(rot_scan, molecule, rotsym=1, even=True)
pf = PartFun(nma, [ExtTrans(), ExtRot(), rotor])
print "Heat capacity at 300K, constant pressure [J/(mol*K)]:", pf.heat_capacity(300*kelvin)/(joule/mol/kelvin)
# Write some general information about the molecule and the partition function
# to a file.
pf.write_to_file("partfun.txt")
rotor.plot_levels("energy_levels.png", 300, do_levels=False, do_data=False)

I = rotor.reduced_moment
for angle, energy in array(rot_scan.potential).transpose():
    K = rotor.hb.eval_deriv2(array([angle]), rotor.v_coeffs)[0]
    nu = sqrt(abs(K/I))
    hertz = 1/second
    print "Angle in deg:", angle/deg
    print "Energy in kJ/mol", energy/joule
    print "Frequency in THz:", nu/hertz/1e12/(2*pi)
    print "Wavenumber in cm-1:", nu/(lightspeed/centimeter)/(2*pi)
    print
