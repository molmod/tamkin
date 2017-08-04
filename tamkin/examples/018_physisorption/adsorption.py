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


import numpy as np

from tamkin import *
from molmod import *


# define some parameters
fixed = [6,7,8,9,10,11]
level = "m062x"
temp = 87.0*kelvin

# load molecules
mol_a = load_molecule_g03fchk("%s/argon/gaussian.fchk" % level)
mol_b = load_molecule_g03fchk("%s/benzene/freq/gaussian.fchk" % level)
mol_c = load_molecule_g03fchk("%s/complex/freq/gaussian.fchk" % level)

# do the normal mode analysis
nma_a = NMA(mol_a, ConstrainExt()) # does not much, but is mandatory due to API
nma_b = NMA(mol_b, PHVA(fixed)) # use PHVA to model the surface
nma_c = NMA(mol_c, PHVA(fixed)) # use PHVA to model the surface

# Add a small correction to the energy from the m062x computation as to get a
# perfect correspondence between the experimental and theoretical isotherm in
# the ideal gas limit There are several things that can cause errors in the
# m062x computation: poor graphite model (i.e. benzene), lack of BSSE
# correction, imperfections in the DFT functional. The first and the second
# error will probably compensate partially. The benzene approximation is rather
# dramatic because the Van der Waals sphere of Argon is larger than benzene.
nma_c.energy += -2.4*kjmol
# Note that kT at 87K is about 0.723 kJ/mol. This is example a very sensitive
# to the accuracy of the adsorption energy.

# setup an NpT partition function for 3D gas
pf_a = PartFun(nma_a, [ExtTrans(cp=True)])
# benzene acts a a substrate fixed in space, no translation or rotation
pf_b = PartFun(nma_b, [])
# setup an NVT partition function for the 2D gas. (The amount of surface is
# constant.) Specify that atom with index 12 (Argon) is the mobile one.
pf_c = PartFun(nma_c, [ExtTrans(cp=False, mobile=[12], dim=2)])
# create a model for the thermodynamic equilibrium
tm = ThermodynamicModel([pf_a, pf_b], [pf_c])
tm.write_to_file("adsorption_%s.txt" % level)
# reference data from http://dx.doi.org/10.1021/jp011745+, table 4
#   The argon cross section in the paper is 0.138 nm**2
#   This example is based on data at temperature of 87 K.
# first column: p/p0
#   with p0 = 0.96848*atm
# second column: alpha_s
#   Alpha_s is the ratio of the adsorbed amount over a reference, i.e. 3.108
#   cm**3 STP g**-1. The adsorbent has a specific surface area of 6.01 m**2
#   g**-1 and can be treated as a macroscopic graphite surface. We model this
#   surface as a graphene layer, which is further approximated by a single
#   benzene molecule. ouch.
ref_data = np.array([
    (1.28e-6, 9.17e-5),
    (1.28e-4, 9.17e-3),
    (2.49e-4, 0.0184),
    (3.54e-4, 0.0274),
    (4.50e-4, 0.0364),
    (5.38e-4, 0.0456),
    (6.20e-4, 0.0548),
    (6.95e-4, 0.0641),
    (7.65e-4, 0.0735),
    (8.30e-4, 0.0829),
    (8.91e-4, 0.0924),
    (9.49e-4, 0.102),
    (1.00e-3, 0.111),
    (1.06e-3, 0.121),
    (1.11e-3, 0.131),
    (1.16e-3, 0.14),
    (1.20e-3, 0.15),
    (1.25e-3, 0.16),
    (1.30e-3, 0.169),
    (1.34e-3, 0.179),
    (1.39e-3, 0.188),
    (1.43e-3, 0.198),
    (1.48e-3, 0.208),
    (1.52e-3, 0.217),
    (1.57e-3, 0.227),
    (1.62e-3, 0.237),
    (1.67e-3, 0.246),
    (1.72e-3, 0.256),
    (1.77e-3, 0.265),
    (1.83e-3, 0.275),
    (1.89e-3, 0.284),
    (1.96e-3, 0.294),
    (2.03e-3, 0.303),
    (2.10e-3, 0.312),
    (2.18e-3, 0.322),
    (2.27e-3, 0.331),
    (2.36e-3, 0.34),
    (2.47e-3, 0.349),
    (2.58e-3, 0.358),
    (2.71e-3, 0.366),
    (2.86e-3, 0.375),
    (3.02e-3, 0.383),
    (3.19e-3, 0.391),
    (3.39e-3, 0.399),
    (3.62e-3, 0.406),
    (3.86e-3, 0.414),
    (4.13e-3, 0.42),
    (4.44e-3, 0.427),
    (4.81e-3, 0.434),
    (5.23e-3, 0.441),
    (5.69e-3, 0.447),
    (6.20e-3, 0.452),
    (6.75e-3, 0.458),
    (7.55e-3, 0.465),
    (8.13e-3, 0.469),
    (8.75e-3, 0.473),
    (9.40e-3, 0.477),
    (0.01, 0.48),
    (0.0113, 0.486),
    (0.0125, 0.491),
    (0.015, 0.498),
    (0.0175, 0.505),
    (0.0202, 0.511),
    (0.025, 0.518),
    (0.0301, 0.525),
    (0.0351, 0.531),
    (0.0401, 0.536),
    (0.0452, 0.541),
    (0.0501, 0.545),
    (0.0573, 0.551),
    (0.0695, 0.56),
    (0.0793, 0.567),
    (0.089, 0.573),
    (0.0987, 0.58),
    (0.117, 0.591),
    (0.137, 0.602),
    (0.156, 0.616),
    (0.176, 0.633),
    (0.195, 0.649),
    (0.215, 0.671),
    (0.235, 0.694),
    (0.256, 0.719),
    (0.276, 0.746),
    (0.297, 0.778),
    (0.317, 0.817),
    (0.337, 0.861),
    (0.36, 0.914),
    (0.381, 0.959),
    (0.4, 1.0),
    (0.42, 1.04),
    (0.44, 1.07),
    (0.46, 1.1),
    (0.48, 1.13),
    (0.5, 1.15),
    (0.52, 1.18),
    (0.54, 1.21),
    (0.56, 1.23),
    (0.58, 1.26),
    (0.6, 1.3),
    (0.62, 1.34),
    (0.64, 1.38),
    (0.66, 1.44),
    (0.68, 1.5),
    (0.7, 1.57),
    (0.72, 1.65),
    (0.74, 1.72),
    (0.76, 1.78),
    (0.78, 1.84),
    (0.8, 1.91),
    (0.819, 1.98),
    (0.84, 2.07),
    (0.86, 2.17),
    (0.879, 2.3),
    (0.899, 2.45),
    (0.919, 2.63),
    (0.94, 2.9),
    (0.954, 3.27),
    (0.964, 3.63),
    (0.975, 4.22),
    (0.984, 5.09),
    (0.988, 5.9),
])[:63]

# convert first column to pressure
p0 = 0.981313*bar
exp_pressure = ref_data[:,0]*p0
# convert second column to occupation, i.e. the number of argon molecules per
# site. A site is defined as a part of the area that has the size of the Argon
# cross-section.
# STD stands for standard temperature and pressure. IUPAC says:
#    standard pressure = 1 bar
#    standard temperature = 273.15 K
alpha_ref = 3.108*centimeter**3/gram # from paper
stp_dens_ar = 44.072*mol/meter**3    # density of argon under STP conditions
spec_area = 6.01*meter**2/gram       # specific area of the graphite sample in the paper
cross_sec_ar = 0.138*nanometer**2    # the cross section of one argon site on the surface
occ_reference = alpha_ref*stp_dens_ar/spec_area*cross_sec_ar
exp_occupation = ref_data[:,1]*occ_reference

# compute the ideal gas estimate of the occupation using the thermodynamic model
mod_occupation = []
for p in exp_pressure:
    # compute the density of argon in the gas phase (3D)
    dens_ar = p/(boltzmann*temp)
    # use the equilibrium constant of the thermodynamic model to compute the
    # density on the 2D surface.
    K = tm.equilibrium_constant(temp)
    dens2d = K*dens_ar
    # from density to occupation using the cross section
    occupation = dens2d*cross_sec_ar
    mod_occupation.append(occupation)
mod_occupation = np.array(mod_occupation)

# Make a plot of the adsorption isotherm with the ideal gas and with the
# Langmuir model. Remaining deviations with experiment are due to (relatively
# strong) Ar-Ar interactions on the graphite surface, i.e. the 2D gas is far
# from ideal. One must do Monte-Carlo simulations to estimate these effects.
import matplotlib.pyplot as pt
pt.clf()
pt.plot(exp_pressure/bar, exp_occupation, 'k-', lw=2, label='experiment')
pt.plot(exp_pressure/bar, mod_occupation, 'b-', lw=1, label='TAMkin ideal gas')
pt.plot(exp_pressure/bar, mod_occupation/(1+mod_occupation), 'g-', lw=1, label='TAMkin Langmuir')
pt.ylim(0,1.0)
pt.xlim(0,5e-3)
pt.ylabel("Fraction of occupied mono-layer\nArgon sites on graphite at T=87K")
pt.xlabel("Pressure [bar]")
pt.legend(loc=0)
pt.savefig("adsorption_%s_isotherm.png" % level)
