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


# Import the tamkin library.
from tamkin import *
from molmod import * # for the units
from molmod.isotopes import ame2003


def create_kinetic_model(mol_react, mol_trans):
    # Perform normal mode analysis on the three molecules
    nma_react = NMA(mol_react, ConstrainExt())
    nma_trans = NMA(mol_trans, ConstrainExt())
    # Construct the two partition functions.
    pf_react = PartFun(nma_react, [ExtTrans(), ExtRot(1), Vibrations(classical=True, freq_scaling=0.9085)])
    pf_trans = PartFun(nma_trans, [ExtTrans(), ExtRot(1), Vibrations(classical=True, freq_scaling=0.9085)])

    return KineticModel([pf_react], pf_trans)


old_mol_react = load_molecule_g03fchk("reactant.fchk")
old_mol_trans = load_molecule_g03fchk("trans.fchk")
old_km = create_kinetic_model(old_mol_react, old_mol_trans)
print "Original rate constant at 303K =", old_km.rate_constant(303)/old_km.unit, old_km.unit_name

# does not work:
# mol_react1.masses[0] = 13*amu

new_masses = old_mol_react.masses.copy()
new_masses[20] = ame2003.masses[7][15]
new_mol_react = old_mol_react.copy_with(masses=new_masses)
new_mol_trans = old_mol_trans.copy_with(masses=new_masses)
new_km = create_kinetic_model(new_mol_react, new_mol_trans)
print "New rate constant at 303K =", new_km.rate_constant(303)/new_km.unit, new_km.unit_name

print "Ratio at 303K =", old_km.rate_constant(303)/new_km.rate_constant(303)
