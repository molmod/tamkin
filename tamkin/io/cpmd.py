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


from tamkin.data import Molecule

from molmod import angstrom, amu
from molmod.periodic import periodic

import numpy as np


__all__ = ["load_molecule_cpmd"]


def load_molecule_cpmd(fn_out, fn_geometry, fn_hessian, multiplicity=1, is_periodic=True):
    """Load a molecule with the Hessian from a CPMD computation

       Arguments:
        | fn_out  --  The filename of the output containing the total energy.
        | fn_geometry   --  The filename of geometry and the gradient of the
                            (partially) optimized system. (This filename is
                            typically GEOMETRY.xyz.)
        | fn_hessian  --  The filename of the of the file containing the
                          Hessian. (This filename is typically MOLVIB.)

       Optional arguments:
        | multiplicity  --  The spin multiplicity of the electronic system
                            [default=1]
        | is_periodic  --  True when the system is periodic in three dimensions.
                           False when the system is aperiodic. [default=True]
    """
    # go through the output file: grep the total energy
    energy = None
    f = file(fn_out)
    while True:
        line = f.readline()
        if line == "":
            raise IOError("Could not find final results in %s. Is the output file truncated?" % fn_out)
        if line == " *                        FINAL RESULTS                         *\n":
            break
    while True:
        line = f.readline()
        if line == "":
            raise IOError("Could not find total energy in %s. Is the output file truncated?" % fn_out)
        if line.startswith(" (K+E1+L+N+X)           TOTAL ENERGY ="):
            words= (line.strip()).split()
            energy = float(words[4])
            break
    f.close()

    # load the optimal geometry
    f = file(fn_geometry)
    num_atoms = int(f.readline())
    numbers = np.zeros(num_atoms, int)
    coordinates = np.zeros((num_atoms,3), float)
    gradient = np.zeros((num_atoms,3), float)

    f.readline()
    i = 0
    while True:
        line = f.readline()
        if line == "":
            break # end of file
        words = (line.strip()).split()
        if len(words) == 7:
            numbers[i] = periodic[words[0]].number
            coordinates[i][0] = float(words[1])*angstrom
            coordinates[i][1] = float(words[2])*angstrom
            coordinates[i][2] = float(words[3])*angstrom
            gradient[i][1] = float(words[4])
            gradient[i][1] = float(words[5])
            gradient[i][2] = float(words[6])
            i += 1
        else:
            raise IOError("Expecting seven words at each atom line in %s." % fn_geometry)
    if i != num_atoms:
        raise IOError("The number of atoms is incorrect in %s." % fn_geometry)
    f.close()

    # go trhough the freq file: hessian
    f = file(fn_hessian)

    line = f.readline()
    if not line.startswith(" &CART"):
        raise IOError("File %s does not start with &CART." % fn_hessian)
    masses = np.zeros(num_atoms, float)
    for i in xrange(num_atoms):
        line = f.readline()
        words = line.split()
        masses[i] = float(words[4])*amu
    f.readline() # &END

    line = f.readline()
    if not line.startswith(" &FCON"):
        raise IOError("File %s does not contain section &FCON." % fn_hessian)
    num_cart = num_atoms*3
    hessian = np.zeros((num_cart, num_cart), float)
    for i in xrange(num_cart):
        line = f.readline()
        words = line.split()
        for j in xrange(num_cart):
            hessian[i,j] = float(words[j])

    f.close()

    return Molecule(
        numbers, coordinates, masses, energy, gradient, hessian, multiplicity,
        0, is_periodic
    )
