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

from molmod.periodic import periodic
from molmod.unit_cells import UnitCell
from molmod.units import angstrom, amu

import numpy as np


__all__ = ["load_molecule_cp2k"]


def load_molecule_cp2k(fn_sp, fn_freq, multiplicity=1, is_periodic=True):
    """Load a molecule with the Hessian from a CP2K computation

       Arguments:
        | fn_sp   --  The filename of the single point .out file containing the
                      energy and the forces.
        | fn_freq  --  The filename of the frequency .out file containing the
                       hessian

       Optional arguments:
        | multiplicity  --  The spin multiplicity of the electronic system
                            [default=1]
        | is_periodic  --  True when the system is periodic in three dimensions.
                           False when the systen is aperiodic. [default=True]
        | unit_cell  --  The unit cell vectors for periodic structures
    """
    # auxiliary routine to read atoms
    def atom_helper(f):
        # skip some lines
        for i in xrange(3):
            f.readline()
        # read the atom lines until an empty line is encountered
        numbers = []
        coordinates = []
        masses = []
        while True:
            line = f.readline()
            if len(line.strip()) == 0:
                break
            symbol = line[14:19].strip()[:2]
            atom = periodic[symbol]
            if atom is None:
                symbol = symbol[:1]
                atom = periodic[symbol]
            if atom is None:
                numbers.append(0)
            else:
                numbers.append(atom.number)
            coordinates.append([float(line[22:33]), float(line[34:45]), float(line[46:57])])
            masses.append(float(line[72:]))

        numbers = np.array(numbers)
        coordinates = np.array(coordinates)*angstrom
        masses = np.array(masses)*amu
        return numbers, coordinates, masses


    # auxiliary routine to read forces
    def force_helper(f, skip, offset):
        # skip some lines
        for i in xrange(skip):
            f.readline()
        # Read the actual forces
        tmp = []
        while True:
            line = f.readline()
            if line == "\n":
                break
            if line == "":
                raise IOError("End of file while reading gradient (forces).")
            words = line.split()
            try:
                tmp.append([float(words[offset]), float(words[offset+1]), float(words[offset+2])])
            except StandardError:
                break
        return -np.array(tmp) # force to gradient

    # go through the single point file: energy and gradient
    energy = None
    gradient = None
    f = file(fn_sp)
    while True:
        line = f.readline()
        if line == "":
            break
        if line.startswith(" ENERGY|"):
            energy = float(line[58:])
        elif line.startswith(" MODULE") and "ATOMIC COORDINATES" in line:
            numbers, coordinates, masses = atom_helper(f)
        elif line.startswith(" FORCES|"):
            gradient = force_helper(f, 0, 1)
            break
        elif line.startswith(' ATOMIC FORCES in [a.u.]'):
            gradient = force_helper(f, 2, 3)
            break
    if energy is None or gradient is None:
        raise IOError("Could not read energy and/or gradient (forces) from single point file.")
    f.close()

    # go through the freq file: lattic vectors and hessian
    f = file(fn_freq)
    vectors = np.zeros((3,3),float)
    while True:
        line = f.readline()
        if line.startswith(" CELL"): break
    for axis in range(3):
        line = f.readline()
        vectors[:,axis] = np.array( [float(line[29:39]), float(line[39:49]), float(line[49:59])] )
    unit_cell = UnitCell(vectors*angstrom)

    hessian = None
    while True:
        line = f.readline()
        if line.startswith(" VIB| Hessian in cartesian coordinates"):
            block_len = coordinates.size
            tmp = np.zeros((block_len,block_len), float)
            i2 = 0
            while i2 < block_len:
                num_cols = min(5, block_len-i2)
                f.readline() # skip two lines
                f.readline()
                for j in xrange(block_len):
                    line = f.readline()
                    if line == "":
                        raise IOError("End of file while reading hessian.")
                    words = line.split()
                    for i1 in xrange(num_cols):
                        tmp[i2+i1,j] = float(words[i1+2])
                i2 += num_cols
            hessian = tmp
            break
    f.close()
    if hessian is None:
        raise IOError("Could not read hessian from freq file.")

    # symmetrize
    hessian = 0.5*(hessian+hessian.transpose())
    # cp2k prints a transformed hessian, here we convert it back to the normal
    # hessian in atomic units.
    conv = 1e-3*np.array([masses, masses, masses]).transpose().ravel()**0.5
    hessian *= conv
    hessian *= conv.reshape((-1,1))

    return Molecule(
        numbers, coordinates, masses, energy, gradient,
        hessian, multiplicity, 0, is_periodic, unit_cell=unit_cell
    )
