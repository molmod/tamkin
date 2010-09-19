# -*- coding: utf-8 -*-
# TAMkin is a post-processing toolkit for thermochemistry and kinetics analysis.
# Copyright (C) 2008-2010 Toon Verstraelen <Toon.Verstraelen@UGent.be>,
# Matthias Vandichel <Matthias.Vandichel@UGent.be> and
# An Ghysels <An.Ghysels@UGent.be>, Center for Molecular Modeling (CMM), Ghent
# University, Ghent, Belgium; all rights reserved unless otherwise stated.
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
# parts of this program are required to cite the following five articles:
#
# "TAMkin: A Versatile Package for Vibrational Analysis and Chemical Kinetics",
# An Ghysels, Toon Verstraelen, Karen Hemelsoet, Michel Waroquier and Veronique
# Van Speybroeck, Journal of Chemical Information and Modeling, Articles ASAP
# (As Soon As Publishable)
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
# --


from tamkin.data import Molecule

from molmod.molecules import Molecule as BaseMolecule
from molmod.periodic import periodic
from molmod.unit_cells import UnitCell
from molmod.units import angstrom

import numpy


__all__ = ["load_molecule_cp2k"]


def load_molecule_cp2k(fn_xyz, fn_sp, fn_freq, multiplicity=1, is_periodic=True):
    """Load a molecule with the Hessian from a CP2K computation

       Arguments:
        | fn_xyz  --  The filename of the xyz file containing the (partially)
                      optimized geometry.
        | fn_sp   --  The filename of the single point .out file containing the
                      energy.
        | fn_freq  --  The filename of the frequency .out file containing the
                       hessian

       Optional arguments:
        | multiplicity  --  The spin multiplicity of the electronic system
                            [default=1]
        | is_periodic  --  True when the system is periodic in three dimensions.
                           False when the systen is aperiodic. [default=True]
        | unit_cell  --  The unit cell vectors for periodic structures
    """
    molecule = BaseMolecule.from_file(fn_xyz)
    masses = numpy.array([periodic[number].mass for number in molecule.numbers])

    # go through the single point file: energy and gradient
    energy = None
    gradient = None
    f = file(fn_sp)
    while True:
        line = f.readline()
        if line == "":
            break
        if line.startswith(" ENERGY|"):
            energy = float(line[60:])
        elif line.startswith(" FORCES|"):
            tmp = []
            while True:
                line = f.readline()
                if line == "\n":
                    break
                if line == "":
                    raise IOError("End of file while reading gradient (forces).")
                words = line.split()
                tmp.append([float(words[1]), float(words[2]), float(words[3])])
            gradient = -numpy.array(tmp) # force to gradient
            break
    if energy is None or gradient is None:
        raise IOError("Could not read energy and/or gradient (forces) from single point file.")
    f.close()

    # go through the freq file: lattic vectors and hessian
    f = file(fn_freq)
    vectors = numpy.zeros((3,3),float)
    while True:
        line = f.readline()
        if line.startswith(" CELL"): break
    for axis in range(3):
        line = f.readline()
        vectors[:,axis] = numpy.array( [float(line[29:39]), float(line[39:49]), float(line[49:59])] )
    unit_cell = UnitCell(vectors*angstrom)

    hessian = None
    while True:
        line = f.readline()
        if line.startswith(" VIB| Hessian in cartesian coordinates"):
            block_len = molecule.size*3
            tmp = numpy.zeros((block_len,block_len), float)
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
    conv = 1e-3*numpy.array([masses, masses, masses]).transpose().ravel()**0.5
    hessian *= conv
    hessian *= conv.reshape((-1,1))

    return Molecule(
        molecule.numbers, molecule.coordinates, masses, energy, gradient,
        hessian, multiplicity, 0, is_periodic, unit_cell=unit_cell
    )
