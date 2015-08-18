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

from molmod import electronvolt, angstrom, amu
from molmod.periodic import periodic
from molmod.unit_cells import UnitCell

import numpy as np


__all__ = ["load_molecule_vasp", "load_fixed_vasp"]


def load_molecule_vasp(vaspfile_xyz, vaspfile_out, energy, multiplicity=1, is_periodic=True):
    """Load a molecule from VASP 4.6.X and 5.3.X output files

       Arguments:
        | vaspfile_xyz  --  Filename of xyz-file containing the (partially)
                            optimized structure. Generate this file with the
                            con2xyz tool from VASP
        | vaspfile_out  --  Filename of VASP output file (OUTCAR) containing
                            the Hessian.
        | energy  --  The potential energy. (It is not printed in the OUTCAR
                      file of a Hessian calculation.)

       Optional arguments:
        | multiplicity  --  The spin multiplicity of the electronic system
                            [default=1]
        | is_periodic  --  True when the system is periodic in three dimensions.
                           False when the systen is nonperiodic. [default=True].
    """
    # Units:
    #   * VASP Atomic coordinates in angstrom, TAMkin internally in atomic units
    #   * VASP Hessian in eV/angstrom**2, TAMkin internally all in atomic units

    # Read symbols of the atoms from xyz-VASP-file
    # format:  one atom per line
    #          Symbol (Si, C, ...)   x-cor  y-cor  z-cor
    symbols = []
    coordinates = []
    with open(vaspfile_xyz) as f:
        # read number of atoms
        natom = int(f.next())
        # skip title
        f.next()
        # Read from each line that contains four words.
        for line in f:
            words = line.split()
            if len(words) == 4:
                symbols.append(words[0])
                coordinates.append([
                    float(words[1])*angstrom,
                    float(words[2])*angstrom,
                    float(words[3])*angstrom,
                ])
    if len(symbols) != natom:
        raise IOError('Inconsistent number of atoms in XYZ file.')
    coordinates = np.array(coordinates)

    # Read other data from out-VASP-file OUTCAR
    with open(vaspfile_out) as f:
        # Read lattice vectors. VASP writes each vector as a column, TAMkin
        # stores each vector as a row.
        for line in f:
            if line.startswith("      direct lattice vectors"):
                break
        rvecs = np.ones((3,3),float)
        for axis in xrange(3):
            words = f.next().split()
            # Tranpose as we read.
            rvecs[:,axis] = np.array([float(word)*angstrom for word in words[:3]])
        unit_cell = UnitCell(rvecs)

        # hessian, not symmetrized, useful to find indices of Hessian elements
        hessian = np.zeros((3*natom, 3*natom), float)
        for line in f:
            if line.startswith(" SECOND DERIVATIVES (NOT SYMMETRIZED)"):  break
        # skip one line
        f.next()
        # number of free atoms (not fixed in space)
        keys = f.next().split()
        nfree_dof = len(keys)
        indices_free = [3*int(key[:-1])+{'X': 0, 'Y': 1, 'Z': 2}[key[-1]]-3 for key in keys]
        assert nfree_dof % 3 == 0
        for ifree0 in xrange(nfree_dof):
            line = f.next()
            irow = indices_free[ifree0]
            # skip first col
            words = line.split()[1:]
            assert len(words) == nfree_dof
            for ifree1 in xrange(nfree_dof):
                icol = indices_free[ifree1]
                hessian[irow, icol] = -float(words[ifree1])*electronvolt/angstrom**2
        hessian = 0.5*(hessian + hessian.T)

    # Masses
    masses = np.array([periodic[symbol].mass for symbol in symbols])

    # Get corresponding atomic numbers
    numbers = np.array([periodic[symbol].number for symbol in symbols])

    # No decent gradient information available
    gradient = np.zeros(coordinates.shape)

    return Molecule(
        numbers, coordinates, masses, energy, gradient, hessian,
        multiplicity=multiplicity, periodic=is_periodic, unit_cell=unit_cell)


def load_fixed_vasp(filename):
    """ Load list of fixed atoms from VASP output file

    Argument:
     | filename  --  Filename of VASP output file (OUTCAR)

    VASP can calculate partial Hessians: only a submatrix of the complete Hessian
    is computed with numerical differentiation. The rest of the Hessian elements
    is put to zero. This function determines which atoms have zero rows/cols in
    the Hessian, or, in other words, which were fixed.
    """
    # Read data from out-VASP-file OUTCAR
    f = open(filename)

    # number of atoms (N)
    for line in f:
        if line.strip().startswith("Dimension of arrays:"):  break
    f.next()
    for line in f:
        words = line.split()
        N = int(words[-1])
        break

    # hessian, not symmetrized, useful to find indices of Hessian elements
    for line in f:
        if line.strip().startswith("SECOND DERIVATIVES (NOT SYMMETRIZED)"):  break
    f.next()
    for line in f:
        Nfree = len(line.split())/3   # nb of non-fixed atoms
        break
    # find the non-fixed atoms
    atoms_free = []
    row = 0
    mu = 0
    for line in f:
        if mu==0:
            atom = int(line.split()[0][:-1])
            atoms_free.append(atom-1)
        mu+=1
        row+=1
        if mu >= 3: mu=0
        if row >= 3*Nfree: break
    f.close()

    fixed_atoms = [at for at in xrange(N) if at not in atoms_free]
    return np.array(fixed_atoms)
