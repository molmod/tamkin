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


def load_molecule_vasp(vaspfile_xyz, vaspfile_out, energy = 0.0, multiplicity=1, is_periodic=True):
    """Load a molecule from VASP output files

       Arguments:
        | vaspfile_xyz  --  Filename of xyz-file containing the (partially)
                            optimized structure.
        | vaspfile_out  --  Filename of VASP output file (OUTCAR) containing
                            eg the Hessian.

       Optional arguments:
        | energy  --  The electronic energy. [default=0.0].
        | multiplicity  --  The spin multiplicity of the electronic system
                            [default=1]
        | is_periodic  --  True when the system is periodic in three dimensions.
                           False when the systen is nonperiodic. [default=True].
    """
    # TODO: read energy from VASP file?
    # Units: VASP gradient in eV/angstrom, TAMkin internally all in atomic units

    # Read atomtypes from xyz-VASP-file
    # format:  one atom per line
    #          atomtype (Si, C, ...)   x-cor  y-cor  z-cor
    f = open(vaspfile_xyz)
    atomtypes = []
    for line in f:
        words = line.split()
        if len(words) == 4:
            atomtypes.append(words[0])
    f.close()
    atomtypes = np.array(atomtypes)

    # Read other data from out-VASP-file OUTCAR
    f = open(vaspfile_out)

    # number of atoms (N)
    for line in f:
        if line.strip().startswith("Dimension of arrays:"):  break
    f.next()
    for line in f:
        words = line.split()
        N = int(words[-1])
        break

    # read lattice vectors: store in columns
    vectors = np.ones((3,3),float)
    for line in f:
        if line.startswith("      direct lattice vectors"): break
    axis = 0
    for line in f:
        words = line.split()
        vectors[:,axis] = np.array([ float(word)*angstrom for word in words[:3]])
        axis += 1
        if axis >= 3: break
    unit_cell = UnitCell(vectors)

    # masses
    # TODO: should be made more general?
    masses = np.zeros((N),float)
    for at,atomtype in enumerate(atomtypes):
        table = { "H": 1.000,   "C": 12.011, "O": 16.000,
                  "Al": 26.982, "Si": 28.085,
                }
        masses[at] = table[atomtype]*amu

    # get corresponding atomic numbers
    atomicnumbers = np.zeros(N, int)
    mass_table = np.zeros(len(periodic))
    for i in xrange(1, len(mass_table)):
        m1 = periodic[i].mass
        if m1 is None:
            m1 = 200000.0
        m2 = periodic[i+1].mass
        if m2 is None:
            m2 = 200000.0
        mass_table[i] = 0.5*(m1+m2)
    for i,mass in enumerate(masses):
        atomicnumbers[i] = mass_table.searchsorted(mass)

    # positions, gradient
    positions = np.zeros((N,3),float)
    gradient = np.zeros((N,3),float)
    for line in f:          # go to first time POSITION is printed (reference point)
        if line.strip().startswith("POSITION"): break
    f.next()
    row = 0
    for line in f:
        words = line.split()
        positions[row,:] = [ float(word)*angstrom for word in words[:3] ]
        gradient[row,:]  = [-float(word)*electronvolt/angstrom for word in words[3:6] ]
        row += 1
        if row >= N: break

    # hessian, not symmetrized, useful to find indices of Hessian elements
    hessian = np.zeros((3*N,3*N),float)
    for line in f:
        if line.strip().startswith("SECOND DERIVATIVES (NOT SYMMETRIZED)"):  break
    f.next()
    for line in f:
        Nfree = len(line.split())/3   # nb of non-fixed atoms
        break
    # find the (cartesian) indices of non-fixed atoms
    indices_free = []
    row = 0
    mu = 0
    for line in f:
        if mu==0:
            atom = int(line.split()[0][:-1])
        indices_free.append(3*(atom-1) + mu)
        mu+=1
        row+=1
        if mu >= 3: mu=0
        if row >= 3*Nfree: break

    # hessian, symmetrized, somehow with a negative sign
    hessian = np.zeros((3*N,3*N),float)
    for line in f:
        if line.strip().startswith("SECOND DERIVATIVES (SYMMETRYZED)"):  break
    f.next()
    f.next()
    row = 0
    for line in f:   # put Hessian elements at appropriate place
        index_row = indices_free[row]
        words = line.split()
        for col in xrange(3*Nfree):
            index_col = indices_free[col]
            hessian[index_row,index_col] = -float(words[col+1])*electronvolt/angstrom**2 # skip first col
        row += 1
        if row >= 3*Nfree: break
    f.close()

    return Molecule(
        atomicnumbers, positions, masses, energy, gradient,
        hessian, multiplicity, periodic=is_periodic, unit_cell=unit_cell )


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
