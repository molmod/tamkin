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
# --


from __future__ import print_function, division

from tamkin.data import Molecule

from molmod import amu

import numpy as np


__all__ = [
    "load_molecule_molpro",
]

def load_molecule_molpro(filename):
    """Load a molecule from Molpro 2012 output file.

       Argument:
         | ``filename`` -- the molpro 2012 output

    """

    lines = open(filename,"r").read().splitlines()
    # locate positions
    beginxyz = None
    begincoordinate = None
    beginmass = None
    beginhessian = None
    begingradient = None
    for lineno, line in enumerate(lines):
        if "Current geometry" in line:
            beginxyz = lineno + 2
        if "FREQUENCIES * CALCULATION OF NORMAL MODES" in line:
            begincoordinate = lineno + 7
        if "Atomic Masses" in line:
            beginmass = lineno
        if "Force Constants" in line:
            beginhessian = lineno
        # if "GRADIENT FOR" in line:
            # begingradient = lineno + 4
    # xyz read, only meta data and energy
    atomnumber = int(lines[beginxyz].split()[0])
    title = lines[beginxyz+1].split()[0]
    energy = float(lines[beginxyz+1].split("=")[1])
    # print(atomnumber, title, energy)

    # coordinate read
    numbers=[]
    coordinates=[]

    for line in lines[begincoordinate : begincoordinate+atomnumber]:
        words = line.split()
        charge = int(float(words[2]))
        x = float(words[3])
        y = float(words[4])
        z = float(words[5])
        numbers.append(charge)
        coordinates.append([x, y, z])
        # print(charge, x, y, z)

    # masses
    masses = []
    for line in lines[beginmass+1:]:
        if len(line) > 8:
            words = line[8:].split()
            for word in words:
                masses.append(float(word)*amu)
        else:
            break
    # print(masses)

    # gradient
    gradient = np.zeros((atomnumber, 3))
    # TODO, Gradient is more difficult than I thought...
    # format of gradient from CCSD(T)-F12 is different from DFT
    '''
    if begingradient != None:
        for i, line in enumerate(lines[begingradient: begingradient + atomnumber]):
            words = line.split()
            gradient[i,0] = float(words[1])
            gradient[i,1] = float(words[2])
            gradient[i,2] = float(words[3])
    '''
    # print(gradient)

    # hessian
    hessian = np.ndarray((3*atomnumber, 3*atomnumber), dtype= float)
    headerline = beginhessian+1
    rowbegin = 0
    while rowbegin < 3*atomnumber:
        for r, line in zip(\
            range(rowbegin, 3*atomnumber),\
            lines[headerline+1 : headerline+1+3*atomnumber-rowbegin]\
            ):
            words = line[16:].split()
            # print(words)
            for c in range(rowbegin, rowbegin+min(5, 3*atomnumber-rowbegin, len(words))):
                # print(r,c)
                hessian[r,c] = float(words[c-rowbegin])
                hessian[c,r]=hessian[r,c]
        headerline += 3*atomnumber-rowbegin+1
        rowbegin += 5

    # print(hessian)

    return Molecule(
        np.array(numbers),
        np.array(coordinates),
        np.array(masses),
        energy,
        gradient,
        hessian,
        title=title,
        multiplicity=1,
    )
