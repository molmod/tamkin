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
"""Analysis of molecular geometries"""

import numpy as np


__all__ = ["transrot_basis", "rank_linearity"]


def transrot_basis(coordinates, rot=True):
    """Constructs 6 vectors which represent global translations/rotations.

       Argument:
        | coordinates  --  The atom coordinates (float numpy array with shape
                           Nx3)

       Optional argument:
        | rot  --  When True the rotations are included [default=True]

       The return value is a numpy array with 3*N columns and 6 (``rot==True``)
       or 3 (``rot==False``) rows. The rows are not mass weighted.
    """
    if not rot:
        result = np.zeros((3, coordinates.size), float)
    else:
        result = np.zeros((6, coordinates.size), float)
    # translation
    result[0, 0::3] = 1
    result[1, 1::3] = 1
    result[2, 2::3] = 1
    if rot:
        # rotation
        result[3, 1::3] =  coordinates[:,2]
        result[3, 2::3] = -coordinates[:,1]
        result[4, 2::3] =  coordinates[:,0]
        result[4, 0::3] = -coordinates[:,2]
        result[5, 0::3] =  coordinates[:,1]
        result[5, 1::3] = -coordinates[:,0]
    return result


def rank_linearity(coordinates, svd_threshold=1e-5, masses3=None):
    """Test the linearity of the given coordinates

       Arguments:
        | coordinates  --  The atom coordinates (float numpy array with shape
                           Nx3)

       Optional argument:
        | svd_threshold  --  Defines the sensitivity for deviations from
                             linearity [default=1e-5]
        | masses3  --  array with atomic masses (each mass is repeated trice)

       Returns the number of external degrees of freedom (dof) and the
       corresponding basis. The dof can be

       * 6 dof if atoms of system are non-collinear
       * 5 dof if atoms of system are collinear, e.g. when the subsystem
         contains only 2 atoms
       * 3 dof if system contains just 1 atom

       If masses3 is given, the analysis is performed in mass-weighted
       coordinates and the returned basis is in mass-weighted coordinates.
       If masses3 is not given everything is done in non-mass-weighted Cartesian
       coordinates.
    """
    if masses3 is None:
        center = coordinates.mean(0)  # geometrical center
        transrot = transrot_basis(coordinates - center)
    else:
        center = (coordinates*masses3.reshape((-1,3))).sum(0)/(masses3.sum()/3)  # center of mass
        #print ((coordinates-center)*masses3.reshape((-1,3))).sum(0)  # check
        transrot = transrot_basis(coordinates - center)
        transrot *= np.sqrt(masses3)

    U, W, Vt = np.linalg.svd(transrot, full_matrices=False)
    rank = (abs(W) > abs(W[0])*svd_threshold).sum()
    external_basis = Vt[:rank]

    return rank, external_basis
