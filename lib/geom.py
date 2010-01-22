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
# "Vibrational Modes in partially optimized molecular systems.", An Ghysels,
# Dimitri Van Neck, Veronique Van Speybroeck, Toon Verstraelen and Michel
# Waroquier, Journal of Chemical Physics, Vol. 126 (22): Art. No. 224102, 2007
# DOI:10.1063/1.2737444
#
# "Cartesian formulation of the Mobile Block Hesian Approach to vibrational
# analysis in partially optimized systems", An Ghysels, Dimitri Van Neck and
# Michel Waroquier, Journal of Chemical Physics, Vol. 127 (16), Art. No. 164108,
# 2007
# DOI:10.1063/1.2789429
#
# "Calculating reaction rates with partial Hessians: validation of the MBH
# approach", An Ghysels, Veronique Van Speybroeck, Toon Verstraelen, Dimitri Van
# Neck and Michel Waroquier, Journal of Chemical Theory and Computation, Vol. 4
# (4), 614-625, 2008
# DOI:10.1021/ct7002836
#
# "Mobile Block Hessian approach with linked blocks: an efficient approach for
# the calculation of frequencies in macromolecules", An Ghysels, Veronique Van
# Speybroeck, Ewald Pauwels, Dimitri Van Neck, Bernard R. Brooks and Michel
# Waroquier, Journal of Chemical Theory and Computation, Vol. 5 (5), 1203-1215,
# 2009
# DOI:10.1021/ct800489r
#
# "Normal modes for large molecules with arbitrary link constraints in the
# mobile block Hessian approach", An Ghysels, Dimitri Van Neck, Bernard R.
# Brooks, Veronique Van Speybroeck and Michel Waroquier, Journal of Chemical
# Physics, Vol. 130 (18), Art. No. 084107, 2009
# DOI:10.1063/1.3071261
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


import numpy


__all__ = ["transrot_basis","rank_linearity"]


def transrot_basis(coordinates, rot=True):
    """Constructs 6 vectors which represent global translations/rotations. (dim: 6xsize)
    Vectors are not mass-weighted.
    If rot=False, only translational vectors.
    """
    if not rot:
        result = numpy.zeros((3, coordinates.size), float)
    else:
        result = numpy.zeros((6, coordinates.size), float)
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


def rank_linearity(coordinates,svd_threshold=1e-5):
    """
    Linearity of system with given coordinates (degrees of freedom = dof)
    6 dof if atoms of system are non-collinear
    5 dof if atoms of system are collinear, e.g. when the subsystem contains only 2 atoms
    3 dof if system contains just 1 atom

    method:
    Construct a kind of inertia matrix (6x6) of the system
              A = transrot_basis . transrot_basis**T
    Diagonalize A. The rank is the number of expected zero freqs.

    transrot_basis contains the 6 global translations and rotations of the subsystem (6 x 3*molecule.size)
    """
    # TODO clean up printing statements
    #print "coooor", coordinates
    corr = numpy.sum(coordinates,0)/len(coordinates)  # geometrical center
    #print "corr", corr
    #print numpy.resize(corr,(len(coordinates),3))
    #print "coooor", coordinates - numpy.resize(corr,(len(coordinates),3))
    transrot = transrot_basis(coordinates - numpy.resize(corr,(len(coordinates),3)) )
    #print "transrot", transrot

    A    = numpy.dot( transrot, transrot.transpose() )
    eigv = numpy.linalg.eigvalsh(A)
    #print "svd: eigv", eigv
    rank = (abs(eigv) > abs(eigv[-1])*svd_threshold).sum()
    #print "rank", rank
    return rank
