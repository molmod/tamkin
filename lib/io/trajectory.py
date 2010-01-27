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

from molmod import angstrom


__all__ = ["write_modes_for_VMD"]


def write_modes_for_VMD(nma, index, filename=None, A=50.0, frames=36):
    """This function selects calls the function write_modes_for_VMD_2,
    where the mode trajectory is actually written.
    The function selects the relevant attributes.
    Numbering modes starts at 0."""

    if filename is None: filename = "mode"+str(index)+".txt"

    # Select mode from the nma.modes and undo mass-weighting
    mode = nma.modes[:,index]
    for i in range(len(mode)):
        mode[i] /= numpy.sqrt(nma.masses3[i])

    write_modes_for_VMD_2(nma.numbers, nma.coordinates, mode, filename=filename, A=A, frames=frames)


def write_modes_for_VMD_2(atomicnumbers, coordinates, mode, filename=None, A=50.0, frames=36):

    import math
    if filename is None: filename = "mode.xyz"

    nbatoms = len(atomicnumbers)
    positions = coordinates/angstrom   # units for VMD

    f = open(filename,'w')

    for time in range(frames+1):
        factor = A * math.sin( 2*math.pi * float(time)/frames)
        print >> f, nbatoms
        print >> f, 'i='+str(time)
        for at in range(nbatoms):
            coor = positions[at,:] + factor * mode[3*at:3*at+3]
            print >> f, "%-5i  %10.4f  %10.4f  %10.4f" %(atomicnumbers[at], coor[0],coor[1],coor[2])
    f.close


