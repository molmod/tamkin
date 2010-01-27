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


from tamkin.nma import NMA

from molmod import angstrom
from molmod.io import XYZWriter
from molmod.periodic import periodic

import numpy


__all__ = ["dump_modes_xyz"]


def dump_modes_xyz(nma, indexes, prefix=None, amplitude=5.0*angstrom, frames=36):
    """Write XYZ trajectory file(s) that vizualize internal mode(s)

       Arguments:
         nma  --  an object that specifies the normal modes, several formats
                  are supported: (i) a Tamkin NMA object, (ii) a 3-tuple with
                  reference coordinates, mass-unweighted modes and atom numbers
                  or (iii) a 4-tuple with reference coordinates, mass-weighted
                  modes, atom numbers and a masses3 vector. the latter is a
                  vector with 3*N elements containing the masses of the
                  atoms.
         indexes  --  the index or a list of indexes of modes that must be
                      written to trajectory files

       Optional arguments:
         prefix  --  a prefix used for the output files. the generated
                     trajectory filenames have the format prefix.index.xyz
                     [default="mode"]
         amplitude  --  the amplitude of the normal mode vibration in atomic
                        untis [default=5*angstrom]
         frames  --  the number of frames written to the trajectory file
                     [default=36]
    """

    if isinstance(nma, NMA):
        coordinates = nma.coordinates
        modes = nma.modes
        numbers = nma.numbers
        masses3 = nma.masses3
    elif hasattr(nma, "__len__") and len(nma==3):
        coordinates, modes, numbers = nma
        masses3 = None
    elif hasattr(nma, "__len__") and len(nma==4):
        coordinates, modes, numbers, masses3 = nma
    else:
        raise TypeError("Could not understand first argument. Check documentation.")

    if not hasattr(indexes, "__len__"):
        indexes = [indexes]

    if prefix is None:
        prefix = "mode"

    symbols = [periodic[n].symbol for n in numbers]

    for index in indexes:
        filename = "%s.%i.xyz" % (prefix, index)
        mode = nma.modes[:,index]
        if masses3 is not None:
            mode /= numpy.sqrt(masses3)
        xyz_writer = XYZWriter(filename, symbols)
        for frame in xrange(frames):
            factor = amplitude*numpy.sin(2*numpy.pi*float(frame)/frames)
            xyz_writer.dump("frame %i" % frame, coordinates + factor*mode.reshape((-1,3)))
        del xyz_writer

