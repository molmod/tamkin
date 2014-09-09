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


from tamkin.nma import NMA

from molmod import angstrom, lightspeed, centimeter
from molmod.io import XYZWriter
from molmod.periodic import periodic

import numpy as np


__all__ = ["dump_modes_xyz", "dump_modes_molden", "dump_modes_gaussian"]


def dump_modes_xyz(nma, indexes=0, prefix="mode", amplitude=5.0*angstrom, frames=36):
    """Write XYZ trajectory file(s) that vizualize internal mode(s)

       Arguments:
         | nma  --  an object that specifies the normal modes, several formats
                    are supported: (i) a Tamkin NMA object, (ii) a 3-tuple with
                    reference coordinates, mass-unweighted mode(s) and atom
                    numbers or (iii) a 4-tuple with reference coordinates, mass-
                    weighted mode(s), atom numbers and a masses3 vector. the
                    latter is a vector with 3*N elements containing the masses of
                    the atoms in groups of three.

       Optional arguments:
         | indexes  --  the index or a list of indexes of modes that must be
                        written to trajectory files [default=0]
         | prefix  --  a prefix used for the output files. the generated
                       trajectory filenames have the format prefix.index.xyz
                       [default="mode"]
         | amplitude  --  the amplitude of the normal mode vibration in atomic
                          untis [default=5*angstrom]
         | frames  --  the number of frames written to the trajectory file
                       [default=36]
    """

    if isinstance(nma, NMA):
        coordinates = nma.coordinates
        modes = nma.modes
        numbers = nma.numbers
        masses3 = nma.masses3
    elif hasattr(nma, "__len__") and len(nma)==3:
        coordinates, modes, numbers = nma
        masses3 = None
    elif hasattr(nma, "__len__") and len(nma)==4:
        coordinates, modes, numbers, masses3 = nma
    else:
        raise TypeError("Could not understand first argument. Check documentation.")

    if not hasattr(indexes, "__len__"):
        indexes = [indexes]

    if len(modes.shape) == 1:
        modes = modes.reshape((-1,1))

    symbols = [periodic[n].symbol for n in numbers]

    for index in indexes:
        filename = "%s.%i.xyz" % (prefix, index)
        mode = modes[:,index]
        if masses3 is not None:
            mode /= np.sqrt(masses3)
        mode /= np.linalg.norm(mode)
        xyz_writer = XYZWriter(filename, symbols)
        for frame in xrange(frames):
            factor = amplitude*np.sin(2*np.pi*float(frame)/frames)
            xyz_writer.dump("frame %i" % frame, coordinates + factor*mode.reshape((-1,3)))
        del xyz_writer


def dump_modes_gaussian(filename, nma, selected=None):
    """Write freqs and modes to a file in the Gaussian log format.

       Arguments:
         | filename  --  modes are written to this file,
                         can be read by Molden (visualization program)
         | nma  --  an NMA object or a tuple or list with five elements: modes,
                    frequencies, masses, numbers, coordinates. See
                    _make_moldenfile for details.

       Optional argument:
         | selected  --  Selection of modes for which to make
                         trajectories. This can be a list or array
                         of mode indices (length <= N), or an
                         array of booleans (length = N).

       The output file will look like a stripped Gaussian03 or Gaussian09 log
       file. It is sufficient to visualize the modes in Molden, Openbabel or
       Avogadro.
    """
    ### A) Parse the NMA argument
    if isinstance(nma, NMA):
        # NMA object
        modes, freqs, masses, numbers, coordinates = \
            nma.modes, nma.freqs, nma.masses, nma.numbers, nma.coordinates
    elif hasattr(nma, "__len__") and len(nma) == 5 and not isinstance(nma, np.ndarray):
        # [modes, freqs, ...] or (modes, freqs, ...)
        modes, freqs, masses, numbers, coordinates = nma
    else:
        raise TypeError("nma argument has wrong type")

    ### B) Select some modes
    if selected is not None:
        modes = np.take(modes, selected, 1)  # modes in columns
        freqs = freqs[selected]
        masses = masses[selected]
        numbers = numbers[selected]
        coordinates = coordinates[selected]

    ### C) convert modes to the right convention
    masses3_sqrt1 = np.array(sum([[1/m,1/m,1/m] for m in np.sqrt(masses)],[]))
    nmode = modes.shape[1]
    modes = modes.copy() # avoid modifying the given modes
    for imode in xrange(nmode):
        modes[:,imode] *= masses3_sqrt1
        modes[:,imode] /= np.linalg.norm(modes[:,imode])

    ### D) Define some multiline text blobs that are used below.
    header = """\
 Entering Gaussian System, Link 0=g09

 This file is generated from the dump_modes_gaussian function in the file
 trajectory.py of TAMkin. Please note, that this is a "fake" output; TAMkin
 doesn't compute intensities. This file should be readable by Molden, Openbabel
 and Avogadro.

 Gaussian, Inc

 # Fake method line for openbabel
 """

    header_coordinates = """\
                         Standard orientation:
 ---------------------------------------------------------------------
 Center     Atomic      Atomic             Coordinates (Angstroms)
 Number     Number       Type             X           Y           Z
 ---------------------------------------------------------------------"""

    header_basisfunctions = """\
 ---------------------------------------------------------------------
       basis functions          primitive gaussians
       alpha electrons          beta electrons
 **********************************************************************"""

    header_freq = """\
 Harmonic frequencies (cm**-1), IR intensities (KM/Mole),
 Raman scattering activities (A**4/AMU), Raman depolarization ratios,
 reduced masses (AMU), force constants (mDyne/A) and normal coordinates:"""

    txt_freqbelow_1 = """\
 Red. masses --     0.0000
 Frc consts  --     0.0000
 IR Inten    --     0.0000
  Atom  AN      X      Y      Z"""
    txt_freqbelow_2 = """\
 Red. masses --     0.0000                 0.0000
 Frc consts  --     0.0000                 0.0000
 IR Inten    --     0.0000                 0.0000
  Atom  AN      X      Y      Z        X      Y      Z"""
    txt_freqbelow_3 = """\
 Red. masses --     0.0000                 0.0000                 0.0000
 Frc consts  --     0.0000                 0.0000                 0.0000
 IR Inten    --     0.0000                 0.0000                 0.0000
  Atom  AN      X      Y      Z        X      Y      Z        X      Y      Z"""
    txt_freqbelow = [txt_freqbelow_1, txt_freqbelow_2, txt_freqbelow_3]

    ### E) Actual writing of the fake log file
    with open(filename, "w") as f:
        print >> f, header

        ### E1) ATOM PART
        assert modes.shape[0] % 3 == 0
        natom = modes.shape[0]/3
        print >> f, header_coordinates
        for iatom in range(natom):
           print >> f, '%5d %10d %13s %15f %11f %11f' %(
                       iatom + 1, numbers[iatom], "0",
                       coordinates[iatom,0]/angstrom,
                       coordinates[iatom,1]/angstrom,
                       coordinates[iatom,2]/angstrom)

        ### E2) ORBITAL PART
        print >> f, header_basisfunctions
        print >> f, " "   #this part is just empty

        ### D3) FREQUENCY PART
        # Modes are written in sections with three columns each.
        print >> f, header_freq
        # istart: the current mode index, will be incremented by 3 in each iteration.
        istart = 0
        while istart < nmode:
            # select modes and freqs for this block
            iend = min(istart+3, nmode)
            ncol = iend - istart # number of columns in section
            # print stuff to file
            #  - mode indexes
            for imode in xrange(istart, iend):
                print >> f, '%22d' % (imode + 1),
            print >> f
            #  - (fake) symmetry info
            print >> f, ' '.join(["                    ?A"]*ncol)
            #  - frequencies converted to inverse centimeters
            print >> f, ' Frequencies --',
            for imode in xrange(istart, iend):
                print >> f, '%10.4f' % (freqs[imode]/lightspeed*centimeter),
                if imode != iend-1:
                    print >> f, '           ',
            print >> f
            #  - a blob of text below the frequencies
            print >> f, txt_freqbelow[ncol-1]
            #  - the modes
            for iatom in range(natom):
                print >> f, '%6d %3d' % (iatom + 1, numbers[iatom]),
                for imode in xrange(istart, iend):
                    print >> f, '%8.2f %6.2f %6.2f' % (
                       modes[3*iatom, imode], modes[3*iatom+1, imode],
                       modes[3*iatom+2, imode]),
                print >> f
            # prepare for next iteration
            istart = iend

        print >> f, " Normal termination of Gaussian 09."

# For backward compatibility with TAMkin before version 1.0.2
dump_modes_molden = dump_modes_gaussian
