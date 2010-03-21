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

from molmod import angstrom, lightspeed, centimeter
from molmod.io import XYZWriter
from molmod.periodic import periodic

import numpy


__all__ = ["dump_modes_xyz", "dump_modes_molden"]


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
            mode /= numpy.sqrt(masses3)
        mode /= numpy.linalg.norm(mode)
        xyz_writer = XYZWriter(filename, symbols)
        for frame in xrange(frames):
            factor = amplitude*numpy.sin(2*numpy.pi*float(frame)/frames)
            xyz_writer.dump("frame %i" % frame, coordinates + factor*mode.reshape((-1,3)))
        del xyz_writer



#======================================
#    write logfile as Gaussian03 does
#======================================

def dump_modes_molden(filename, nma, selected=None):
    """Write freqs and modes to a file in Gaussian-output-format.

       Arguments:
         | filename  --  modes are written to this file,
                         can be read by Molden (visualization program)
         | nma  --  modes information (see below)

       Optional argument:
         | selected  --  Selection of modes for which to make
                         trajectories. This can be a list or array
                         of mode indices (length <= N), or an
                         array of booleans (length = N).

       The nma arguments can have different formats:

       1) an NMA object
       2) a tuple or list with five elements: modes,
          frequencies, masses, numbers, coordinates
    """
    def parse_nma(nma):
        if isinstance(nma, NMA):
            # NMA object
            return nma.modes, nma.freqs, nma.masses, nma.numbers, nma.coordinates
        elif hasattr(nma, "__len__") and len(nma) == 5 and not isinstance(nma, numpy.ndarray):
            # [modes,freqs,...] or (modes,freqs,...)
            return nma
        else:
            raise TypeError("nma argument has wrong type")

    modes, freqs, masses, numbers, coordinates = parse_nma(nma)

    if selected is not None:
        modes = numpy.take(modes, selected, 1)  # modes in columns
        freqs = freqs[selected]
        masses = masses[selected]
        numbers = numbers[selected]
        coordinates = coordinates[selected]

    _make_moldenfile(filename, masses, numbers, coordinates, modes, freqs)


def _make_moldenfile(filename, masses, atomicnumbers, positions, modes, ev):
    """This function produces a molden-readable file: coordinates + frequencies + modes

    | positions  -- coordinates, convert to angstrom
    | modes  -- each col is a mode in mass weighted Cartesian coordinates
             un-mass-weighting necessary and renormalization (in order to see some movement)
    | ev  -- eigenvalues (freqs), convert to cm-1
    """
    masses3_sqrt1 = numpy.array(sum([[1/m,1/m,1/m] for m in numpy.sqrt(masses)],[]))
    HEAD, head_coordinates, head_basisfunctions, \
    head_freq0, head_freq1, head_freq2, head_freq3, head_end = _make_molden_texts()

    [rows,cols] = modes.shape
    number_of_atoms = rows/3
    number_of_modes = cols
    number_of_iterations = number_of_modes/3    # organisation of file: per 3 modes

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # start writing
    f = file(filename,"w+")

    print >> f, HEAD

    # ATOM PART
    print >> f, head_coordinates
    for at in range(number_of_atoms):
       print >> f, '%5d %10d %13s %15f %11f %11f' %(
                   at+1,atomicnumbers[at],"0",
                   positions[at,0]/angstrom,
                   positions[at,1]/angstrom,
                   positions[at,2]/angstrom)

    # ORBITAL PART
    print >> f, head_basisfunctions
    print >> f, " "   #this part is just empty

    # FREQUENCY PART
    print >> f, head_freq0

    for iteration in range(number_of_iterations):
        nb = 3*iteration   #number of mode
        mode1 = modes[:,nb]  *masses3_sqrt1
        mode2 = modes[:,nb+1]*masses3_sqrt1
        mode3 = modes[:,nb+2]*masses3_sqrt1
        mode1 = mode1/numpy.linalg.norm(mode1)
        mode2 = mode2/numpy.linalg.norm(mode2)
        mode3 = mode3/numpy.linalg.norm(mode3)
        print >> f, '%22d %22d %22d' %(nb+1,nb+2,nb+3)
        print >> f, head_freq1[2]
        print >> f, '%s %10.4f %22.4f %22.4f' %(head_freq2,
                                ev[nb]/lightspeed*centimeter,
                                ev[nb+1]/lightspeed*centimeter,
                                ev[nb+2]/lightspeed*centimeter)
        print >> f, head_freq3[2]
        for atomnb in range(number_of_atoms):
            i = 3*atomnb
            print >> f, '%4d %3d %8.2f %6.2f %6.2f %8.2f %6.2f %6.2f %8.2f %6.2f %6.2f' %(
                   atomnb+1 , atomicnumbers[atomnb],
                   mode1[i], mode1[i+1], mode1[i+2],
                   mode2[i], mode2[i+1], mode2[i+2],
                   mode3[i], mode3[i+1], mode3[i+2])

    rest = number_of_modes - 3*number_of_iterations

    if rest == 1:
        nb = number_of_modes-1   #number of mode: the last one
        mode1 = modes[:,nb]*masses3_sqrt1
        mode1 = mode1/numpy.linalg.norm(mode1)
        print >> f, '%22d' %(nb+1)
        print >> f, head_freq1[0]
        print >> f, '%s %10.4f' %(head_freq2, ev[nb]/lightspeed*centimeter)
        print >> f, head_freq3[0]
        for atomnb in range(number_of_atoms):
            i = 3*atomnb
            print >> f, '%4d %3d %8.2f %6.2f %6.2f' %(
                   atomnb+1 , atomicnumbers[atomnb],
                   mode1[i], mode1[i+1], mode1[i+2])

    elif rest == 2:
        nb = number_of_modes-2   #number of mode: the 2 last ones
        mode1 = modes[:,nb]  *masses3_sqrt1
        mode2 = modes[:,nb+1]*masses3_sqrt1
        mode1 = mode1/numpy.linalg.norm(mode1)
        mode2 = mode2/numpy.linalg.norm(mode2)
        print >> f, '%22d %22d' %(nb+1,nb+2)
        print >> f, head_freq1[1]
        print >> f, '%s %10.4f %22.4f' %(head_freq2,
                            ev[nb]  /lightspeed*centimeter,
                            ev[nb+1]/lightspeed*centimeter)
        print >> f, head_freq3[1]
        for atomnb in range(number_of_atoms):
            i = 3*atomnb
            print >> f, '%4d %3d %8.2f %6.2f %6.2f %8.2f %6.2f %6.2f' %(
                   atomnb+1 , atomicnumbers[atomnb],
                   mode1[i], mode1[i+1], mode1[i+2],
                   mode2[i], mode2[i+1], mode2[i+2],)

    elif rest != 0:
         print "error?! in number of iterations/number of atoms (writing molden file)"
    print >> f, head_end

    f.close()

def _make_molden_texts():

   HEAD = """ Entering Gaussian System
 this file is generated from the MLDGAU subroutine in the file secder.F
 Please note, that this is a "faked" output;
 there are no intensities computed in CPMD."""

   head_coordinates = """ Standard orientation:
 ---------------------------------------------------------------------
Center     Atomic     Atomic              Coordinates (Angstroms)
Number     Number      Type              X           Y           Z
 ---------------------------------------------------------------------"""

   head_basisfunctions = """ ---------------------------------------------------------------------
       basis functions          primitive gaussians
       alpha electrons          beta electrons
 **********************************************************************"""

   head_end = " Normal termination of Gaussian 98."

   head_freq0= """ Harmonic frequencies (cm**-1), IR intensities (KM/Mole),
 Raman scattering activities (A**4/AMU), Raman depolarization ratios,
 reduced masses (AMU), force constants (mDyne/A) and normal coordinates:"""

   head_freq1_1 =  "                    ?A"
   head_freq1_2 =  "                    ?A                     ?A"
   head_freq1_3 =  "                    ?A                     ?A                     ?A"
   head_freq1 = [head_freq1_1,head_freq1_2,head_freq1_3]

   head_freq2 =  " Frequencies --"

   head_freq3_1 = """ Red. masses --     0.0000
 Frc consts  --     0.0000
 IR Inten    --     0.0000
 Raman Activ --     0.0000
 Depolar     --     0.0000
 Atom AN      X      Y      Z"""
   head_freq3_2 = """ Red. masses --     0.0000                 0.0000
 Frc consts  --     0.0000                 0.0000
 IR Inten    --     0.0000                 0.0000
 Raman Activ --     0.0000                 0.0000
 Depolar     --     0.0000                 0.0000
 Atom AN      X      Y      Z        X      Y      Z"""
   head_freq3_3 = """ Red. masses --     0.0000                 0.0000                 0.0000
 Frc consts  --     0.0000                 0.0000                 0.0000
 IR Inten    --     0.0000                 0.0000                 0.0000
 Raman Activ --     0.0000                 0.0000                 0.0000
 Depolar     --     0.0000                 0.0000                 0.0000
 Atom AN      X      Y      Z        X      Y      Z        X      Y      Z"""
   head_freq3 = [head_freq3_1,head_freq3_2,head_freq3_3]

   return HEAD, head_coordinates, head_basisfunctions, \
          head_freq0, head_freq1, head_freq2, head_freq3, head_end


#======================================
#    (END) write logfile as Gaussian03 does
#======================================


