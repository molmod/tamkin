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


from tamkin.data import Molecule
from tamkin.nma import NMA
from tamkin.io.charmm import load_peptide_info_charmm

from molmod import lightspeed, angstrom, amu, centimeter

import numpy, pylab


__all__ = [
    "compute_overlap", "write_overlap",
    "compute_delta", "compute_sensitivity_freq",
    "create_blocks_peptide_charmm", "create_subs_peptide_charmm",
    "plot_spectrum_lines", "plot_spectrum_dos",
    "create_enm_molecule",
]


invcm = lightspeed/centimeter


def compute_overlap(nma1, nma2, filename=None):
    """Compute overlap of modes and print to file if requested

       Arguments:
         nma1  --  modes and frequencies (see below)
         nma2  --  modes and frequencies (see below)

       Optional argument:
         filename  --  when given, the overlap is written to file by the
                       function write_overlap

       The nma arguments can have different formats:
       1) an NMA object
       2) a tuple or list with two elements: modes and frequencies
       3) a numpy array with the mass-weighted modes
       4) a numpy array with one mass-weighted mode
    """

    def parse_nma(nma):
        if isinstance(nma, NMA):
            # NMA object
            return nma.modes, nma.freqs
        elif hasattr(nma, "__len__") and len(nma) == 2 and not isinstance(nma, numpy.ndarray):
            # [modes,freqs] or (modes,freqs)
            return nma
        elif isinstance(nma, numpy.ndarray) and len(nma.shape) == 2:
            # modes only
            return nma, numpy.zeros(nma.shape[1], float)
        elif isinstance(nma, numpy.ndarray) and len(nma.shape) == 1:
            # one mode only
            return nma.reshape(-1,1), numpy.zeros(1, float)
        else:
            raise TypeError("nma argument has wrong type")

    modes1, freqs1 = parse_nma(nma1)
    modes2, freqs2 = parse_nma(nma2)

    # check dimensions
    if modes1.shape[0] != modes2.shape[0] :
        raise ValueError("Length of columns in modes1 and modes2 should be equal, but found %i and %i." % (modes1.shape[0], modes2.shape[0]))
    # compute overlap
    overlap = numpy.dot(numpy.transpose(modes1), modes2)
    if filename is not None:
        write_overlap(freqs1, freqs2, overlap, filename=filename)
    return overlap


def write_overlap(freqs1, freqs2, overlap, filename="overlap.csv"):
    """Write the overlap matrix to a csv file

       Arguments:
         freqs1  --  the list of frequencies associated with the rows of the
                     overlap matrix
         freqs2  --  the list of frequencies associated with the columns of the
                     overlap matrix
         overlap  --  the overlap matrix

       Optional arguments:
         filename  --  the file to write to [default="overlap.csv"]

       The table is contains the following blocks:
        ------------------------
                 | freqs2
        ------------------------
        freqs1^T | mat1^T . mat2
        ------------------------
    """
    #freqs1 = freqs1 / invcm
    #freqs2 = freqs2 / invcm

    to_append="w+"   # not append, just overwrite
    f = file(filename,to_append)

    [rows,cols] = overlap.shape

    # 1. row of freqs2
    print >> f, ";"+";".join(str(g) for g in freqs2)  #this is the same

    # 2. start each row with freq of freqs1 and continue with overlaps
    for r in range(rows):
        print >> f, str(freqs1[r])+";"+";".join(str(g) for g in overlap[r,:].tolist())
    f.close()


def compute_delta(coor1, coor2, masses=None, normalize=False):
    """Compute mass weighted delta vector between two conformations

       It is assumed that the structures have been aligned (center of mass,
       orientation) previously.

       Arguments:
         coor1  --  coordinates of structure 1 in a numpy array with shape (N,3)
         coor2  --  coordinates of structure 2 in a numpy array with shape (N,3)

       Optional arguments:
         masses  --  when given, the mass-weighted delta vector is computed
         normalize  --  whether delta vector should be normalized [default=False]
    """
    # check consistency
    if len(coor1) != len(coor2):
        raise ValueError("coordinates should have same length: found %i and %i" % (len(coor1), len(coor2)))

    delta = numpy.ravel(coor1 - coor2)
    if not masses is None:  #Mass-weighting delta vector
        for i,mass in enumerate(masses):
            delta[3*i:3*(i+1)] *=  numpy.sqrt(mass)
    if normalize:   #Normalizing delta vector
        norm = numpy.sum(delta**2)
        delta /= numpy.sqrt(norm)
    return numpy.reshape(delta, (-1,1))


def compute_sensitivity_freq(nma, index, symmetric=False, massweight=True):
    """Compute the sensity of the index-th frequency to changes in
       the mass-weighted Hessian elements.

       Arguments:
         nma  --  an NMA object

       Optional argumets:
         symmetric  --  when True, a slightly different formula is used to take
                        into account the symmetry of the Hessian [default=False]
         massweight  --  when True, a mass-weighted hessian is considered
    """
    L = 3*len(nma.masses)
    mode = nma.modes[:,index]
    if not massweight: # un-mass-weight the mode
        for at,mass in enumerate(nma.masses):
            mode[3*at:3*(at+1)] /= numpy.sqrt(mass)
        mode = mode / numpy.sqrt(numpy.sum(mode**2))  # renormalization necessary

    mat = numpy.dot( numpy.reshape(mode,(L,1)), numpy.reshape(mode,(1,L)) )
    if symmetric:
        mat *= 2
        for i in range(L):
             mat[i,i] -= mode[i]**2
    return mat


def create_blocks_peptide_charmm(filename, label="normal", blocksize=1):
    """Create blocks list for CHARMM peptides

       Argument:
         filename  --  the CHARMM coordinate file (typically extension .crd or
                       .cor)

       Optional argument:
         label  --  type of MBH blocks: RTB, dihedral, RHbending, normal
                    [default=normal]
         blocksize  --  when using the RTB scheme, blocksize defines the number
                        of residues in a block

       TODO: referenties
    """
    N, calpha, proline, carbon, oxygen, nitrogen = load_peptide_info_charmm(filename)

    if label is "RTB":
        return _calc_blocks_RTB(blocksize, N, calpha, proline, carbon, oxygen, nitrogen)
    elif label is "dihedral":
        return _calc_blocks_dihedral(N, calpha, proline, carbon, oxygen, nitrogen)
    elif label is "RHbending":
        return _calc_blocks_RHbending(N, calpha, proline, carbon, oxygen, nitrogen)
    elif label is "normal":
        return _calc_blocks_normal(N, calpha, proline, carbon, oxygen, nitrogen)
    else:
        raise NotImplementedError


def _get_pept_linked(N, calpha, proline):
    """Defines the peptide blocks [CA,C,O,N,NH,CA]

       The alpha carbons are shared between neighboring blocks. The routine
       relies on the conventional atom order in CHARMM peptide files.
    """
    # PEPT bonds = calpha + CONH + calpha = 6 atoms
    pept = []
    for i in xrange(1,len(calpha)):
        if calpha[i] in proline:
            pept.append( [calpha[i-1]] + range(calpha[i]-6,calpha[i]+1) )
            #print "proline found!", calpha[i]
        else:
            pept.append( [calpha[i-1]] + range(calpha[i]-4,calpha[i]+1) )
    return pept


def _calc_blocks_RTB(blocksize, N, calpha, proline, carbon, oxygen, nitrogen):
    """Rotation-Translation Blocks scheme of Tama et al

       Arguments: see functions 'create_blocks_peptide_charmm' and
                  'load_peptide_info_charmm'

       This amounts to applying the Mobile Block Hessian with
       one or more residues per block.
    """
    pept = []
    k = blocksize               # number of residues per block
    n = len(calpha)/blocksize   # number of complete RTB-blocks

    # do for first block : ... (calpha,R) x size - CO
    if calpha[k-1] > 1:
        if calpha[k-1] in proline:
            pept.append(range(1,calpha[k-1]-4))
        else:
            pept.append(range(1,calpha[k-1]-2))

    # for next blocks :  (N - calpha,R - CO) x size
    for i in xrange(1,n):   # do n-1 times
        # from first N till next N
        if calpha[i*k-1] in proline:
            if calpha[(i+1)*k-1] in proline:
                pept.append(range(calpha[i*k-1]-4,calpha[(i+1)*k-1]-4))
                #print "(side chain) proline found! (1,2)", calpha[i*k-1],calpha[(i+1)*k-1]
            else:
                pept.append(range(calpha[i*k-1]-4,calpha[(i+1)*k-1]-2))
                #print "(side chain) proline found! (1)", calpha[i*k-1]
        else:
            if calpha[(i+1)*k-1] in proline:
                pept.append(range(calpha[i*k-1]-2,calpha[(i+1)*k-1]-4))
                #print "(side chain) proline found! (2)", calpha[(i+1)*k-1]
            else:
                pept.append(range(calpha[i*k-1]-2,calpha[(i+1)*k-1]-2))

    # for last block : N - calpha,R
    if n*k-1 < len(calpha):
        if calpha[n*k-1] in proline:
            pept.append(range(calpha[n*k-1]-4,N+1))
            #print "(side chain) proline found! (1)", calpha[n*k-1]
        else:
            pept.append(range(calpha[n*k-1]-2,N+1))
    return pept


def _calc_blocks_dihedral(N, calpha, proline, carbon, oxygen, nitrogen):
    """MBH scheme with linked blocks that selects the characteristic Phi and Psi
       dihedral angles of the protein backbone as variables, while other
       degrees of freedom are fixed.

       Arguments: see function 'load_peptide_info_charmm'
    """
    pept = _get_pept_linked(N, calpha, proline)
    res  = []

    # start with an ending :   ... calpha_0,R - C
    if calpha[1] in proline:
        res.append( range(1,calpha[1]-5) )
        #print "proline found!", calpha[1]
    else:
        res.append( range(1,calpha[1]-3) )
    # continue with normal residues : N - calpha,R - C
    for i in xrange(1,len(calpha)-1):
        if calpha[i] in proline:
            if calpha[i+1] in proline:
                res.append( [calpha[i]-4]+range(calpha[i],calpha[i+1]-5) )
                #print "(side chain) proline found! (1,2)", calpha[i],calpha[i+1]
            else:
                res.append( [calpha[i]-4]+range(calpha[i],calpha[i+1]-3) )
                #print "(side chain) proline found! (1)", calpha[i]
        else:
            if calpha[i+1] in proline:
                res.append( [calpha[i]-2]+range(calpha[i],calpha[i+1]-5) )
                #print "(side chain) proline found! (2)", calpha[i+1]
            else:
                res.append( [calpha[i]-2]+range(calpha[i],calpha[i+1]-3) )
    # finish with another ending : N - calpha,R
    if calpha[-1] in proline:
        res.append([calpha[-1]-4]+range(calpha[-1],N+1))
        #print "(side chain) proline found!", calpha[-1]
    else:
        res.append([calpha[-1]-2]+range(calpha[-1],N+1))
    return res + pept


def _calc_blocks_RHbending(N, calpha, proline, carbon, oxygen, nitrogen):
    """MBH scheme in which the calpha-H bond is considered as a seperate block.

       Arguments: see function 'load_peptide_info_charmm'
    """
    pept = _get_pept_linked(N, calpha, proline)
    res  = []
    CH   = []

    # start with an ending
    if calpha[1] in proline:
        res.append( range(1,calpha[1]-6) )
        #print "(side chain) proline found!", calpha[1]
    else:
        res.append( range(1,calpha[1]-4) )
    # continue with normal residues
    for i in xrange(1,len(calpha)-1):  # calpha and HA
        CH.append( [calpha[i],calpha[i]+1] )
    for i in xrange(1,len(calpha)-1):  # calpha and rest of residue
        if calpha[i+1] in proline:
            res.append( [calpha[i]]+range(calpha[i]+2,calpha[i+1]-6) )
            #print "(side chain) proline found!", calpha[i+1]
        else:
            res.append( [calpha[i]]+range(calpha[i]+2,calpha[i+1]-4) )
    # finish with another ending
    CH.append([calpha[-1],calpha[-1]+1]) # calpha and HA
    res.append([calpha[-1]]+range(calpha[-1]+2,N+1)) # calpha and rest of residue

    return res + pept + CH


def _calc_blocks_normal(N, calpha, proline, carbon, oxygen, nitrogen):
    """MBH scheme with linked blocks where each side chain is considered as a
       block attached to the calpha atom.

       Arguments: see function 'load_peptide_info_charmm'
    """
    pept = _get_pept_linked(N, calpha, proline)
    res  = []

    # start with an ending
    if calpha[1] in proline:
        res.append( range(1,calpha[1]-6) )
        #print "(side chain) proline found!", calpha[1]
    else:
        res.append( range(1,calpha[1]-4) )
    # continue with normal residues
    for i in xrange(1,len(calpha)-1):
        if calpha[i+1] in proline:
            res.append( range(calpha[i],calpha[i+1]-6) )
            #print "(side chain) proline found!", calpha[i+1]
        else:
            res.append( range(calpha[i],calpha[i+1]-4) )
    # finish with another ending
    res.append(range(calpha[-1],N+1))

    return res + pept


def create_subs_peptide_charmm(filename, atomtypes=["CA"], frequency=1):
    """Create subsystem selection for CHARMM peptides

       Argument:
         filename  --  the CHARMM coordinate file (typically extension .crd or
                       .cor)

       Optional argument:
         atomtypes  --  list of strings. Let only these atom types be part of
                        the subsystem.
         frequency  --  let only one out of every *frequency* residues be part
                        of the subsystem.
    """
    N, calpha, proline, carbon, oxygen, nitrogen = load_peptide_info_charmm(filename)

    subs = []
    if "C" in atomtypes:
        subs.extend(carbon[::frequency])
    if "O" in atomtypes:
        subs.extend(oxygen[::frequency])
    if "N" in atomtypes:
        subs.extend(nitrogen[::frequency])
    if "CA" in atomtypes:
        subs.extend(calpha[::frequency])
    if len(subs) == 0:
        raise NotImplementedError("No matching atom types found.")
    return subs


def plot_spectrum_lines(filename, all_freqs, low=None, high=None, title=None):
    """Plot multiple spectra in a comparative line plot

       Arguments:
         filename  --  the filename to write the figure too (the extension and
                       the matplotlib settings determine the file format)
         all_freqs  --  a list with spectra, each item in the list is an array
                       with multiple frequencies that represent one spectrum

       Optional arguments:
         low  --  minimum on x-axis, in atomic units
         high  --  maximum on x-axis, in atomic units
         title  --  title for plot (a string)
    """
    pylab.clf()
    for i, freqs in enumerate(all_freqs):
        for freq in freqs:
            if (freq>low or low is None) and (freq<high or high is None):
                pylab.plot([i+0.75,i+1.25],[freq/invcm,freq/invcm],"k-")
    pylab.xticks(range(1,len(all_freqs)+1))
    if low is not None:
        pylab.ylim(ymin=low/invcm)
    if high is not None:
        pylab.ylim(ymax=high/invcm)
    pylab.ylabel("Frequency in cm-1")
    if title is not None:
        pylab.title(title)
    pylab.savefig(filename)


def plot_spectrum_dos(filename, all_freqs, low=None, high=None, imax=None,
                      step=1.0*invcm, width=10.0*invcm, all_amps=None, title=None):
    """Plot multiple spectra in a comparative density of states plot

       Arguments:
         filename  --  the filename to write the figure too (the extension and
                       the matplotlib settings determine the file format)
         all_freqs  --  a list with spectra, each item in the list is an array
                       with multiple frequencies that represent one spectrum

       Optional arguments:
         low  --  minimum on x-axis, in atomic units
         high  --  maximum on x-axis, in atomic units
         imax  --  maximum intensity on y-axis, no unit
         step  --  resulotion of plot, in atomic units
         width  --  width of Gaussian, in atomic units
         all_amps  --  list of arrays in the same format as all_freqs with an
                       amplitude for each individual frequency
         title  --  title for plot (a string)
    """

    def plot_single_dos(freqs, low, high, step, width, amps):
        if low is None:
            low=freqs[0]-3*width   # to avoid unnecessary plotting
        if high is None:
            high=freqs[-1]+3*width

        nu = numpy.arange(low, high, step)
        intensity = numpy.zeros(nu.shape, float)

        s2 = width**2 / ( 8*numpy.log(2) )  # standard deviation squared
        for i in xrange(len(freqs)):
            if freqs[i] > low and freqs[i] < high:
                tmp = numpy.exp(-(nu-freqs[i])**2/(2*s2))
                if amps is not None:
                    if hasattr(amps, "__len__"):
                        tmp *= amps[i]
                    else:
                        tmp *= amps
                intensity += tmp
        pylab.plot(nu/invcm, intensity)

    pylab.clf()
    for i, freqs in enumerate(all_freqs):
        if all_amps is None:
            amps = None
        else:
            amps = all_amps[i]
        plot_single_dos(freqs, low, high, step, width, amps)
    pylab.ylim(0.0, imax)
    if low is not None:
        pylab.xlim(xmin=low/invcm)
    if high is not None:
        pylab.xlim(ymax=high/invcm)
    pylab.legend([str(nb) for nb in range(1,len(all_freqs)+1)])
    pylab.xlabel("Frequency in cm-1")
    pylab.ylabel("Intensity")
    if title is not None:
        pylab.title(title)
    pylab.savefig(filename)


def create_enm_molecule(molecule, selected=None, numbers=None, masses=None,
                        rcut=8.0*angstrom, K=1.0, periodic=None):
    """Create a molecule according to the Elastic Network Model

       Argument:
         molecule  --  The molecule to start from. can be two types: (i) a
                       Molecule object or (ii) a numpy array with shape (N,3)
                       with coordinates in atomic units.

       When a Molecule object is given, atom numbers, masses and periodic are
       inherited from the molecule, unless they are specified explicitly in the
       optional arguments.

       Optional arguments:
         selected  --  Selection of atoms to include in the ENM model. This can
                       be a list or array of atom indices (length <= N), or an
                       array of booleans (length = N).
         numbers  --  atom numbers in the ENM model (length = N). default is
                      array of ones or the numbers from the molecule object.
         masses  --  atomic masses in atomic units in the ENM model (length = N).
                     default is array of hydrogen masses or the masses from the
                     molecule object.
         rcut  --  cutoff distance between interacting pairs in atomic units
         K  --  strength of the interaction in atomic units (Hartree/Bohr**2).
                The interaction strength is the same for all interacting pairs.
    """
    if isinstance(molecule, Molecule):
        coordinates = molecule.coordinates
        if numbers is None:
            numbers = molecule.numbers
        if masses is None:
            masses = molecule.masses
        if periodic is None:
            periodic = molecule.periodic
    else:
        coordinates = numpy.array(molecule, copy=False)
        if numbers is None:
            # pretend being hydrogens
            numbers = numpy.ones(len(coordinates))
        if masses is None:
            # pretend being hydrogens
            masses  = numpy.ones(len(coordinates),float)*amu
        if periodic is None:
            periodic = False

    if selected is not None:
        coordinates = coordinates[selected]
        numbers = numbers[selected]
        masses = masses[selected]

    rcut2 = rcut**2
    N = len(coordinates)
    hessian = numpy.zeros((3*N,3*N),float)
    for i in range(N):
        for j in range(i+1,N):
            x = numpy.reshape( (coordinates[i,:]-coordinates[j,:]), (3,1))
            dist2 = numpy.sum(x**2)
            if dist2 < rcut2:
                corr = K * numpy.dot( x, x.transpose() ) / dist2
                hessian[3*i:3*(i+1), 3*i:3*(i+1)] += corr
                hessian[3*i:3*(i+1), 3*j:3*(j+1)] -= corr
                hessian[3*j:3*(j+1), 3*i:3*(i+1)] -= corr
                hessian[3*j:3*(j+1), 3*j:3*(j+1)] += corr

    return Molecule(
        numbers,
        coordinates,
        masses,
        0.0, #energy
        numpy.zeros((len(coordinates),3), float),  #gradient
        hessian,
        1, # multiplicity
        1, # rotational symmetry number
        periodic,
    )



