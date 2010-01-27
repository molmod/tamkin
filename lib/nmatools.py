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
from tamkin.io import load_coordinates_charmm

from molmod import lightspeed, angstrom, amu, centimeter

import numpy, matplotlib, pylab


__all__ = [
           "calculate_overlap", "write_overlap",
           "compute_delta",
           "calculate_sensitivity_freq",
           "create_blocks_peptide_charmm", "create_subs_peptide_charmm",
           "BlocksPeptideMBH", "SubsPeptideVSA",
           "blocks_write_to_file", "selectedatoms_write_to_file",
           "plot_spectrum",
           "create_enm_molecule",
           "write_modes_for_VMD", "make_moldenfile_nma",
          ]




def calculate_overlap(nma1, nma2, filename=None):
    """Calculate overlap of modes and print to file if requested

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
    # calculate overlap
    overlap = numpy.dot(numpy.transpose(modes1), modes2)
    if filename is not None:
        write_overlap(freqs1, freqs2, overlap, filename=filename)
    return overlap


def write_overlap(freqs1, freqs2, overlap, filename=None):
    """Write overlap matrix to a file, default is overlap.csv. Format:
    ------------------------
           | freqs2
    ------------------------
    freqs1 | mat1^T . mat2
    ------------------------
    """
    #freqs1 = freqs1 /lightspeed*centimeter
    #freqs2 = freqs2 /lightspeed*centimeter

    # write to file
    if filename==None:
        filename="overlap.csv"   # TODO sys.currentdir

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
    """Calculate mass weighted delta vector between two conformations

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


def calculate_sensitivity_freq(nma, index, symmetric = False, massweight = True):
    """Calculate the sensity of the index-th frequency to changes of
    the mass-weighted Hessian elements.
    Optional:
    symmetric  --  Slightly different formula if symmetry of matrix is taken into account. Default False.
    massweight  --  Whether mass-weighted or un-mass-weighted Hessian is considered."""
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


#############################################################################
# In the following part: create atoms in blocks or atoms in subsystem for a
# standard peptide chain created in CHARMM
#
# Warning: first and last blocks are to be checked, because the atom order
# of the first and last residue are less predictable.


class BlocksPeptideMBH(object):
    # TODO add references

    def __init__(self, label=None, blocksize=1):
        self.label = label
        self.blocksize = blocksize    # the nb of residues per block, in RTB

    def calc_blocks(self, N, CA, PRO, Carbon, Oxigen, Nitrogen):
        if self.label is "RTB":
            return self.calc_blocks_RTB(N, CA, PRO, Carbon, Oxigen, Nitrogen)
        elif self.label is "dihedral":
            return self.calc_blocks_dihedral(N, CA, PRO, Carbon, Oxigen, Nitrogen)
        elif self.label is "RHbending":
            return self.calc_blocks_RHbending(N, CA, PRO, Carbon, Oxigen, Nitrogen)
        elif self.label is "normal" or self.label is None:
            return self.calc_blocks_normal(N, CA, PRO, Carbon, Oxigen, Nitrogen)
        else:
            raise NotImplementedError

    def calc_blocks_RTB(self, N, CA, PRO, Carbon, Oxigen, Nitrogen):
        """Rotation-Translation Blocks scheme of Tama et al.
        This amounts to applying the Mobile Block Hessian with
        one or more residues per block."""
        pept = []
        k = self.blocksize           # number of residues per block
        n = len(CA)/self.blocksize   # number of complete RTB-blocks

        # do for first block : ... (CA,R) x self.size - CO
        if CA[k-1] > 1:
            if CA[k-1] in PRO:
                pept.append(range(1,CA[k-1]-4))
            else:
                pept.append(range(1,CA[k-1]-2))

        # for next blocks :  (N - CA,R - CO) x self.size
        for i in xrange(1,n):   # do n-1 times
            # from first N till next N
            if CA[i*k-1] in PRO:
                if CA[(i+1)*k-1] in PRO:
                    pept.append(range(CA[i*k-1]-4,CA[(i+1)*k-1]-4))
                    #print "(side chain) proline found! (1,2)", CA[i*k-1],CA[(i+1)*k-1]
                else:
                    pept.append(range(CA[i*k-1]-4,CA[(i+1)*k-1]-2))
                    #print "(side chain) proline found! (1)", CA[i*k-1]
            else:
                if CA[(i+1)*k-1] in PRO:
                    pept.append(range(CA[i*k-1]-2,CA[(i+1)*k-1]-4))
                    #print "(side chain) proline found! (2)", CA[(i+1)*k-1]
                else:
                    pept.append(range(CA[i*k-1]-2,CA[(i+1)*k-1]-2))

        # for last block : N - CA,R
        if n*k-1 < len(CA):
            if CA[n*k-1] in PRO:
                pept.append(range(CA[n*k-1]-4,N+1))
                #print "(side chain) proline found! (1)", CA[n*k-1]
            else:
                pept.append(range(CA[n*k-1]-2,N+1))
        return pept


    def calc_blocks_dihedral(self, N, CA, PRO, Carbon, Oxigen, Nitrogen):
        """MBH scheme with linked blocks that selects the characteristic Phi and Psi
        dihedral angles of the protein backbone as variables, while other
        degrees of freedom are fixed. """
        pept = get_pept_linked(N, CA, PRO)
        res  = []

        # start with an ending :   ... CA_0,R - C
        if CA[1] in PRO:
            res.append( range(1,CA[1]-5) )
            #print "proline found!", CA[1]
        else:
            res.append( range(1,CA[1]-3) )
        # continue with normal residues : N - CA,R - C
        for i in xrange(1,len(CA)-1):
            if CA[i] in PRO:
                if CA[i+1] in PRO:
                    res.append( [CA[i]-4]+range(CA[i],CA[i+1]-5) )
                    #print "(side chain) proline found! (1,2)", CA[i],CA[i+1]
                else:
                    res.append( [CA[i]-4]+range(CA[i],CA[i+1]-3) )
                    #print "(side chain) proline found! (1)", CA[i]
            else:
                if CA[i+1] in PRO:
                    res.append( [CA[i]-2]+range(CA[i],CA[i+1]-5) )
                    #print "(side chain) proline found! (2)", CA[i+1]
                else:
                    res.append( [CA[i]-2]+range(CA[i],CA[i+1]-3) )
        # finish with another ending : N - CA,R
        if CA[-1] in PRO:
            res.append([CA[-1]-4]+range(CA[-1],N+1))
            #print "(side chain) proline found!", CA[-1]
        else:
            res.append([CA[-1]-2]+range(CA[-1],N+1))
        return res + pept


    def calc_blocks_RHbending(self, N, CA, PRO, Carbon, Oxigen, Nitrogen):
        """MBH scheme in which the CA-H bond is considered as a seperate block."""
        pept = get_pept_linked(N, CA, PRO)
        res  = []
        CH   = []

        # start with an ending
        if CA[1] in PRO:
           res.append( range(1,CA[1]-6) )
           #print "(side chain) proline found!", CA[1]
        else:
           res.append( range(1,CA[1]-4) )
        # continue with normal residues
        for i in xrange(1,len(CA)-1):  # CA and HA
            CH.append( [CA[i],CA[i]+1] )
        for i in xrange(1,len(CA)-1):  # CA and rest of residue
            if CA[i+1] in PRO:
                res.append( [CA[i]]+range(CA[i]+2,CA[i+1]-6) )
                #print "(side chain) proline found!", CA[i+1]
            else:
                res.append( [CA[i]]+range(CA[i]+2,CA[i+1]-4) )
        # finish with another ending
        CH.append([CA[-1],CA[-1]+1]) # CA and HA
        res.append([CA[-1]]+range(CA[-1]+2,N+1)) # CA and rest of residue

        return res + pept + CH

    def calc_blocks_normal(self, N, CA, PRO, Carbon, Oxigen, Nitrogen):
        """MBH scheme with linked blocks where each side chain is
        considered as a block attached to the CA atom."""
        pept = get_pept_linked(N, CA, PRO)
        res  = []

        # start with an ending
        if CA[1] in PRO:
           res.append( range(1,CA[1]-6) )
           #print "(side chain) proline found!", CA[1]
        else:
           res.append( range(1,CA[1]-4) )
        # continue with normal residues
        for i in xrange(1,len(CA)-1):
            if CA[i+1] in PRO:
                res.append( range(CA[i],CA[i+1]-6) )
                #print "(side chain) proline found!", CA[i+1]
            else:
                res.append( range(CA[i],CA[i+1]-4) )
        # finish with another ending
        res.append(range(CA[-1],N+1))

        return res + pept


def get_pept_linked(N, CA, PRO):
        # PEPT bonds = CA + CONH + CA = 6 atoms
        pept = []
        for i in xrange(1,len(CA)):
            if CA[i] in PRO:
                pept.append( [CA[i-1]] + range(CA[i]-6,CA[i]+1) )
                #print "proline found!", CA[i]
            else:
                pept.append( [CA[i-1]] + range(CA[i]-4,CA[i]+1) )
        return pept



class SubsPeptideVSA(object):

    def __init__(self, atomtype=["CA"], frequency=1):
        """VSA subsystem/environment for peptides.
        Optional:
        atomtype  --  list of strings. Let only these atom types be part of the subsystem.
        frequency  --  let only one out of every *frequency* residues participate"""
        self.atomtype = atomtype
        self.frequency = frequency

    def calc_subs(self, N, CA, PRO, Carbon, Oxigen, Nitrogen):
        subs = []
        if "C" in self.atomtype:
            subs.extend(Carbon)
        if "O" in self.atomtype:
            subs.extend(Oxigen)
        if "N" in self.atomtype:
            subs.extend(Nitrogen)
        if "CA" in self.atomtype:
            subs.extend( numpy.take(CA, range(0,len(CA),self.frequency)).tolist() )
        else:
            raise NotImplementedError
        return subs


def create_subs_peptide_charmm(charmmfile_crd, blockchoice):
    N, CA, PRO, Carbon, Oxigen, Nitrogen = select_info_peptide_charmm(charmmfile_crd)
    return blockchoice.calc_subs(N, CA, PRO, Carbon, Oxigen, Nitrogen)


def create_blocks_peptide_charmm(charmmfile_crd, blockchoice):
    N, CA, PRO, Carbon, Oxigen, Nitrogen = select_info_peptide_charmm(charmmfile_crd)
    return blockchoice.calc_blocks(N, CA, PRO, Carbon, Oxigen, Nitrogen)


def select_info_peptide_charmm(charmmfile_crd):
    # Reading from charmmfile
    f = file(charmmfile_crd)
    # nb of atoms
    for i,line in enumerate(f):
        words = line.split()
        if words[0]!="*":
            N = int(words[0])
            break
    # find CA atoms, PRO residues, Carbons, Oxigens, Nitrogens
    PRO = []
    CA = []
    Carbon = []
    Oxigen = []
    Nitrogen = []
    for i,line in enumerate(f):
        words = line.split()
        if words[3].startswith("CA"):
            CA.append(int(words[0]))
            if words[2]=="PRO":
                PRO.append(int(words[0]))
        if words[3]=="C":
            Carbon.append(int(words[0]))
        if words[3]=="O":
            Oxigen.append(int(words[0]))
        if words[3]=="N":
            Nitrogen.append(int(words[0]))
    f.close()
    return N, CA, PRO, Carbon, Oxigen, Nitrogen


def blocks_write_to_file(blocks, filename, shift=1):
    """write atoms in blocks to file.
    One atom per line, a blank line starts a new block.
    Optional
    shift  --  write atom+shift to file.
               Default is 1, because default shift in load_subs_txt and load_fixed_txt is -1."""
    f = file(filename, "w")
    for bl in blocks:
        for at in bl:
            print >> f, at+shift
        print >> f, ""
    f.close()

def selectedatoms_write_to_file(selected, filename, shift=1):
    """write selected atoms to file, e.g. subsystem or environment atoms.
    One atom per line.
    Optional:
    shift  --  write atom+shift to file.
               Default is 1, because default shift in load_subs_txt and load_fixed_txt is -1."""
    f = file(filename, "w")
    for at in selected:
        print >> f, at+shift
    print >> f, ""
    f.close()



def plot_spectrum(label, filename, freqlist, min=None, max=None, Imax=None,
                  step=1.0, width=10.0, amplitudes=None, title=None):
    """
    Plot the spectra of the frequencies in the freqlist: each item in freqlist
    leads to an additional spectrum.
    * label dos (density of states): sum of Gaussians each centered around
                                     a certain freq, with certain width
    * label lines (line spectrum): draw a horizontal line for each freq.

    Options
    min  --  min on x-axis, in cm-1
    max  --  max on x-axis, in cm-1
    Imax  --  maximum intensity on y-axis, no unit
    step  --  resulotion of plot, in cm-1
    width  --  width of Gaussian, in cm-1
    amplitudes  --  if different amplitudes for the items in freqlist
    title  --  title for plot (a string)
    """
    if label=="dos":   # plot density of states spectrum
        fig = pylab.figure()
        fig.hold(True)
        ax = fig.add_subplot(1,1,1)
        # do plotting
        for i,freqs in enumerate(freqlist):
            if amplitudes is not None: ampl = amplitudes[i]
            else: ampl = 1.0
            do_plot_dos(ax, freqs /lightspeed*centimeter, min=min, max=max,
                              step=step, width=width, amplitude=ampl)
        # clean up
        if Imax is not None:  ax.set_ylim(0.0,IMax)
        if min is not None:   ax.set_xlim(xmin=min)
        if max is not None:   ax.set_xlim(xmax=max)
        pylab.legend([str(nb) for nb in range(1,len(freqlist)+1)])
        pylab.ylabel("intensity")
        if title is not None: pylab.title(title)
        fig.savefig(filename)
        pylab.close()

    if label=="lines":   # plot line spectrum
        pylab.figure()
        pylab.hold(True)
        # do plotting
        for i,freqs in enumerate(freqlist):
            for freq in freqs /lightspeed*centimeter :
                if (freq>min or min is None) and (freq<max or max is None):
                    pylab.plot([i+0.75,i+1.25],[freq,freq],"k-")
        # clean up
        pylab.xticks(range(1,len(freqlist)+1))
        if min is not None:  pylab.ylim(ymin=min)
        if max is not None:  pylab.ylim(ymax=max)
        pylab.ylabel("freq in cm-1")
        if title is not None: pylab.title(title)
        pylab.savefig(filename)
        pylab.close()

def do_plot_dos(ax, freqs, min, max,
                step=1.0, width=10.0, amplitude=1.0):

    if min is None:  min=freqs[0]-3*width   # to avoid unnecessary plotting
    if max is None:  max=freqs[-1]+3*width

    freqs = numpy.array(freqs)
    nu = numpy.arange(min, max, step)
    intensity = numpy.zeros(nu.shape,float)

    s2 = width**2 / ( 8*numpy.log(2) )  # standard deviation squared
    for freq in freqs:
        if freq>min and freq<max:
            intensity += numpy.exp(-(nu-freq)**2/(2*s2))
    ax.plot(nu,intensity*amplitude)


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
                       be a list or array of atom indices (length < N), or an
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


def write_modes_for_VMD(nma, index, filename=None,
                        A = 50.0, frames = 36):
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


#======================================
#    write logfile as Gaussian03 does
#======================================

# this file contains all texts needed to generate
# a .log file that can be read by molden (visualization program)
# - def make_molden_texts()
# - def make_moldenfile(filename, masses, atomicnumbers, positions, modes, ev)

def make_molden_texts():

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


def make_moldenfile_nma(filename, nma):
    """Write a logfile with the modes and frequencies
    in the way Gaussian03 does, such that molden can
    read this logfile.
    """
    if filename is None:
        filename = "molden.log"
    if nma.modes is None:
        raise ValueError("No modes available (do_modes=False), cannot write logfile with modes.")
    make_moldenfile(filename, nma.masses, nma.numbers, nma.coordinates,
                    nma.modes, nma.freqs )



def make_moldenfile(filename, masses, atomicnumbers, positions, modes, ev):
    """This function produces a molden-readable file: coordinates + frequencies + modes

    positions  -- coordinates, convert to angstrom
    modes  -- each col is a mode in mass weighted Cartesian coordinates
             un-mass-weighting necessary and renormalization (in order to see some movement)
    ev  -- eigenvalues (freqs), convert to cm-1
    """
    masses3_sqrt1 = numpy.array(sum([[1/m,1/m,1/m] for m in numpy.sqrt(masses)],[]))
    HEAD, head_coordinates, head_basisfunctions, \
    head_freq0, head_freq1, head_freq2, head_freq3, head_end = make_molden_texts()

    [rows,cols]=modes.shape
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

#======================================
#    (END) write logfile as Gaussian03 does
#======================================

