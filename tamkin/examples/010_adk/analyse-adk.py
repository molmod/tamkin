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

from __future__ import print_function
from tamkin import *


def do_NMA(filecor, filehess, filechk):
    # Load the data.
    molecule = load_molecule_charmm(filecor,filehess)

    # Perform the normal mode analysis
    nma = NMA(molecule)

    # Write to file
    nma.write_to_file(filechk)


def do_MBH(filecor, filehess, filefixed, filechk):
    # Load the data.
    molecule = load_molecule_charmm(filecor,filehess)
    fixed = load_blocks_txt(filefixed)

    # Perform the normal mode analysis
    nma = NMA(molecule, MBH(fixed))

    # Write to file
    nma.write_to_file(filechk)


def do_VSA(filecor, filehess, filefixed, filechk):
    # Load the data.
    molecule = load_molecule_charmm(filecor,filehess)
    subs = load_subs_txt(filefixed)

    # Perform the normal mode analysis
    nma = NMA(molecule, VSA(subs))

    # Write to file
    nma.write_to_file(filechk)


def take_cut(filechk):
    # Read the NMA object
    nma = NMA.read_from_file(filechk)


def do_gibbs(filechk, filethermo):
    nma = NMA.read_from_file(filechk)

    pf = PartFun(nma, [Vibrations(), ExtTrans(), ExtRot(1)] )

    ta = ThermoAnalysis(pf, [300.0])
    ta.write_to_file(filethermo)


def overlap(filecor1, filecor2, filechk, fileoverlap):

    # Get delta vector
    delta = get_delta_vector_charmmcor(filecor1, filecor2, normalize=True)

    # Read the NMA object
    nma = NMA.read_from_file(filechk)

    # Calculate overlap
    calculate_overlap(nma.modes,nma.freqs, delta,[0.0],filename=fileoverlap)


def get_delta(filecor1,filecor2,filedelta):
    get_delta_vector_charmmcor(filecor1, filecor2, normalize=True)



#-----------------------------------------
#------------ MAIN -----------------------

if __name__ == '__main__':
    from optparse import OptionParser

    usage = """
    python %prog FREQFILE [options]
    python %prog [options]

    This script ...
"""

    parser = OptionParser(usage=usage)

    parser.add_option('-o','--fig',dest='filename_fig',
                      default="intensity.png",
                      metavar="FIGNAME",
                      help="The name of the file that contains the plot, will be written.")

    parser.add_option('--job',dest='job',
                      default=None,
                      metavar="JOB",
                      help="Jobtype: nma, cut, ...")

    parser.add_option('--filecor', dest='filecor',
                      default = None,
                      metavar = "FILECOR",
                      help = "Name of file that contains the charmm coords.")

    parser.add_option('--filehess', dest='filehess',
                      default = None,
                      metavar = "FILEHESS",
                      help = "Name of file that contains the charmm Hessian.")

    parser.add_option('--filechk', dest='filechk',
                      default = None,
                      metavar = "FILECHK",
                      help = "Name of file that contains the checkpoint file made with Tamkin.")

    parser.add_option('--filecor2', dest='filecor2',
                      default = None,
                      metavar = "FILECOR2",
                      help = "Name of second file that contains the charmm coords.")

    parser.add_option('--filedelta', dest='filedelta',
                      default = None,
                      metavar = "FILEDELTA",
                      help = "Name of file that will contain the Delta (pacmanvector) between filecor and filecor2.")

    parser.add_option('--fileoverlap', dest='fileoverlap',
                      default = None,
                      metavar = "FILEOVERLAP",
                      help = "Name of file that will contain the overlap between the Delta (pacmanvector) and filecor (first cor file).")

    parser.add_option('--filefixed', dest='filefixed',
                      default = None,
                      metavar = "FILEFIXED",
                      help = "Name of file with partioning information eg blockchoice, environment atoms, ...")

    parser.add_option('--filethermo', dest='filethermo',
                      default = None,
                      metavar = "FILETHERMO",
                      help = "Name of file that will contain the thermodynamical analysis.")



    options, args = parser.parse_args()


    if options.job == "nma":
        if None in [options.filecor, options.filehess, options.filechk]:
            print("at least one of the variables filecor, filehess or filechk is empty")
            print("Doing nothing...")
        else:
            do_NMA(options.filecor, options.filehess, options.filechk)

    if options.job == "mbh":
        if None in [options.filecor, options.filehess, options.filefixed, options.filechk]:
            print("at least one of the variables filecor, filehess, filefixed or filechk is empty")
            print("Doing nothing...")
        else:
            do_MBH(options.filecor, options.filehess, options.filefixed, options.filechk)

    if options.job == "vsa":
        if None in [options.filecor, options.filehess, options.filefixed, options.filechk]:
            print("at least one of the variables filecor, filehess, filefixed or filechk is empty")
            print("Doing nothing...")
        else:
            do_VSA(options.filecor, options.filehess, options.filefixed, options.filechk)

    elif options.job == "cut":
        if options.filechk is None:
            print("variable filechk is empty")
            print("Doing nothing...")
        else:
            take_cut(options.filechk)

    elif options.job == "delta":
        if None in [options.filecor, options.filecor2, options.filedelta]:
            print("at least one of the variables filecor, filecor2 or filedelta is empty")
            print("Doing nothing...")
        else:
            get_delta(options.filecor, options.filecor2, options.filedelta)


    elif options.job == "overlap":
        if None in [options.filecor, options.filecor2, options.filechk, options.fileoverlap]:
            print("at least one of the variables filecor, filecor2, filechk or fileoverlap is empty")
            print("Doing nothing...")
        else:
            overlap(options.filecor, options.filecor2,
                    options.filechk, options.fileoverlap)

    elif options.job == "gibbs":
        if None in [options.filechk, options.filethermo]:
            print("at least one of the variables filechk or filethermo is empty")
            print("Doing nothing...")
        else:
            do_gibbs(options.filechk, options.filethermo)

    elif options.job is None:
        print("no job specified")
        print("Doing nothing...")
    else:
        print("The job "+options.job+" is not defined/implemented. Doing nothing...")
