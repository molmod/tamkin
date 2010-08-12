#!/usr/bin/python
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


from molmod.io.psf import PSFFile
from molmod.io.xyz import XYZFile
from molmod.molecular_graphs import MolecularGraph, HasAtomNumber, HasNeighbors, HasNeighborNumbers
from molmod.graphs import CritAnd


import sys


CC30A = CritAnd(HasAtomNumber(6), HasNeighborNumbers(6,6,6,6))
CC31A = CritAnd(HasAtomNumber(6), HasNeighborNumbers(6,6,6,1))
CC32A = CritAnd(HasAtomNumber(6), HasNeighborNumbers(6,6,1,1))
CC33A = CritAnd(HasAtomNumber(6), HasNeighborNumbers(6,1,1,1))

atom_filters = {
    "CC30A": CC30A,
    "CC31A": CC31A,
    "CC32A": CC32A,
    "CC33A": CC33A,
    "HCA1": CritAnd(HasAtomNumber(1), HasNeighbors(CC31A)),
    "HCA2": CritAnd(HasAtomNumber(1), HasNeighbors(CC32A)),
    "HCA3": CritAnd(HasAtomNumber(1), HasNeighbors(CC33A)),
}

def get_atom_type(index, graph):
    for atom_type, atom_filter in atom_filters.iteritems():
        if atom_filter(index, graph):
            return atom_type
    raise ValueError("Unrecognized atom (index %i)." % index)

args = sys.argv[1:]

molecule = XYZFile(args[0]).get_molecule()
graph = MolecularGraph.from_geometry(molecule)
atom_types = [get_atom_type(index, graph) for index in xrange(molecule.size)]

psf = PSFFile()
psf.add_molecular_graph(graph, atom_types=atom_types)
psf.write_to_file(args[0].replace(".xyz", ".psf"))
