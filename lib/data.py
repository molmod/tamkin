# TAMkin is a post-processing toolkit for thermochemistry and kinetics analysis.
# Copyright (C) 2008-2009 Toon Verstraelen <Toon.Verstraelen@UGent.be>,
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


from tamkin.geom import transrot_basis

from molmod.molecules import Molecule as BaseMolecule
from molmod.molecular_graphs import MolecularGraph
from molmod.data.periodic import periodic
from molmod.graphs import cached

import numpy


__all__ = ["Molecule", "BareNucleus", "Proton", "RotScan"]


class Molecule(BaseMolecule):
    def __init__(self, numbers, coordinates, masses, energy, gradient, hessian, multiplicity, symmetry_number, periodic):
        BaseMolecule.__init__(self, numbers, coordinates)
        self._masses = numpy.array(masses, float)
        self._masses.setflags(write=False)
        self._energy = energy
        self._gradient = numpy.array(gradient, float)
        self._gradient.setflags(write=False)
        self._hessian = numpy.array(hessian, float)
        self._hessian.setflags(write=False)
        self._multiplicity = multiplicity
        self._symmetry_number = symmetry_number
        self._periodic = periodic

    masses = property(lambda self: self._masses)
    energy = property(lambda self: self._energy)
    gradient = property(lambda self: self._gradient)
    hessian = property(lambda self: self._hessian)
    multiplicity = property(lambda self: self._multiplicity)
    symmetry_number = property(lambda self: self._symmetry_number)
    periodic = property(lambda self: self._periodic)

    @cached
    def mass(self):
        return self.masses.sum()

    @cached
    def com(self):
        return (self.coordinates*self.masses.reshape((-1,1))).sum(axis=0)/self.mass

    @cached
    def inertia_tensor(self):
        return sum(
            m*(numpy.identity(3)*(r**2).sum()-numpy.outer(r,r))
            for m, r in zip(self.masses, (self.coordinates-self.com))
        )

    @cached
    def chemical_formula(self):
        counts = {}
        for number in self.numbers:
            counts[number] = counts.get(number, 0)+1
        return " ".join(
            "%s%i" % (periodic[number].symbol, count)
            for number, count in sorted(counts.iteritems(), reverse=True)
        )

    @cached
    def external_basis(self):
        """The basis for small displacements in the external degrees of freedom.
        """
        result = transrot_basis(self.coordinates, not self.periodic)
        result *= numpy.sqrt(self.masses3) # transform basis to mass weighted coordinates
        return result

    @cached
    def masses3(self):
        return numpy.array([self.masses, self.masses, self.masses]).transpose().ravel()


class BareNucleus(Molecule):
    def __init__(self, number, mass=None):
        if mass is None:
            mass = periodic[number].mass
        Molecule.__init__(self, [number], [[0,0,0]], [mass], 0.0, [[0,0,0]], [[0,0,0],[0,0,0],[0,0,0]], 1, None, False)


class Proton(BareNucleus):
    def __init__(self, mass=None):
        BareNucleus.__init__(self, 1, mass)


class RotScan(object):
    def __init__(self, dihedral, molecule=None, top_indexes=None, potential=None):
        """Initialize a rotational scan object

           Arguments
             dihedral  --  the index of the atoms that define the dihedral angle

           Optional arguments
             molecule  --  a molecule object. required when top_indexes is not
                           given
             top_indexes  --  a list of atom indexes involved in the rotor.
                              required when molecule is not given
             potential  --  rotational potential info (if this is a hindered
                            rotor). must be a two-tuple containing the angles
                            and the corresponding energies.

        """
        if len(dihedral) != 4:
            raise ValueError("The first argument must be a list of 4 integers")
        self.dihedral = dihedral
        if top_indexes is None:
            # try to deduce the top indexes
            if molecule is None:
                raise ValueError("missing arguemt: top_indexes or molecule from which top_indexes can be derived")
            atom0, atom1 = dihedral[1:3]
            graph = MolecularGraph.from_geometry(molecule)
            half0, half1 = graph.get_halfs(atom0, atom1)
            if len(half1) > len(half0):
                top_indexes = half0
                top_indexes.discard(atom0)
            else:
                top_indexes = half1
                top_indexes.discard(atom1)
            self.top_indexes = list(top_indexes)
        else:
            self.top_indexes = top_indexes
        if len(self.top_indexes) < 1:
            raise ValueError("A rotational scan must have at least one atom rotating")
        self.potential = potential


