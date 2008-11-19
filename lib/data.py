# TAMkin is a post-processing toolkit for thermochemistry and kinetics analysis.
# Copyright (C) 2008 Toon Verstraelen <Toon.Verstraelen@UGent.be>,
# Matthias Vandichel <Matthias.Vandichel@UGent.be> and
# An Ghysels <An.Ghysels@UGent.be>
#
# This file is part of TAMkin.
#
# TAMkin is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
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


from molmod.molecules import Molecule as BaseMolecule
from molmod.data.periodic import periodic
from molmod.graphs import cached

import numpy


__all__ = ["Molecule"]


class Molecule(BaseMolecule):
    def __init__(self, numbers, coordinates, masses, energy, gradient, hessian, multiplicity):
        BaseMolecule.__init__(self, numbers, coordinates)
        self._masses = numpy.array(masses, float)
        self._masses.setflags(write=False)
        self._energy = energy
        self._gradient = numpy.array(gradient, float)
        self._gradient.setflags(write=False)
        self._hessian = numpy.array(hessian, float)
        self._hessian.setflags(write=False)
        self._multiplicity = multiplicity

    masses = property(lambda self: self._masses)
    energy = property(lambda self: self._energy)
    gradient = property(lambda self: self._gradient)
    hessian = property(lambda self: self._hessian)
    multiplicity = property(lambda self: self._multiplicity)

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

