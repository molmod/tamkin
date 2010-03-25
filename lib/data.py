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
"""Data containers for input from QM or MM simulation codes

   Most IO routines in :mod:`tamkin.io` return instances of the classes
   defined here. These objects are just read-only containers for the QM or MM
   output with a standardized interface.
"""

from tamkin.geom import transrot_basis

from molmod import Molecule as BaseMolecule, MolecularGraph, ReadOnly
from molmod.periodic import periodic
from molmod.graphs import cached

import numpy


__all__ = ["Molecule", "BareNucleus", "Proton", "RotScan"]


class Molecule(BaseMolecule):
    """A container for a Hessian computation output from QM or MM codes."""

    def __init__(self, numbers, coordinates, masses, energy, gradient, hessian, multiplicity, symmetry_number=0, periodic=False, title=None, graph=None, symbols=None):
        """
           Arguments:
            | numbers  --  The atom numbers (integer numpy array with shape N)
            | coordinates  --  the atom coordinates in Bohr (float numpy array
                               with shape Nx3)
            | masses  --  The atomic masses in atomic units (float numpy array
                          with shape N)
            | energy  --  The molecular energy in Hartree
            | gradient  --  The gradient of the energy, i.e. the derivatives
                            towards Cartesian coordinates, in atomic units
                            (float numpy array with shape Nx3)
            | hessian  --  The hessian of the energy, i.e. the matrix with
                           second order derivatives towards Cartesian
                           coordinates, in atomic units (float numpy array with
                           shape 3Nx3N)
            | multiplicity  --  The spin multiplicity of the electronic system

           Optional arguments:
            | symmetry_number  --  The rotational symmetry number, 0 when not
                                   known or computed [default=0]
            | periodic  --  True when the system is periodic in three dimensions
                            [default=False]
            | title  --  The title of the system
            | graph  --  The molecular graph of the system
            | symbols  --  A list with atom symbols
        """
        ReadOnly.__init__(self)
        mandatory = {
            "numbers": numpy.array(numbers, int),
            "coordinates": numpy.array(coordinates, float),
            "masses": numpy.array(masses, float),
            "energy": energy,
            "gradient": numpy.array(gradient, float),
            "hessian": numpy.array(hessian, float),
            "multiplicity": multiplicity,
            "symmetry_number": symmetry_number,
            "periodic": periodic,
        }
        optional = {
            "title": title,
            "graph": graph,
            "symbols": symbols,
        }
        self._init_attributes(mandatory, optional)

    @cached
    def external_basis(self):
        """The basis for small displacements in the external degrees of freedom.

           The basis is expressed in mass-weighted Cartesian coordinates.
           The three translations are along the x, y, and z-axis. The three
           rotations are about an axis parallel to the x, y, and z-axis
           through the center of mass.
           The returned result depends on the periodicity of the system:

           * When the system is periodic, only the translation external degrees
             are included. The result is an array with shape (3,3N)
           * When the system is not periodic, only the translation external
             degrees are included. The result is an array with shape (6,3N). The
             first three rows correspond to translation, the latter three rows
             correspond to rotation.
        """
        center = (self.coordinates*self.masses3.reshape((-1,3))).sum(0)/self.mass  # center of mass
        result = transrot_basis(self.coordinates - center, not self.periodic)
        result *= numpy.sqrt(self.masses3) # transform basis to mass weighted coordinates

        return result

    @cached
    def masses3(self):
        """An array with the diagonal of the mass matrix in Cartesian coordinates.

           Each atom mass is repeated three times. The total length of the
           array is 3N.
        """
        return numpy.array([self.masses, self.masses, self.masses]).transpose().ravel()


class BareNucleus(Molecule):
    """A molecule object for bare nuclei"""

    def __init__(self, number, mass=None):
        """
           Argument:
            | number  --  The atom number

           Optional argument:
            | mass  --  The mass of the atom in atomic units
        """
        if mass is None:
            mass = periodic[number].mass
        Molecule.__init__(
            self, [number], [[0,0,0]], [mass], 0.0, [[0,0,0]],
            [[0,0,0],[0,0,0],[0,0,0]], 1, 0, False
        )


class Proton(BareNucleus):
    """A molecule object for a proton (or one of its isotopes)"""

    def __init__(self, mass=None):
        """
           Optional argument:
            | mass  --  The mass of the proton in atomic units
        """
        BareNucleus.__init__(self, 1, mass)


class RotScan(object):
    """A container for rotational scan data"""

    def __init__(self, dihedral, molecule=None, top_indexes=None, potential=None):
        """
           Arguments
            | dihedral  --  the index of the atoms that define the dihedral
                            angle

           Optional arguments
            | molecule  --  a molecule object. required when top_indexes is not
                            given
            | top_indexes  --  a list of atom indexes involved in the rotor.
                               required when molecule is not given
            | potential  --  rotational potential info (if this is a hindered
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


