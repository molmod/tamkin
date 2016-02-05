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
"""Data containers for input from QM or MM simulation codes

   Most IO routines in :mod:`tamkin.io` return instances of the classes
   defined here. These objects are just read-only containers for the QM or MM
   output with a standardized interface.
"""

from tamkin.geom import transrot_basis

from molmod import Molecule as BaseMolecule, MolecularGraph, ReadOnly, \
    UnitCell, ReadOnlyAttribute
from molmod.periodic import periodic
from molmod.graphs import cached

import numpy as np


__all__ = ["Molecule", "BareNucleus", "Proton", "RotScan",
           "translate_pbc"]


class Molecule(BaseMolecule):
    """A container for a Hessian computation output from QM or MM codes."""

    def check_gradient(self, gradient):
        if len(gradient) != self.size:
            raise TypeError("The size of the gradient does not match the "
                "length of the atomic numbers array.")

    def check_hessian(self, hessian):
        if hessian.shape != (self.size*3, self.size*3):
            raise TypeError("The Hessian must be a 3N by 3N matrix where N is "
                "the number of atoms.")

    energy = ReadOnlyAttribute(float, none=False)
    gradient = ReadOnlyAttribute(np.ndarray, none=False,
        check=check_gradient, npdim=2, npshape=(None, 3), npdtype=float)
    hessian = ReadOnlyAttribute(np.ndarray, none=False,
        check=check_hessian, npdim=2, npdtype=float)
    multiplicity = ReadOnlyAttribute(int)
    symmetry_number = ReadOnlyAttribute(int)
    periodic = ReadOnlyAttribute(bool)
    fixed = ReadOnlyAttribute(np.ndarray, npdim=1, npdtype=int)

    def __init__(self, numbers, coordinates, masses, energy, gradient, hessian, multiplicity=None, symmetry_number=None, periodic=False, title=None, graph=None, symbols=None, unit_cell=None, fixed=None):
        """
           Arguments:
            | ``numbers`` -- The atom numbers (integer numpy array with shape N)
            | ``coordinates`` -- the atom coordinates in Bohr (float numpy array
                                 with shape Nx3)
            | ``masses`` -- The atomic masses in atomic units (float numpy array
                            with shape N)
            | ``energy`` -- The molecular energy in Hartree
            | ``gradient`` -- The gradient of the energy, i.e. the derivatives
                              towards Cartesian coordinates, in atomic units
                              (float numpy array with shape Nx3)
            | ``hessian`` -- The hessian of the energy, i.e. the matrix with
                             second order derivatives towards Cartesian
                             coordinates, in atomic units (float numpy array
                             with shape 3Nx3N)
            | ``multiplicity`` -- The spin multiplicity of the electronic system

           Optional arguments:
            | ``symmetry_number`` -- The rotational symmetry number, None when
                                     not known or computed [default=None]
            | ``periodic`` -- True when the system is periodic in three
                              dimensions [default=False]
            | ``title`` -- The title of the system
            | ``graph`` -- The molecular graph of the system
            | ``symbols`` -- A list with atom symbols
            | ``unit_cell`` -- The unit cell vectors for periodic structures
            | ``fixed`` -- An array with indices of fixed atoms
        """
        BaseMolecule.__init__(self, numbers, coordinates, title, masses, graph, symbols, unit_cell)
        self.energy = energy
        self.gradient = gradient
        self.hessian = hessian
        self.multiplicity = multiplicity
        self.symmetry_number = symmetry_number
        self.periodic = periodic
        self.fixed = fixed

    def get_external_basis_new(self, im_threshold=1.0):
        """Create a robust basis for small displacements in the external degrees of freedom.

           The basis is expressed in mass-weighted Cartesian coordinates.
           The three translations are along the x, y, and z-axis. The three
           rotations are about an axis parallel to the x, y, and z-axis
           through the center of mass.

           The returned result depends on the periodicity of the system and on
           the parameter im_threshold:

           * When the system is periodic, only the translation external degrees
             are included. The result is an array with shape (3,3N)

           * When the system is not periodic, the rotational external
             degrees are also included. The result is an array with shape (5,3N)
             or (6,3N). The first three rows correspond to translation, the
             latter rows correspond to rotation. Rotational modes that have a
             angular moment below ``im_threshold`` are discarded.
        """
        if self.periodic:
            result = np.zeros((3, 3*self.size), float)
        else:
            ims, iaxs = np.linalg.eigh(self.inertia_tensor)
            mask = ims > im_threshold
            ims = ims[mask]
            iaxs = iaxs[:,mask]
            ext_dof = len(ims) + 3
            result = np.zeros((ext_dof, 3*self.size), float)
        # translation
        result[0, 0::3] = 1
        result[1, 1::3] = 1
        result[2, 2::3] = 1
        if not self.periodic:
            # rotation
            counter = 3
            for iax in iaxs.transpose():
                result[counter,0::3] = (iax[1]*(self.coordinates[:,2] - self.com[2])
                                       -iax[2]*(self.coordinates[:,1] - self.com[1])) # x
                result[counter,1::3] = (iax[2]*(self.coordinates[:,0] - self.com[0])
                                       -iax[0]*(self.coordinates[:,2] - self.com[2])) # y
                result[counter,2::3] = (iax[0]*(self.coordinates[:,1] - self.com[1])
                                       -iax[1]*(self.coordinates[:,0] - self.com[0])) # z
                counter += 1
        return result

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
           * When the system is not periodic, the rotational external
             degrees are also included. The result is an array with shape (6,3N). The
             first three rows correspond to translation, the latter three rows
             correspond to rotation.
        """
        center = (self.coordinates*self.masses3.reshape((-1,3))).sum(0)/self.mass  # center of mass
        result = transrot_basis(self.coordinates - center, not self.periodic)
        result *= np.sqrt(self.masses3) # transform basis to mass weighted coordinates

        return result

    @cached
    def masses3(self):
        """An array with the diagonal of the mass matrix in Cartesian coordinates.

           Each atom mass is repeated three times. The total length of the
           array is 3N.
        """
        return np.array([self.masses, self.masses, self.masses]).transpose().ravel()

    def get_submolecule(self, selected, energy=None, multiplicity=None, symmetry_number=None, periodic=None, graph=None, title=None, symbols=None, unit_cell=None):
        """Create a submolecule with a selection of atoms

           Argument:
            | ``selected`` -- list of atom indices, numbering starts at 0

           Optional arguments:
            | ``energy`` -- Molecular electronic energy
            | ``multiplicity`` -- The spin multiplicity of the electronic system
            | ``symmetry_number`` -- The rotational symmetry number. Inherited,
                                     None when not known or computed [default=0]
            | ``periodic`` -- True when the system is periodic in three dimensions
                              [default=False]
            | ``title`` -- The title of the system
            | ``graph`` -- The molecular graph of the system
            | ``symbols`` -- A list with atom symbols
            | ``unit_cell`` -- The unit cell vectors for periodic structures

           The function returns a Molecule object consisting of the
           atoms in the selected list. The numbers, coordinates, masses, gradient,
           and Hessian are reduced in size.
           The energy, multiplicity, symmetry_number, periodic and symbols
           attributes are copies of the original molecule attributes, except if
           they are explicitly specified as optional arguments, then they are
           overwritten.
        """
        selected = np.array(selected)
        selected3 = np.array( sum( [[3*at, 3*at+1, 3*at+2] for at in selected] ,[]) )

        # if the following are none, then use the attributes of the original molecule
        if energy is None: energy = self.energy
        if multiplicity is None: multiplicity = self.multiplicity
        if symmetry_number is None: symmetry_number = self.symmetry_number
        if periodic is None: periodic = self.periodic
        if symbols is None and self.symbols is not None:
            symbols = [self.symbols[at] for at in selected]
        if unit_cell is None:
            unit_cell = self.unit_cell

        return Molecule(
            self.numbers[selected],
            self.coordinates[selected,:],
            self.masses[selected],
            energy,
            self.gradient[selected,:],
            self.hessian[selected3,:][:,selected3],
            multiplicity,
            symmetry_number = symmetry_number,
            periodic = periodic,
            title = title,
            graph = graph,
            symbols = symbols,
            unit_cell = unit_cell,
        )

    def write_to_file(self, filename):
        """Write the molecule to a human-readable checkpoint file.

           Argument:
            | ``filename`` -- the file to write to
        """
        from tamkin.io.internal import dump_chk
        data = {}
        for key in "numbers", "coordinates", "masses", "energy", "gradient", \
                   "hessian", "multiplicity", "symmetry_number", "periodic", \
                   "title", "symbols":
            value = getattr(self, key, None)
            if value is not None:
                data[key] = value
        if self.graph is not None:
            data["edges"] = np.array([tuple(edge) for edge in self.graph.edges])
        if self.unit_cell is not None:
            data["cell_vectors"] = self.unit_cell.matrix
            data["cell_active"] = self.unit_cell.active
        dump_chk(filename, data)

    @classmethod
    def read_from_file(cls, filename):
        """Construct a Molecule object from a previously saved checkpoint file

           Arguments:
            | ``filename`` -- the file to load from

           Usage::

             >>> mol = Molecule.read_from_file("mol.chk")

        """
        from tamkin.io.internal import load_chk
        # load the file
        data = load_chk(filename)
        # check the names of the fields:
        mandatory_fields = set([
            "numbers", "coordinates", "masses", "energy", "gradient", "hessian"
        ])
        if not set(data.iterkeys()).issuperset(mandatory_fields):
            raise IOError("The Checkpoint file does not contain the mandatory fields.")
        # take the mandatory fields
        constructor_args = {}
        for mfield in mandatory_fields:
            constructor_args[mfield] = data[mfield]
        # take the optional arguments if present
        opt_fields = ["multiplicity", "symmetry_number", "periodic", "title", "symbols"]
        for ofield in opt_fields:
            if ofield in data:
                constructor_args[ofield] = data[ofield]
        # take the special optional arguments that need conversion
        if "edges" in data:
            graph = MolecularGraph(data["edges"], data["numbers"])
            constructor_args["graph"] = graph
        if "cell_vectors" in data:
            unit_cell = UnitCell(data["cell_vectors"], data.get("cell_active"))
            constructor_args["unit_cell"] = unit_cell
        # construct the molecule object
        return Molecule(**constructor_args)

    def raise_ext(self, shift=1e0):
        """Raise the eigenvalues of the global translations and rotations
        to a high value, such that their coupling with the internal vibrations
        becomes negligible, and they can easily be isolated from the vibrations."""

        # Construct basis for global translations and rotations
        D = self.external_basis.transpose() # mass-weighted
        # make it orthonormal
        svd_threshold = 1e-5
        U, W, Vt = np.linalg.svd(D, full_matrices=False)
        rank = (abs(W) > abs(W[0])*svd_threshold).sum()
        D = U[:,:rank]
        proj = np.dot(D,D.transpose())

        proj1 = proj * (self.masses3.reshape((1,-1)))**(0.5)
        proj2 = proj1* (self.masses3.reshape((-1,1)))**(0.5)

        # Add the shift to the Hessian. The gradient is not changed I guess TODO check this.
        hessian = self.hessian + shift*proj2

        # Use the attributes of the original molecule if they exist
        if hasattr(self,"title"): # check if attribute exists
            title = self.title
        else: title = None
        if hasattr(self,"graph"): # check if attribute exists
            graph = self.graph
        else: graph = None
        if hasattr(self,"symbols"): # check if attribute exists
            symbols = self.symbols
        else: symbols = None
        if hasattr(self,"unit_cell"): # check if attribute exists
            unit_cell = self.unit_cell
        else: unit_cell = None

        return Molecule(
            self.numbers,
            self.coordinates,
            self.masses,
            self.energy,
            self.gradient,
            hessian,             # replace Hessian by new one
            self.multiplicity,
            symmetry_number = self.symmetry_number,
            periodic = self.periodic,
            title = title,
            graph = graph,
            symbols = symbols,
            unit_cell = unit_cell,
        )

    def constrain_ext(self):
        """Project the global translational and rotational vectors
        out of the Hessian and the gradient and return a new Molecule instance."""

        # Construct projector
        D = self.external_basis.transpose() # mass-weighted
        # make it orthonormal
        svd_threshold = 1e-5
        U, W, Vt = np.linalg.svd(D, full_matrices=False)
        rank = (abs(W) > abs(W[0])*svd_threshold).sum()
        D = U[:,:rank]
        proj = np.dot(D,D.transpose())
        #print np.sum((proj-np.dot(proj,proj))**2)

        proj1 = proj * (self.masses3.reshape((1,-1)))**(0.5)
        proj2 = proj1* (self.masses3.reshape((-1,1)))**(-0.5)
        projL = np.identity(len(D)) - proj2.transpose()
        projR = np.identity(len(D)) - proj2
        # test
        #print np.sum((projR-np.dot(projR,projR))**2)
        #print np.sum((projL-np.dot(projL,projL))**2)

        # Project hessian and gradient
        hessian = np.dot(np.dot(projL,self.hessian),projR)
        gradient = np.dot(projL, self.gradient.reshape((-1,1))).reshape((-1,3))

        # Use the attributes of the original molecule if they exist
        if hasattr(self,"title"): # check if attribute exists
            title = self.title
        else: title = None
        if hasattr(self,"graph"): # check if attribute exists
            graph = self.graph
        else: graph = None
        if hasattr(self,"symbols"): # check if attribute exists
            symbols = self.symbols
        else: symbols = None
        if hasattr(self,"unit_cell"): # check if attribute exists
            unit_cell = self.unit_cell
        else: unit_cell = None

        return Molecule(
            self.numbers,
            self.coordinates,
            self.masses,
            self.energy,
            gradient,       # replace gradient by new one
            hessian,             # replace Hessian by new one
            self.multiplicity,
            symmetry_number = self.symmetry_number,
            periodic = self.periodic,
            title = title,
            graph = graph,
            symbols = symbols,
            unit_cell = unit_cell,
        )


class BareNucleus(Molecule):
    """A molecule object for bare nuclei"""

    def __init__(self, number, mass=None):
        """
           Argument:
            | ``number`` -- The atom number

           Optional argument:
            | ``mass`` -- The mass of the atom in atomic units
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
            | ``mass`` -- The mass of the proton in atomic units
        """
        BareNucleus.__init__(self, 1, mass)


class RotScan(ReadOnly):
    """A container for rotational scan data"""

    def check_top_indexes(self, top_indexes):
        if len(top_indexes) < 1:
            raise TypeError("A rotational scan must have at least one atom "
                "rotating")

    dihedral = ReadOnlyAttribute(np.ndarray, none=False, npdim=1,
        npshape=(4,), npdtype=int)
    top_indexes = ReadOnlyAttribute(np.ndarray, none=False,
        check=check_top_indexes, npdim=1, npdtype=int)
    potential = ReadOnlyAttribute(np.ndarray, npdim=2, npshape=(2, None),
        npdtype=float)

    def __init__(self, dihedral, molecule=None, top_indexes=None, potential=None):
        """
           Arguments
            | ``dihedral`` -- the index of the atoms that define the dihedral
                              angle

           Optional arguments
            | ``molecule`` -- a molecule object. required when top_indexes is not
                              given
            | ``top_indexes`` -- a list of atom indexes involved in the rotor.
                                 required when molecule is not given or when the
                                 top_indexes can not be derived automatically.
            | ``potential`` -- rotational potential info (if this is a hindered
                               rotor). must be a two-tuple containing the angles
                               and the corresponding energies.

        """
        self.dihedral = dihedral
        if top_indexes is None:
            # try to deduce the top indexes
            if molecule is None:
                raise ValueError("missing arguemt: top_indexes or molecule from which top_indexes can be derived")
            atom0, atom1 = dihedral[1:3]
            graph = MolecularGraph.from_geometry(molecule)
            half0, half1 = graph.get_halfs(atom0, atom1)
            if (len(half0) + len(half1) != molecule.size) or len(half0) == 1 or len(half1) == 1:
                raise ValueError("The rotating top could not be assigned properly. Specify the top_indexes manually.")
            if len(half1) > len(half0):
                top_indexes = half0
                top_indexes.discard(atom0)
            else:
                top_indexes = half1
                top_indexes.discard(atom1)
            top_indexes = tuple(top_indexes)
        else:
            # check the sanity of the top indexes
            if not ((self.dihedral[0] in top_indexes) ^ (self.dihedral[3] in top_indexes)):
                raise ValueError('The top must contain either first or the last atom of the dihedral angle.')
            if ((self.dihedral[1] in top_indexes) | (self.dihedral[2] in top_indexes)):
                raise ValueError('The top may not contain atoms from the central bond of the dihedral angle.')
        self.top_indexes = top_indexes
        self.potential = potential


def translate_pbc(molecule, selected, displ, vectors = None):
    """Translate the structure along the lattice vectors

       This method is meant to be used in periodic structures, where periodic
       boundary conditions apply (pbc).

       Arguments:
        | ``molecule`` -- A Molecule instance.
        | ``selected`` -- A list of indices of the atoms that will be displaced.
        | ``displ`` -- A list of 3 integers: [i0,i1,i2]. The selected atoms will
                       be displaced over i0 lattice distances in the 0-axis
                       direction, similarly for i1 and i2.

       Optional argument:
        | ``vectors`` -- The lattice vectors, one in each column. If not
                         specified, the vectors in the unit_cell attribute of
                         the molecule is used.
    """
    coordinates = molecule.coordinates.copy()
    if vectors is None:
        if molecule.unit_cell is None:
            raise ValueError("No vectors for axes known")
        else:
            vectors = molecule.unit_cell.matrix

    # translate
    for atom in selected:
        for axis in range(3):
            coordinates[atom,:] += displ[axis]*vectors[:,axis]

    return molecule.copy_with(coordinates=coordinates)
