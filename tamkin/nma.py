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
"""Normal mode analysis with default and extended schemes.

A normal mode analysis is carried out by constructing an NMA object. The first
argument is a molecule object created by one of the IO routines in
:mod:`tamkin.io`. ::

>>> nma = NMA(molecule)

This leads to a standard normal mode analysis in 3*N degrees of freedom.
The results, including those relevant for the construction of the molecular
partition function, are stored as attributes of the NMA object. For example::

>>> print nma.freqs

prints the frequencies of the normal modes. Note that all data is stored in
atomic units and that the freqs array contains really frequencies, not
wavenumbers. If you want to print the wavenumbers in cm**-1, use the unit
conversion constants from the ``molmod`` package::

>>> from molmod import centimeter, lightspeed
>>> invcm = lightspeed/centimeter
>>> print nma.freqs/invcm

One can also use modified schemes by giving a second argument to the NMA
constructor. The following example computes the normal modes in 3*N-6 degrees
of freedom::

>>> nma = NMA(molecule, ConstrainExt())

The second argument is an instance of a class that derives from the
:class:`Treatment` class. Other treatments include: :class:`Full` (the default),
:class:`PHVA`, :class:`VSA`, :class:`VSANoMass`, :class:`MBH`,
:class:`PHVA_MBH`, :class:`Constrain`, and :class:`MBHConstrainExt`.
"""

# A few conventions for the variables names:
#
# * The suffix _small refers to the small set of degrees of freedom that
#   obey constraints implicitly. The absence of the suffix refers to
#   conventional Cartesian coordinates, e.g. hessian is the 3Nx3N Hessian
#   and hessian_small is the Hessian in the new coordinates.
#
# * The suffix _mw stands for mass-weighted, e.g. hessian_mw is the mass
#   weighted Hessian in Cartesian coordinates. hessiam_small_mw is the mass
#   weighted Hessian in the new coordinates.


from tamkin.data import Molecule
from tamkin.geom import transrot_basis, rank_linearity
from tamkin.io.internal import load_chk, dump_chk

import numpy as np


__all__ = [
    "NMA", "AtomDivision", "Transform", "MassMatrix", "Treatment",
    "Full", "ConstrainExt", "PHVA", "VSA", "VSANoMass", "MBH",
    "Blocks","PHVA_MBH", "Constrain", "MBHConstrainExt",
]


class NMA(object):
    """A generic normal mode analysis class.

       This class gathers the functionality that is common between all types of
       NMA variations, i.e. computation of frequencies and modes, once the
       problem is transformed to reduced coordinates. The actual nature of the
       reduced coordinates is determined by the treatment argument.
    """

    def __init__(self, molecule, treatment=None, do_modes=True):
        """
           Arguments:
            | ``molecule`` -- a molecule object obtained from a routine in
                              :mod:`tamkin.io`

           Optional arguments:
            | ``treatment`` -- an instance of a Treatment subclass
                               [default=Full()]
            | ``do_modes`` -- When False, only the frequencies are computed.
                              When True, also the normal modes are computed.
                              [default=True]

           Referenced attributes of molecule:
              ``mass``, ``masses``, ``masses3``, ``numbers``, ``coordinates``,
              ``inertia_tensor``, ``multiplicity``, ``symmetry_number``,
              ``periodic``, ``energy``

           Extra attributes:
            | ``freqs`` -- array of frequencies
            | ``modes`` -- array of mass-weighted Cartesian modes (if do_modes
                           is True). Each column corresponds to one mode. One
                           has to divide a column the square root of the masses3
                           attribute to obtain the mode in non-mass-weighted
                           coordinates.
            | ``zeros`` -- list of indices of zero frequencies

        """
        if treatment == None:
            treatment = Full()

        # the treatment object will store the results as attributes
        treatment(molecule, do_modes)
        # treatment.hessian_small:
        #    the Hessian in reduced coordinates
        # treatment.mass_matrix_small:
        #    the mass matrix in reduced coordinates (see MassMatrix class)
        # treatment.transform: (None if do_modes=False)
        #    the transformation from small displacements in reduced coordinates
        #    to small displacements in Cartesian coordinates. (see Transform class)
        # treatment.num_zeros:
        #    the number of zero eigenvalues to expect
        # treatment.external_basis: (None if do_modes=False)
        #    the basis of external degrees of freedom. number of basis vectors
        #    matches the number of zeros.
        #
        # For the implementation of certain treatments, it is easier to produce
        # a mass-weighted small Hessian immediately. In such cases, the
        # transform is also readily mass-weighted and mass_matrix_small is
        # None.

        # the conventional frequency computation in the reduced coordinates
        if treatment.mass_matrix_small is None:
            hessian_small_mw = treatment.hessian_small
        else:
            hessian_small_mw = treatment.mass_matrix_small.get_weighted_hessian(treatment.hessian_small)
        del treatment.hessian_small # save memory

        if hessian_small_mw.size == 0:
            self.freqs = np.array([])
            self.modes = np.array([])
            self.zeros = []
        else:
            if do_modes:
                evals, modes_small_mw = np.linalg.eigh(hessian_small_mw)
            else:
                evals = np.linalg.eigvalsh(hessian_small_mw)
                modes_small_mw = None

            # frequencies
            self.freqs = np.sqrt(abs(evals))/(2*np.pi)
            # turn imaginary frequencies into negative frequencies
            self.freqs *= (evals > 0)*2-1

            if do_modes:
                # At this point the transform object transforms unweighted reduced
                # coordinates into Cartesian coordinates. Now we will alter it, so that
                # it transforms from weighted reduced coordinates to Cartesian
                # coordinates.
                if treatment.mass_matrix_small is not None:
                    treatment.transform.make_weighted(treatment.mass_matrix_small)
                # transform the modes to unweighted Cartesian coordinates.
                self.modes = treatment.transform(modes_small_mw)
                # transform the modes to weighted Cartesian coordinates.
                self.modes *= molecule.masses3.reshape((-1,1))**0.5
            else:
                self.modes = None

            # guess which modes correspond to the zero frequencies
            if treatment.num_zeros == 0:
                # don't bother
                self.zeros = []
            else:
                if do_modes:
                    # take the 20 lowest modes and compute the overlap with the
                    # external basis
                    num_try = 20
                    to_try = abs(self.freqs).argsort()[:num_try]   #indices of lowest 20 modes
                    overlaps = np.zeros(num_try, float)
                    for counter, i in enumerate(to_try):
                        components = np.dot(treatment.external_basis, self.modes[:,i])
                        overlaps[counter] = np.linalg.norm(components)
                    self.zeros = to_try[overlaps.argsort()[-treatment.num_zeros:]]
                else:
                    self.zeros = abs(self.freqs).argsort()[:treatment.num_zeros]

        # a few more attributes that are worth keeping
        self.mass = molecule.mass
        self.masses = molecule.masses
        self.masses3 = molecule.masses3
        self.numbers = molecule.numbers
        self.coordinates = molecule.coordinates
        self.inertia_tensor = molecule.inertia_tensor
        self.multiplicity = molecule.multiplicity
        self.symmetry_number = molecule.symmetry_number
        self.periodic = molecule.periodic
        self.energy = molecule.energy
        self.title = molecule.title
        self.chemical_formula = molecule.chemical_formula

    def write_to_file(self, filename, fields='all'):
        """Write the NMA results to a human-readable checkpoint file.

           Argument:
            | ``filename`` -- the file to write to

           Optional argument:
            | ``fields`` -- define the selection of attributes to be written to
                            file. This is one of 'all' (all attributes), 'modes'
                            (only attributes required for nmatools.py), or
                            'partf' (only attributes required for the
                            construction of a partition function)
        """
        if fields == 'all':
            data = dict((key, val) for key, val in self.__dict__.iteritems())
        elif fields == 'modes':
            keys = ["freqs", "modes", "masses", "numbers", "coordinates", "zeros", "title"]
            data = dict((key, self.__dict__[key]) for key in keys)
        elif fields == 'partf':
            keys = [
                "freqs", "mass", "masses3", "inertia_tensor", "multiplicity",
                "symmetry_number", "periodic", "energy", "zeros", "title",
                "chemical_formula",
            ]
            data = dict((key, self.__dict__[key]) for key in keys)
        dump_chk(filename, data)

    @classmethod
    def read_from_file(cls, filename):
        """Construct an NMA object from a previously saved checkpoint file

           Arguments:
            | ``filename`` -- the file to load from

           Usage::

             >>> nma = NMA.read_from_file("foo.chk")

        """
        # ugly way to bypass the default constructor
        result = cls.__new__(cls)
        # load the file
        data = load_chk(filename)
        # check the names of the fields:
        possible_fields = set([
            "freqs", "modes", "mass", "masses", "masses3", "numbers",
            "coordinates", "inertia_tensor", "multiplicity", "symmetry_number",
            "periodic", "energy", "zeros", "title", "chemical_formula",
        ])
        if not set(data.iterkeys()).issubset(possible_fields):
            raise IOError("The Checkpoint file does not contain the correct fields.")
        # assign the attributes
        result.__dict__.update(data)
        return result


class AtomDivision(object):
    """A division of atoms into transformed, free and fixed."""

    def __init__(self, transformed, free, fixed):
        """
           Arguments:
            | ``transformed`` -- the atom indices of the atoms whose coordinates
                                 are transformed into non-Cartesian coordinates.
            | ``free`` -- the atom indices that are not transformed and retained
                          as Cartesian coordinates in the new set of coordinates
            | ``fixed`` -- the atoms that are not used for the new coordinates,
                           i.e. their positions are constrained.
        """
        self.transformed = np.array(transformed, int)
        self.free = np.array(free, int)
        self.fixed = np.array(fixed, int)

        self.num_cartesian = 3*(len(self.transformed)+len(self.free)+len(self.fixed))
        self.to_cartesian_order = np.zeros(self.num_cartesian, int)
        self.to_reduced_order = np.zeros(self.num_cartesian, int)
        counter = 0
        for l in self.transformed, self.free, self.fixed:
            for i in l:
                # index corresponds to Cartesian index
                # value corresponds to reduced index
                self.to_cartesian_order[3*i] = counter
                self.to_cartesian_order[3*i+1] = counter+1
                self.to_cartesian_order[3*i+2] = counter+2
                self.to_reduced_order[counter] = 3*i
                self.to_reduced_order[counter+1] = 3*i+1
                self.to_reduced_order[counter+2] = 3*i+2
                counter += 3


class Transform(object):
    """A clever transformation object. It is sparse when atom coordinates remain
       Cartesian in the reduced coordinates.

       This object transforms small displacements (first order) in reduced
       internal coordinates (can be mass weighted) into plain Cartesian
       coordinates.

       It is assumed that the reduced coordinates are always split into two
       parts (in order):

       1) the coordinates that are non-Cartesian
       2) the free coordinates that are Cartesian

    """

    def __init__(self, matrix, atom_division=None):
        """
           Arguments:
             | ``matrix`` -- the linear transformation from the transformed
                             displacements to Cartesian coordinates.

           Optional argument
             | ``atom_division`` -- an AtomDivision instance, when not given all
                                    atom coordinates are `transformed`

           Attributes:
             | ``matrix`` -- see above
             | ``scalars`` -- diagonal part of the linear transformation (only
                              used with mass-weighted transformations)
        """
        if matrix is None:
            matrix = np.zeros((0,0), float)
        if atom_division is None:
            # internal usage only:
            self._num_reduced = matrix.shape[1]
        else:
            # Quality Assurance:
            if matrix.shape[0] != 3*len(atom_division.transformed):
                raise ValueError("The matrix must have %i columns (matching the number of transformed atoms), got %i." %
                    3*len(atom_division.transformed), matrix.shape[0]
                )
            # internal usage only:
            self._num_reduced = matrix.shape[1] + 3*len(atom_division.free)

        self.matrix = matrix
        self.atom_division = atom_division
        # as long as we do not work with mass weighted coordinates, the
        # following remains None. In case of weighted coordinates, this
        # becomes a scaling vector with 3*len(free) floats:
        self.scalars = None
        self._weighted = False

    def get_weighted(self):
        """Return True when the transform is already mass-weighted"""
        return self._weighted

    weighted = property(get_weighted)

    def __call__(self, modes):
        """Transform small displacement vectors from new to Cartesian coordinates.

           Argument:
            | ``modes`` -- Small (mass-weighted) displacements (or modes) in
                           internal coordinates (float numpy array with shape
                           KxM, where K is the number of internal coordinates
                           and M is the number of modes)

           Returns:
              Small non-mass-weighted displacements (or modes) in Cartesian
              coordinates (float numpy array with shape 3NxM, where N is the
              number of Cartesian coordinates and M is the number of modes)

           Usage::

             >>> transform = Transform(...)
             >>> modes_cartesian = transform(modes_internal)

        """
        # Quality Assurance:
        if len(modes.shape) != 2:
            raise ValueError("Modes must be a two-dimensional array.")
        if modes.shape[0] != self._num_reduced:
            raise ValueError("The modes argument must be an array with %i rows, got %i." %
                (self._num_reduced, modes.shape[0])
            )
        # Computation
        if self.atom_division is None:
            return np.dot(self.matrix, modes)
        else:
            result = np.zeros((self.atom_division.num_cartesian, modes.shape[1]), float)  # 3NxM
            i1 = 3*len(self.atom_division.transformed)
            i2 = i1 + 3*len(self.atom_division.free)
            result[:i1] = np.dot(self.matrix, modes[:self.matrix.shape[1]])
            if self.weighted:
                result[i1:i2] = modes[self.matrix.shape[1]:]*self.scalars
            else:
                result[i1:i2] = modes[self.matrix.shape[1]:]
            #    result[:,i2:] remains zero because these atoms are fixed
            # Reorder the atoms and return the result
            tmp = result[self.atom_division.to_cartesian_order]
            return tmp

    def make_weighted(self, mass_matrix):
        """Include mass-weighting into the transformation.

           The original transformation is from non-mass-weighted new coordinates
           to non-mass-weighted Cartesian coordinates and becomes a transform
           from mass-weighted new coordinates to non-mass-weighted Cartesian
           coordinates.

           Argument:
            | ``mass_matrix`` -- A MassMatrix instance for the new coordinates
        """
        # modifies the transformation matrix in place:
        # the transformation matrix always transforms to non-mass-weighted Cartesian coords
        if self.weighted:
            raise Exception("The transformation is already weighted.")
        self.matrix = np.dot(self.matrix, mass_matrix.mass_block_inv_sqrt)
        self.scalars = mass_matrix.mass_diag_inv_sqrt.reshape((-1,1))
        self._weighted = True


class MassMatrix(object):
    """A clever mass matrix object. It is sparse when atom coordinates remain
       Cartesian in the reduced coordinates.
    """

    def __init__(self, *args):
        """
           Arguments, if one is given and it is a two-dimensional matrix:
            | ``mass_block`` -- the mass matrix associated with the transformed
                                coordinates

           Arguments, if one is given and it is a one-dimensional matrix:
            | ``mass_diag`` -- the diagonal of the mass matrix associated with
                               the free atoms (each mass appears three times)

           Arguments, if two are given:  ! Attention for order of arguments.
            | ``mass_block`` -- the mass matrix associated with the transformed
                                coordinates
            | ``mass_diag`` -- the diagonal of the mass matrix associated with
                               the free atoms (each mass appears three times)

           The mass of the fixed atoms does not really matter here.
        """
        if len(args) == 1:
            if len(args[0].shape) == 1:
                self.mass_diag = args[0]
                self.mass_block = np.zeros((0,0), float)
            elif len(args[0].shape) == 2:
                self.mass_diag = np.zeros((0,), float)
                self.mass_block = args[0]
            else:
                raise TypeError("When MassMatrix.__init__ gets one argument, it must be a one- or two-dimensional array.")
        elif len(args) == 2:
            self.mass_block = args[0]  #mass_block is first argument
            self.mass_diag  = args[1]  #mass_diag is second argument
        else:
            raise TypeError("MassMatrix.__init__ takes one or two arguments, %i given." % len(args))

        # the square root of the inverse
        if len(self.mass_block) == 0:
            self.mass_block_inv_sqrt = np.zeros((0,0), float)
        else:
            evals, evecs = np.linalg.eigh(self.mass_block)
            self.mass_block_inv_sqrt = np.dot(evecs/np.sqrt(evals), evecs.transpose())
        self.mass_diag_inv_sqrt = 1/np.sqrt(self.mass_diag)

    def get_weighted_hessian(self, hessian):
        hessian_mw = np.zeros(hessian.shape,float)
        n = len(self.mass_block)
        # transform block by block:
        hessian_mw[:n,:n] = np.dot(np.dot(self.mass_block_inv_sqrt, hessian[:n,:n]), self.mass_block_inv_sqrt)
        hessian_mw[:n,n:] = np.dot(self.mass_block_inv_sqrt, hessian[:n,n:])*self.mass_diag_inv_sqrt
        hessian_mw[n:,:n] = hessian[:n,n:].transpose()
        hessian_mw[n:,n:] = (hessian[n:,n:]*self.mass_diag_inv_sqrt).transpose()*self.mass_diag_inv_sqrt
        return hessian_mw


class Treatment(object):
    """An abstract base class for the NMA treatments. Derived classes must
       override the __call__ function, or they have to override the individual
       compute_zeros and compute_hessian methods. Parameters specific for the
       treatment are passed to the constructor, see for example the PHVA
       implementation.
    """

    def __init__(self):
        self.hessian_small = None
        self.mass_matrix_small = None
        self.transform = None
        self.num_zeros = None
        self.external_basis = None

    def __call__(self, molecule, do_modes):
        """Calls compute_hessian and compute_zeros (in order) with same arguments

           Arguments:
            | ``molecule`` -- a Molecule instance
            | ``do_modes`` -- a boolean indicates whether the modes have to be
                              computed
        """
        self.compute_hessian(molecule, do_modes)
        self.compute_zeros(molecule, do_modes)

    def compute_hessian(self, molecule, do_modes):
        """To be computed in derived classes

           Arguments:
            | ``molecule`` -- a Molecule instance
            | ``do_modes`` -- a boolean indicates whether the modes have to be

           Attributes to be computed:

           * ``treatment.hessian_small``: the Hessian in reduced coordinates
           * ``treatment.mass_matrix_small``: the mass matrix in reduced
             coordinates (see MassMatrix class)
           * ``treatment.transform``: (None if ``do_modes==False``) the
             transformation from small displacements in reduced coordinates
             to small displacements in Cartesian coordinates. (see Transform
             class)

           For the implementation of certain treatments, it is easier to produce
           a mass-weighted small Hessian immediately. In such cases, the
           transform is readily mass-weighted and mass_matrix_small is None.
        """
        raise NotImplementedError

    def compute_zeros(self, molecule, do_modes):
        """To be computed in derived classes

           Arguments:
            | ``molecule`` -- a Molecule instance
            | ``do_modes`` -- a boolean indicates whether the modes have to be

           Attributes to be computed:

           * ``treatment.num_zeros``: the number of zero eigenvalues to expect
           * ``treatment.external_basis``: (None if ``do_modes=False``) the
             basis of external degrees of freedom. number of basis vectors
             matches the number of zeros. These basis vectors are mass-weighted.
        """
        raise NotImplementedError


class Full(Treatment):
    """A full vibrational analysis, without transforming to a new set of
       coordinates.
    """
    def __init__(self, im_threshold=1.0):
        """
           Optional argument:
            | ``im_threshold`` -- Threshold for detection of deviations from
                                  linearity. When a moment of inertia is below
                                  this threshold, it is treated as a zero.
        """
        self.im_threshold = im_threshold
        Treatment.__init__(self)

    def compute_zeros(self, molecule, do_modes):
        """See :meth:`Treatment.compute_zeros`.

           The number of zeros should be:

           - 3 for a single atom, nonperiodic calculation
           - 5 for a linear molecule, nonperiodic calculation
           - 6 for a nonlinear molecule, nonperiodic calculation
           - 3 in periodic calculations
        """
        # determine nb of zeros
        external_basis = molecule.get_external_basis_new(self.im_threshold)
        self.num_zeros = external_basis.shape[0]

        # check
        if molecule.periodic:
            assert self.num_zeros == 3, "Number of zeros is expected to be 3 "\
                "(periodic calculation), but found %i." % self.num_zeros
        else:
            assert self.num_zeros in [3,5,6], "Number of zeros is expected to "\
                "be 3, 5 or 6, but found %i." % self.num_zeros

        if do_modes:
            # Mass-weighted and orthonormal basis vectors for external degrees
            # of freedom. These are used to detect which vibrational modes match
            # the external degrees of freedom.
            U, W, Vt = np.linalg.svd(molecule.get_external_basis_new(), full_matrices=False)
            self.external_basis = Vt

    def compute_hessian(self, molecule, do_modes):
        """See :meth:`Treatment.compute_hessian`.

        The Hessian is the full 3Nx3N Hessian matrix ``H``.
        The mass matrix is the full 3Nx3N mass matrix ``M``.
        It is assumed that the coordinates are Cartesian coordinates, so the
        mass matrix is diagonal.
        """
        self.hessian_small = molecule.hessian
        self.mass_matrix_small = MassMatrix(molecule.masses3)
        if do_modes:
            atom_division = AtomDivision([], np.arange(molecule.size), [])
            self.transform = Transform(None, atom_division)


class ConstrainExt(Treatment):
    """Almost a full vibrational analysis, but with constrained external degrees
       of freedom.

       Note that the current implementation only works correctly when the
       gradient is zero.
    """

    def __init__(self, gradient_threshold=1e-4, im_threshold=1.0):
        """
           Optional arguments:
            | ``gradient_threshold`` -- The maximum allowed value of the
                                        components of the Cartesian gradient in
                                        atomic units. When the threshold is
                                        exceeded, a ValueError is raised.
                                        [default=1-e4]
            | ``im_threshold`` -- Threshold for detection of deviations from
                                  linearity. When a moment of inertia is below
                                  this threshold, it is treated as a zero.
        """
        self.gradient_threshold = gradient_threshold
        self.im_threshold = im_threshold
        Treatment.__init__(self)

    def compute_zeros(self, molecule, do_modes):
        """See :meth:`Treatment.compute_zeros`.

           The number of zeros is set to 0, because the global translations
           and rotations are already projected out.
        """
        self.num_zeros = 0

    def compute_hessian(self, molecule, do_modes):
        """See :meth:`Treatment.compute_hessian`

           First a basis is constructed for the internal coordinates. The 3N-6
           (or 3N-5) basis vectors of length 3N (matrix B is (3N-6)x3N) are
           mass-weighted.
           The ConstrainExt Hessian is then: ``B^T H B``.
           This matrix is already mass weigthed, such that no ConstrainExt mass
           matrix needs to be specified.
        """
        # Compute the rmsd of the gradient per atom. The maximum of these rmsds
        # should be below the threshold. This is a rotationally invariant test.
        atom_gradient_norms = np.sqrt((molecule.gradient**2).mean(axis=1))
        gradmax = (atom_gradient_norms).max()
        if gradmax > self.gradient_threshold:
            raise ValueError(
                "The rmsd of the gradient on some atoms exceeds the threshold "
                "(%.1e > %.1e). The current implementation of the ConstrainExt "
                "treatment only works on optimized geometries." % (
                    gradmax, self.gradient_threshold
                )
            )
        # project the hessian on the orthogonal complement of the basis of small
        # displacements in the external degrees of freedom.
        external_basis = molecule.get_external_basis_new(self.im_threshold)
        U, W, Vt = np.linalg.svd(molecule.external_basis, full_matrices=True)
        rank = external_basis.shape[0]
        internal_basis_mw = (Vt[rank:]/np.sqrt(molecule.masses3)).transpose()
        # the following hessian is already mass-weighted;
        self.hessian_small = np.dot(internal_basis_mw.transpose(), np.dot(molecule.hessian, internal_basis_mw))
        # we do not define mass_matrix_small since it is useless when the hessian
        # is already mass-weighted
        if do_modes:
            # also mass-weighted because the hessian is mass-weighted too:
            self.transform = Transform(internal_basis_mw)


class PHVA(Treatment):
    """Perform the partial Hessian vibrational analysis.

       Part of the system is fixed during the vibrational analysis: the fixed
       atoms are kept at their reference positions. The rest of the atoms can
       still vibrate.

       See references:

       - J.D. Head, Int. J. Quantum Chem. 65, 827 (1997)
       - J.D. Head and Y. Shi, Int. J. Quantum Chem. 75, 81 (1999)
       - J.D. Head, Int. J. Quantum Chem. 77, 350 (2000)
       - H. Li and J. Jensen, Theor. Chem. Acc. 107, 211 (2002)

    """
    def __init__(self, fixed, svd_threshold=1e-5):
        """
           Argument:
            | ``fixed`` -- a list with fixed atoms, counting starts from zero.

           Optional argument:
            | ``svd_threshold`` -- threshold for detection of deviations for
                                   linearity
        """
        # QA:
        if len(fixed) == 0:
            raise ValueError("At least one fixed atom is required.")
        # Rest of init:
        self.fixed = np.array(fixed)
        self.fixed.sort()
        self.svd_threshold = svd_threshold
        Treatment.__init__(self)

    def compute_zeros(self, molecule, do_modes):
        """See :meth:`Treatment.compute_zeros`.

        This is a bit tricky. Most of the times the number of zero eigenvalues
        is zero, but there are a few exceptions. When there is one fixed
        point, there are in general three zeros. When there are two (or more
        colinear fixed atoms), there is in general one zero. When both the
        fixed and the free atoms are colinear, there are no zeros.
        """
        # An unambigous way to define the 'external' degrees of freedom is as
        # follows: first construct an external basis of the entire systems,
        # TODO: this will fail if the molecule is displaced far from the origin
        # TODO: make it complicated and analyze the inertia tensor
        # TODO: make ext_dof a molecule property
        U, W, Vt = np.linalg.svd(molecule.external_basis, full_matrices=False)
        rank = (abs(W) > abs(W[0])*self.svd_threshold).sum()
        external_basis = Vt[:rank]
        # then project this basis on a subspace of the fixed atoms and try to
        # find linear combinations that do not move the fixed atoms.
        fixed3 = []
        for i in self.fixed:
            fixed3.append(3*i)
            fixed3.append(3*i+1)
            fixed3.append(3*i+2)
        system = external_basis[:,fixed3].transpose()
        # The homogenuous solutions of the system corresponds to remaining
        # degrees of freedom. The number of homogenuous solutions is equal to
        # the nullity of the system, i.e. the number of zero singular values.
        U, W, Vt = np.linalg.svd(system, full_matrices=False)
        self.num_zeros = (abs(W) < abs(W[0])*self.svd_threshold).sum()
        if do_modes and self.num_zeros > 0:
            # TODO: fix this
            #self.exernal_basis = ...
            raise NotImplementedError

    def compute_hessian(self, molecule, do_modes):
        """See :meth:`Treatment.compute_hessian`.

        The Hessian matrix for the PHVA is the submatrix of the full (3Nx3N)
        Hessian, corresponding with the non-fixed atoms: ``H_nonfixed``.
        The mass matrix for the PHVA is the (diagonal) submatrix of the
        full (3Nx3N) mass matrix, corresponding with the non-fixed atoms:
        ``M_nonfixed``.
        So it is a diagonal matrix with the masses of the non-fixed atoms on
        the diagonal.
        """
        free = np.zeros(molecule.size - len(self.fixed), int)
        free3 = np.zeros(len(free)*3, int)
        counter_fixed = 0
        counter_free = 0
        for i in xrange(molecule.size):
            if counter_fixed < len(self.fixed) and self.fixed[counter_fixed] == i:
                counter_fixed += 1
            else:
                free[counter_free] = i
                free3[counter_free*3] = i*3
                free3[counter_free*3+1] = i*3+1
                free3[counter_free*3+2] = i*3+2
                counter_free += 1

        self.hessian_small = molecule.hessian[free3][:,free3]
        masses3_small = molecule.masses3[free3]
        self.mass_matrix_small = MassMatrix(masses3_small)
        if do_modes:
            atom_division = AtomDivision([], free, self.fixed)
            self.transform = Transform(None, atom_division)


class VSA(Treatment):
    """
    Perform a Vibrational Subsystem Analysis.

    Frequencies and modes are computed with the VSA approach, as described in
    the references:

    - W. Zheng, B.R. Brooks, J. Biophys. 89, 167 (2006)
    - H.L. Woodcock, W. Zheng, A. Ghysels, Y. Shao, J. Kong, B.R. Brooks,
      J. Chem. Phys. 129 (21), Art. No. 214109 (2008)

    The system is partitioned into a subsystem and an environment. The subsystem
    atoms are allowed to vibrate, while the environment atoms follow the motions
    of the subsystem atoms. The environment atoms are force free.
    """
    def __init__(self, subs, svd_threshold=1e-5):
        """
           One argument:
            | ``subs`` -- a list with the subsystem atoms, counting starts from
                          zero.

           Optional argument:
            | ``svd_threshold`` -- threshold for detection of deviations for
                                   linearity
        """
        # QA:
        if len(subs) == 0:
            raise ValueError("At least one subsystem atom is required.")
        # Rest of init:
        self.subs = np.array(subs)
        #self.subs.sort()
        self.svd_threshold = svd_threshold
        Treatment.__init__(self)

    def compute_zeros(self, molecule, do_modes):
        """See :meth:`Treatment.compute_zeros`.

        The number of zeros should be:

        - 3 for subsystem = a single atom, nonperiodic calculation
        - 5 for subsystem = a linear molecule, nonperiodic calculation
        - 6 for subsystem = a nonlinear molecule, nonperiodic calculation
        - 3 in periodic calculations
        """
        # determine nb of zeros
        subs3 = sum([[3*at, 3*at+1, 3*at+2] for at in self.subs],[])
        U, W, Vt = np.linalg.svd(np.take(molecule.external_basis,subs3,1), full_matrices=False)
        rank = (abs(W) > abs(W[0])*self.svd_threshold).sum()
        self.num_zeros = rank

        # check
        if self.num_zeros not in [3,5,6] and not molecule.periodic :
            raise ValueError("Number of zeros is expected to be 3, 5 or 6, but found %i." % self.num_zeros)
        if self.num_zeros != 3 and molecule.periodic :
            raise ValueError("Number of zeros is expected to be 3 (periodic calculation), but found %i." % self.num_zeros)

        # return mass-weighted basis vectors for external degrees of freedom
        if do_modes:
            foo = np.dot(U.transpose(), molecule.external_basis)
            self.external_basis = foo[:rank]/W[:rank].reshape((-1,1))

        #---- other implementation ----
        #if not molecule.periodic:
        #    self.num_zeros = rank_linearity(np.take(molecule.coordinates,self.subs,0), svd_threshold = self.svd_threshold)
        #else:
        #    self.num_zeros = 3
        #if do_modes and self.num_zeros > 0:
        #    if self.num_zeros == 3:
        #        self.external_basis = molecule.external_basis[:3,:]  # three translations
        #    elif self.num_zeros == 5:
        #        # Compute direction of the linear SUBSystem (with two atoms) and check for highest alignment with one of the axes.
        #        diff = molecule.coordinates[self.subs[0]] - molecule.coordinates[self.subs[1]]
        #        axis = 3+abs(diff).tolist().index(max(abs(diff)))   # axis to be skipped
        #        alphas = [i for i in range(6) if i is not axis]
        #        self.external_basis = np.take(molecule.external_basis,alphas,0)
        #    elif self.num_zeros == 6:
        #        self.external_basis = molecule.external_basis
        #    else:
        #        raise ValueError("Number of zeros is expected to be 3, 5 or 6, but found %i." % self.num_zeros)


    def compute_hessian(self, molecule, do_modes):
        """See :meth:`Treatment.compute_hessian`.

        The VSA Hessian reads: ``H_ss - H_se (H_ee)**(-1) H_es``,
        the VSA mass matrix reads: ``M_s - H_se (H_ee)**(-1) M_e (H_ee)**(-1) H_es``,
        where the indices ``s`` and ``e`` refer to the subsystem and environment
        atoms respectively.
        """
        # fill lists with subsystem/environment atoms/coordinates
        subs = self.subs.tolist()
        envi = sum([[at] for at in xrange(molecule.size) if at not in subs],[])
        subs3 = sum([[3*at, 3*at+1, 3*at+2] for at in subs],[])
        envi3 = sum([[3*at, 3*at+1, 3*at+2] for at in xrange(molecule.size) if at not in subs],[])

        # 1. Construct Hessian (small: 3Nsubs x 3Nsubs)
        # construct H_ss, H_ee, H_es
        hessian_ss = np.take(np.take(molecule.hessian,subs3,0), subs3,1)
        hessian_ee = np.take(np.take(molecule.hessian,envi3,0), envi3,1)
        hessian_es = np.take(np.take(molecule.hessian,envi3,0), subs3,1)
        # construct H_ee**-1 and H_ee**-1 . H_es
        hessian_e1 = np.linalg.inv(hessian_ee)
        hessian_e1_es = np.dot(hessian_e1,hessian_es)
        # construct H_ss - H_se . H_ee**-1 . H_es
        self.hessian_small = hessian_ss - np.dot( hessian_es.transpose(), hessian_e1_es)

        # 2. Construct mass matrix (small: 3Nsubs x 3Nsubs)
        # with corrected mass matrix
        masses3_subs = molecule.masses3[subs3]               # masses subsystem
        masses3_envi = molecule.masses3[envi3]               # masses environment
        tempmat = np.zeros((len(envi3),len(subs3)),float) # temporary matrix

        # construct   M_e . H_ee**-1 . H_es
        tempmat = masses3_envi.reshape((-1,1))*hessian_e1_es
        # construct   H_se . H_ee**-1 . M_e . H_ee**-1 . H_es
        massmatrixsmall = np.dot(hessian_e1_es.transpose(), tempmat)
        # construct   M_s + H_se . H_ee**-1 . M_e . H_ee**-1 . H_es
        # by adding the diagonal contributions
        massmatrixsmall.ravel()[::len(massmatrixsmall)+1] += masses3_subs
        self.mass_matrix_small = MassMatrix( massmatrixsmall )

        if do_modes:
            atom_division = AtomDivision(envi+subs,[],[])
            self.transform = Transform( np.concatenate( (- hessian_e1_es, np.identity(len(subs3))),0), atom_division)


class VSANoMass(Treatment):
    """
    Perform a Vibrational Subsystem Analysis, without taking into account the
    mass of the environment.

    Frequencies and modes are computed as described in the reference:

    - W. Zheng, B.R. Brooks, Journal of Biophysics 89, 167 (2006)
    - A. Ghysels, V. Van Speybroeck, E. Pauwels, S. Catak, B.R. Brooks,
      D. Van Neck, M. Waroquier, Journal of Computational Chemistry 31 (5),
      94-1007 (2010)

    The system is partitioned into a subsystem and an environment. The subsystem
    atoms are allowed to vibrate, while the environment atoms follow the motions
    of the subsystem atoms. The environment atoms are force free.
    Moreover, the VSA is performed according to the original version of 2006:
    no mass correction for the environment is included.
    This version of VSA corresponds to the approximation
    of zero mass for all environment atoms.
    """
    def __init__(self, subs, svd_threshold=1e-5):
        """
           One argument:
            | ``subs`` -- a list with the subsystem atoms, counting starts from
                          zero.

           Optional argument:
            | ``svd_threshold`` -- threshold for detection of deviations for
                                   linearity
        """
        # QA:
        if len(subs) == 0:
            raise ValueError("At least one subsystem atom is required.")
        # Rest of init:
        self.subs = np.array(subs)
        #self.subs.sort()
        self.svd_threshold = svd_threshold
        Treatment.__init__(self)

    def compute_zeros(self, molecule, do_modes):

        """See :meth:`Treatment.compute_zeros`.

        The number of zeros should be:

        - 3 for subsystem = a single atom, nonperiodic calculation
        - 5 for subsystem = a linear molecule, nonperiodic calculation
        - 6 for subsystem = a nonlinear molecule, nonperiodic calculation
        - 3 in periodic calculations
        """
        # determine nb of zeros
        subs3 = sum([[3*at, 3*at+1, 3*at+2] for at in self.subs],[])
        U, W, Vt = np.linalg.svd(np.take(molecule.external_basis,subs3,1), full_matrices=False)
        rank = (abs(W) > abs(W[0])*self.svd_threshold).sum()
        self.num_zeros = rank

        # check
        if self.num_zeros not in [3,5,6] and not molecule.periodic :
            raise ValueError("Number of zeros is expected to be 3, 5 or 6, but found %i." % self.num_zeros)
        if self.num_zeros != 3 and molecule.periodic :
            raise ValueError("Number of zeros is expected to be 3 (periodic calculation), but found %i." % self.num_zeros)

        # return mass-weighted basis vectors for external degrees of freedom
        if do_modes:
            self.external_basis = np.zeros((rank, molecule.size*3), float)
            self.external_basis[:,subs3] = Vt[:rank]

        #---- other implementation ----
        #if not molecule.periodic:
        #    self.num_zeros = rank_linearity(np.take(molecule.coordinates,self.subs,0), svd_threshold = self.svd_threshold)
        #else:
        #    self.num_zeros = 3
        #if do_modes and self.num_zeros > 0:
        #    if self.num_zeros == 3:
        #        self.external_basis = molecule.external_basis[:3,:]  # three translations
        #    elif self.num_zeros == 5:
        #        # Compute direction of the linear SUBSystem (with two atoms) and check for highest alignment with one of the axes.
        #        diff = molecule.coordinates[self.subs[0]] - molecule.coordinates[self.subs[1]]
        #        axis = 3+abs(diff).tolist().index(max(abs(diff)))   # axis to be skipped
        #        alphas = [i for i in range(6) if i is not axis]
        #        self.external_basis = np.take(molecule.external_basis,alphas,0)
        #    elif self.num_zeros == 6:
        #        self.external_basis = molecule.external_basis
        #    else:
        #        raise ValueError("Number of zeros is expected to be 3, 5 or 6, but found %i." % self.num_zeros)


    def compute_hessian(self, molecule, do_modes):
        """See :meth:`Treatment.compute_hessian`.

        The VSANoMass Hessian reads: ``H_ss - H_se (H_ee)**(-1) H_es``.
        and the VSANoMass mass matrix reads: ``M_s``,
        where the indices ``s`` and ``e`` refer to the subsystem and environment
        atoms respectively.
        """
        # fill lists with subsystem/environment atoms/coordinates
        subs = self.subs.tolist()
        envi = sum([[at] for at in xrange(molecule.size) if at not in subs],[])
        subs3 = sum([[3*at, 3*at+1, 3*at+2] for at in subs],[])
        envi3 = sum([[3*at, 3*at+1, 3*at+2] for at in xrange(molecule.size) if at not in subs],[])

        # 1. Construct Hessian (small: 3Nsubs x 3Nsubs)
        # construct H_ss, H_ee, H_es
        hessian_ss = np.take(np.take(molecule.hessian,subs3,0), subs3,1)
        hessian_ee = np.take(np.take(molecule.hessian,envi3,0), envi3,1)
        hessian_es = np.take(np.take(molecule.hessian,envi3,0), subs3,1)
        # construct H_ee**-1 and H_ee**-1 . H_es
        hessian_e1 = np.linalg.inv(hessian_ee)
        hessian_e1_es = np.dot(hessian_e1,hessian_es)
        # construct H_ss - H_se . H_ee**-1 . H_es
        self.hessian_small = hessian_ss - np.dot( hessian_es.transpose(), hessian_e1_es)

        # 2. Construct mass matrix (small: 3Nsubs x 3Nsubs)
        # with plain submatrix M_s
        self.mass_matrix_small = MassMatrix( np.diag(np.take(molecule.masses3,subs3)) )

        if do_modes:
            atom_division = AtomDivision(envi+subs,[],[])
            self.transform = Transform( np.concatenate( (- hessian_e1_es, np.identity(len(subs3))),0), atom_division)


class MBH(Treatment):
    """
    The Mobile Block Hessian approach.

    Frequencies and modes are computed with the MBH approach, as described in
    the refences:

    * "Vibrational modes in partially optimized molecular systems", An Ghysels,
      Dimitri Van Neck, Veronique Van Speybroeck, Toon Verstraelen and Michel
      Waroquier, Journal of Chemical Physics, Vol. 126 (22), Art. No. 224102,
      2007, http://dx.doi.org/1.2737444

    * "Cartesian formulation of the Mobile Block Hesian Approach to vibrational
      analysis in partially optimized systems", An Ghysels, Dimitri Van Neck and
      Michel Waroquier, Journal of Chemical Physics, Vol. 127 (16), Art. No.
      164108, 2007, http://dx.doi.org/1.2789429

    * "Calculating reaction rates with partial Hessians: validation of the MBH
      approach", An Ghysels, Veronique Van Speybroeck, Toon Verstraelen, Dimitri
      Van Neck and Michel Waroquier, Journal of Chemical Theory and Computation,
      Vol. 4 (4), 614-625, 2008, http://dx.doi.org/ct7002836

    For the Mobile Block Hessian method with linked blocks, please refer to the
    following papers:

    * "Mobile Block Hessian approach with linked blocks: an efficient approach
      for the calculation of frequencies in macromolecules", An Ghysels,
      Veronique Van Speybroeck, Ewald Pauwels, Dimitri Van Neck, Bernard R.
      Brooks and Michel Waroquier, Journal of Chemical Theory and Computation,
      Vol. 5 (5), 1203-1215, 2009, http://dx.doi.org/ct800489r

    * "Normal modes for large molecules with arbitrary link constraints in the
      mobile block Hessian approach", An Ghysels, Dimitri Van Neck, Bernard R.
      Brooks, Veronique Van Speybroeck and Michel Waroquier, Journal of Chemical
      Physics, Vol. 130 (18), Art. No. 084107, 2009, http://dx.doi.org/1.3071261

    The system is partitioned into blocks which are only allowed to move as
    rigid bodies during the vibrational analysis. Atoms that are not part of a
    block can still move individually. The internal geometry of the blocks need
    not be optimized, because the MBH method performs a gradient correction to
    account for the internal forces. Only the position and orientation of each
    block should be optimized. This make MBH an appropriate method to perform
    NMA in partially optimized structures.
    """
    def __init__(self, blocks, do_gradient_correction=True, svd_threshold=1e-5):
        """
           One argument:
            | ``blocks`` -- a list of blocks, each block is a list of atoms,
                            counting starts from zero.

           Optional arguments:
            | ``do_gradient_correction`` -- boolean, whether gradient correction
                                            to MBH should be added
            | ``svd_threshold`` -- threshold for zero singular values in svd
        """
        # QA:
        if len(blocks) == 0:
            raise ValueError("At least one block is required.")
        if not isinstance(do_gradient_correction, bool):
            raise TypeError("Optional argument do_gradient_correction should be boolean.")
        # Rest of init:
        self.blocks = blocks
        self.svd_threshold = svd_threshold
        self.do_gradient_correction = do_gradient_correction
        Treatment.__init__(self)

    def compute_zeros(self, molecule, do_modes):
        """See :meth:`Treatment.compute_zeros`.

        The number of zeros should be:

        - 3 for a single atom, nonperiodic calculation
        - 5 for a linear molecule, nonperiodic calculation
        - 6 for a nonlinear molecule, nonperiodic calculation
        - 3 in periodic calculations
        """
        # determine nb of zeros
        U, W, Vt = np.linalg.svd(molecule.external_basis, full_matrices=False)
        rank = (abs(W) > abs(W[0])*self.svd_threshold).sum()
        self.num_zeros = rank

        # check
        if self.num_zeros not in [3,5,6] and not molecule.periodic :
            raise ValueError("Number of zeros is expected to be 3, 5 or 6, but found %i." % self.num_zeros)
        if self.num_zeros != 3 and molecule.periodic :
            raise ValueError("Number of zeros is expected to be 3 (periodic calculation), but found %i." % self.num_zeros)

        # return mass-weighted basis vectors for external degrees of freedom
        if do_modes:
            self.external_basis = Vt[:rank]

        #---- other implementation ----
        #if not molecule.periodic:
        #    self.num_zeros = rank_linearity(molecule.coordinates, svd_threshold = self.svd_threshold)
        #else:
        #    self.num_zeros = 3
        #if do_modes and self.num_zeros > 0:
        #    if self.num_zeros == 3:
        #        self.external_basis = molecule.external_basis[:3,:]
        #    elif self.num_zeros == 5:
        #        # Compute direction of the linear SYSTEM (with two atoms) and check for highest alignment with one of the axes.
        #        diff = molecule.coordinates[0,:] - molecule.coordinates[1,:]
        #        axis = 3+abs(diff).tolist().index(max(abs(diff)))   # axis to be skipped
        #        alphas = [i for i in range(6) if i is not axis]
        #        self.external_basis = np.take(molecule.external_basis,alphas,0)
        #    elif self.num_zeros == 6:
        #        self.external_basis = molecule.external_basis
        #    else:
        #        raise ValueError("Number of zeros is expected to be 3, 5 or 6, but found %i." % self.num_zeros)

        #---- other implementation ----
        #U, W, Vt = np.linalg.svd(molecule.external_basis, full_matrices=False)
        #rank = (abs(W) > abs(W[0])*self.svd_threshold).sum()
        #self.num_zeros = rank
        #if do_modes:
        #    self.external_basis = Vt[:rank]


    def compute_hessian(self, molecule, do_modes):
        """See :meth:`Treatment.compute_hessian`.

        Gather all information about the block choice in the
        blkinfo attribute. If adjoined blocks (blocks with common atoms)
        are present, an extra STRICT block choice is defined: the partitioning
        where each atom belongs to only one block.

        First, the 3Nxd matrix U is constructed from the column vectors describing
        the block motions. At this point, the strict partitioning is used.
        The Hessian is equal to: ``Hp = U^T H U` + G:C`. The term ``G:C`` is the
        gradient correction. The mass matrix is equal to: ``Mp = U^T M U``.

        Second, impose the link constraints between blocks, if any of the blocks
        are adjoined. The matrix K constains the linking constraints. If
        nullspace is the null space of K, then the final Mobile block Hessian
        reads: ``Hy = nullspace^T Hp nullspace + Gp:Cp``. The term ``Gp:Cp``
        is the gradient correction. The Mobile Block mass matrix is equal to:
        ``My = nullspace^T Mp nullspace``.
        """
        # Notation: b,b0,b1   --  a block index
        #           block  --  a list of atoms, e.g. [at1,at4,at6]
        #           alphas  --  the 6 block parameter indices (or 5 for linear block)

        # Block information
        blkinfo = Blocks(self.blocks, molecule, self.svd_threshold)
        mbhdim1 = 6*blkinfo.nb_nlin + 5*blkinfo.nb_lin + 3*len(blkinfo.free)

        # TRANSFORM from CARTESIAN to BLOCK PARAMETERS
        U = self._construct_U(molecule,mbhdim1,blkinfo)

        # Construct Hessian in block parameters: Hp = U**T . H . U + correction
        Hp = np.dot(np.dot( U.transpose(), molecule.hessian) , U)

        # gradient correction
        if self.do_gradient_correction:
         for b,block in enumerate(blkinfo.blocks_nlin_strict+blkinfo.blocks_lin_strict):   # Nonlinear AND linear blocks
            GP = np.zeros((3,3),float)
            G = np.take(molecule.gradient,block,0)
            P = np.take(molecule.coordinates,block,0)
            for i in range(len(block)):
                GP += np.dot(G[i,:].reshape((3,1)),P[i,:].reshape((1,3)))
            # note: GP is not symmetric
            p = GP[0,0] # Gx.x
            q = GP[1,1] # Gy.y
            r = GP[2,2] # Gz.z
            s = GP[1,0] # Gy.x
            t = GP[2,0] # Gz.x
            u = GP[2,1] # Gz.y

            corr    = np.zeros((6,6),float)
            corr[3:,3:] = np.diag([-q-r, -p-r, -p-q])
            corr[3,4]   =  corr[4,3] = s
            corr[3,5]   =  corr[5,3] = t
            corr[4,5]   =  corr[5,4] = u

            if b < blkinfo.nb_nlin:             # Nonlinear block
                alphas = range(6)
                col = 6*b
                dim = 6
                Hp[col:(col+dim),col:(col+dim)] += corr

            else:                               # Linear block
                b   = b - blkinfo.nb_nlin       # reset
                col = 6*blkinfo.nb_nlin + 5*b   # offset
                dim = 5

                alphas = [index for index in range(6) if index != blkinfo.skip_axis_lin[b]]
                Hp[col:(col+dim),col:(col+dim)] += np.take(np.take(corr,alphas,0),alphas,1)

        # Construct mass matrix in block parameters: Mp = U**T . M . U
        Mp = np.dot(U.transpose(),  U * molecule.masses3.reshape((-1,1)))

        if blkinfo.is_linked:
            # SECOND TRANSFORM: from BLOCK PARAMETERS to Y VARIABLES
            # Necessary if blocks are linked to each other.
            nullspace = self._construct_nullspace_K(molecule,mbhdim1,blkinfo)

            My = np.dot(nullspace.transpose(), np.dot( Mp,nullspace) )
            Hy = np.dot(nullspace.transpose(), np.dot( Hp,nullspace) )

            # TODO
            # gradient correction of the second transform...


        if do_modes:
            if not blkinfo.is_linked:
                self.hessian_small = Hp
                self.mass_matrix_small = MassMatrix(Mp)
                self.transform = Transform(U)
            else:
                self.hessian_small = Hy
                self.mass_matrix_small = MassMatrix(My)
                self.transform = Transform(np.dot(U, nullspace))
        else:
            if not blkinfo.is_linked:
                self.hessian_small = Hp
                self.mass_matrix_small = MassMatrix(Mp)
            else:
                self.hessian_small = Hy
                self.mass_matrix_small = MassMatrix(My)

    def _construct_U(self,molecule,mbhdim1,blkinfo):
        # Construct first transformation matrix
        D = transrot_basis(molecule.coordinates)   # is NOT mass-weighted

        U = np.zeros((3*molecule.size, mbhdim1),float)

        for b,block in enumerate(blkinfo.blocks_nlin_strict):
            for at in block:
                for alpha in range(6):
                    U[3*at:3*(at+1), 6*b+alpha] = D[alpha,3*at:3*(at+1)]

        for b,block in enumerate(blkinfo.blocks_lin_strict):
            for at in block:
                alphas = [index for index in range(6) if index != blkinfo.skip_axis_lin[b]]
                for i,alpha in enumerate(alphas):
                    col = 6*blkinfo.nb_nlin + 5*b + i
                    U[3*at:3*(at+1), col] = D[alpha,3*at:3*(at+1)]

        for i,at in enumerate(blkinfo.free):
            for mu in range(3):
                U[3*at+mu, 6*blkinfo.nb_nlin+5*blkinfo.nb_lin+3*i +mu] = 1.0
        return U

    def _construct_nullspace_K(self,molecule,mbhdim1,blkinfo):
        # SECOND TRANSFORM: from BLOCK PARAMETERS to Y VARIABLES
        # Necessary if blocks are linked to each other.
        # Construct K matrix, with constraints
        D = transrot_basis(molecule.coordinates)   # is NOT mass-weighted
        nbrows = (np.sum(blkinfo.sharenbs)-molecule.size)*3
        K = np.zeros(( nbrows, mbhdim1-3*len(blkinfo.free)), float)
        row = 0
        for (at,apps) in blkinfo.appearances.iteritems():
            if len(apps) >= 2:
                # the first block
                b0 = apps[0]
                D0 = D[:,3*at:3*(at+1)]

                if b0 < blkinfo.nb_nlin:  # if b0 is nonlinear block
                    sta0 = 6*b0
                    end0 = 6*(b0+1)
                    D0 = D[:,3*at:3*(at+1)]

                else:   # if b0 is a linear block
                    b0 = b0 - blkinfo.nb_nlin     # reset
                    sta0 = 6*blkinfo.nb_nlin + 5*b0       # offset
                    end0 = 6*blkinfo.nb_nlin + 5*(b0+1)
                    alphas = [index for index in range(6) if index != blkinfo.skip_axis_lin[b0]]
                    D0 = np.take(D[:,3*at:3*(at+1)],alphas,0)

                for b1 in apps[1:]:
                    # add 3 rows to K, for each block connected to b0
                    K[row:row+3, sta0:end0] = D0.transpose()

                    if b1 < blkinfo.nb_nlin:  # if b1 is nonlinear block
                        sta1 = 6*b1
                        end1 = 6*(b1+1)
                        D1 = D[:,3*at:3*(at+1)]

                    else:   # if b1 is a linear block
                        b1 = b1 - blkinfo.nb_nlin     # reset
                        sta1 = 6*blkinfo.nb_nlin + 5*b1       # offset
                        end1 = 6*blkinfo.nb_nlin + 5*(b1+1)
                        alphas = [index for index in range(6) if index != blkinfo.skip_axis_lin[b1]]
                        D1 = np.take(D[:,3*at:3*(at+1)],alphas,0)

                    K[row:row+3,sta1:end1] = -D1.transpose()
                    row += 3

        # Do SVD of matrix K
        u,s,vh = np.linalg.svd(K)

        # construct nullspace of K
        rank = sum(s>max(s)*self.svd_threshold)
        nullspace = vh[rank:,:].transpose()
        [r_null,c_null] = nullspace.shape
        n = np.zeros((mbhdim1,c_null+3*len(blkinfo.free)),float)
        n[:r_null,:c_null] = nullspace
        n[r_null:,c_null:] = np.identity(3*len(blkinfo.free),float)
        return n


class MBHConstrainExt(MBH):
    """The Mobile Block Hessian approach with the Eckart constraints imposed

       This method is completely similar to the MBH, except that first the
       global translations and rotations are first projected out of the
       Hessian before applying the block partitioning and projecting by
       the MBH. The contribution of the gradient is also adapted.
       In case of a periodic simulation, only the global translations
       are projected out."""
    def __init__(self, blocks, do_gradient_correction=True, svd_threshold=1e-5):
        MBH.__init__(self,blocks)

    def compute_zeros(self, molecule, do_modes):
        MBH.compute_zeros(self,molecule,do_modes)

    def compute_hessian(self, molecule, do_modes):
        # perform projection of Hessian and gradient
        D = transrot_basis(molecule.coordinates, rot=molecule.periodic).transpose()
        for i in range(D.shape[1]):
            D[:,i] /= np.sqrt(np.sum(D[:,i]**2))
        proj = np.identity(D.shape[0]) - np.dot(D,D.transpose())
        hessian = np.dot(proj,molecule.hessian)
        gradient = (np.dot(proj,molecule.gradient.reshape(3*molecule.size,-1))).reshape(molecule.size,3)
        # construct a new Molecule instance
        mol = Molecule(molecule.numbers, molecule.coordinates, molecule.masses,
                       molecule.energy, gradient, hessian, molecule.multiplicity,
                       periodic=molecule.periodic)
        # do the usual MBH
        MBH.compute_hessian(self,mol,do_modes)


class Blocks(object):
    """Object that extracts all information from a block choice."""
    def __init__(self,blocks,molecule,svd_threshold):
        """
        Arguments:
         | ``blocks`` -- a list of lists of atoms
                         [ [at1,at5,at3], [at4,at5], ...]
                         with a list of atoms for each block
         | ``molecule`` -- Molecule object, necessary for N (total nb of atoms)
                           and positions (linearity of blocks).
         | ``svd_trheshold`` -- threshold for zero singular values in svd
        """
        N = molecule.size
        # check for empty blocks and single-atom-blocks
        to_remove = []
        for b,block in enumerate(blocks):
            if len(block)==0:
                to_remove.append(b)
            elif len(block)==1:
                to_remove.append(b)
            elif max(block)>=N:
                raise ValueError("block "+str(b)+": atoms should be in range [0,N-1], N="+str(N))
            elif min(block)<0:
                raise ValueError("block "+str(b)+": atoms should be in range [0,N-1], N="+str(N))
            elif len(block) != len(set(block)):
                raise ValueError("Duplicate atoms encountered in block %i." % b)
        #remove single atoms and empty blocks
        for i in range(len(to_remove)):
            del blocks[to_remove[len(to_remove)-i-1]]   #remove starting from largests b

        # define fixed atoms (in one or more blocks) and free atoms (not in blocks)
        fixed = set([])
        for block in blocks:
            fixed.update(block)
        free  = set(range(N)) - fixed
        # turn them into lists and sort them
        fixed = sorted(fixed)
        free = sorted(free)

        # check for linearity and fill in dimensions
        #D = molecule.external_basis
        #dim_block=np.zeros((len(blocks)),int)
        indices_blocks_nlin = []    # nonlinear blocks
        indices_blocks_lin  = []    # linear blocks
        for b,block in enumerate(blocks):
            rank = rank_linearity(np.take(molecule.coordinates,block,0), svd_threshold=svd_threshold)[0]
            if rank==6:    indices_blocks_nlin.append(b)
            elif rank==5:  indices_blocks_lin.append(b)
            else:          raise ValueError("In principle rank should have been 5 or 6, found "+str(rank))

        # REORDER THE BLOCKS: nonlinear blocks, linear blocks, single-atom-blocks
        blocks_nlin = []
        blocks_lin  = []
        for b in indices_blocks_nlin:
            blocks_nlin.append(blocks[b])
        for b in indices_blocks_lin:
            blocks_lin.append(blocks[b])
        orderedblocks = blocks_nlin + blocks_lin
        nb_nlin = len(blocks_nlin)
        nb_lin  = len(blocks_lin)

        # fill in appearances
        appearances = {}
        for b,block in enumerate(orderedblocks):
            for atom in block:
                # If atom is not already in appearances,
                # setdefault sets appearances[atom] to an empty list.
                appearances.setdefault(atom,[]).append(b)
        for i,atom in enumerate(free):
            appearances.setdefault(atom,[]).append( i+len(orderedblocks) )

        # make a strict partition of the atoms: each atom belongs to one block only
        bA1 = np.zeros((molecule.size),int)
        for (at,apps) in appearances.iteritems():
            bA1[at] = apps[0]

        blocks_nlin_strict = []
        for b,block in enumerate(blocks_nlin):
            atoms = []
            for at in block:
                if bA1[at] == b:
                    atoms.append(at)
            blocks_nlin_strict.append(atoms)

        blocks_lin_strict = []
        for b,block in enumerate(blocks_lin):
            atoms = []
            for at in block:
                if bA1[at] == (b+nb_nlin):
                    atoms.append(at)
            blocks_lin_strict.append(atoms)

        # for linear blocks: axis?
        # Compute direction of the linear block (with two atoms) and check for highest
        # alignment with one of the axes.
        skip_axis_lin = np.zeros((nb_lin),int)
        for b,block in enumerate(blocks_lin):       # do not use strict partition here
            diff = molecule.coordinates[block[0]] - molecule.coordinates[block[1]]
            skip_axis_lin[b] = 3+abs(diff).tolist().index(max(abs(diff)))   # axis to be skipped


        # Check if there are linked blocks
        sharenbs = np.zeros((molecule.size),int)  # share number of each atom
        for (at,apps) in appearances.iteritems():
            sharenbs[at] = len(apps)
        is_linked = False
        for sharenb in np.ravel(sharenbs):
            if sharenb > 1:
                is_linked = True
                break


        self.blocks = orderedblocks
        self.free   = free
        self.nb_nlin   = nb_nlin
        self.nb_lin    = nb_lin
        self.blocks_nlin = blocks_nlin
        self.blocks_lin  = blocks_lin

        self.blocks_nlin_strict = blocks_nlin_strict
        self.blocks_lin_strict  = blocks_lin_strict

        self.appearances = appearances
        self.bA1 = bA1
        self.skip_axis_lin = skip_axis_lin

        self.sharenbs = sharenbs
        self.is_linked = is_linked


class PHVA_MBH(MBH):
    """
    The Mobile Block Hessian combined with the Partial Hessian Vibrational Analysis.

    The system is partitioned into

    - blocks which are only allowed to move as rigid bodies (the MBH concept,
      see :class:`PHVA`),
    - fixed atoms which are not allowed to move at all (the PHVA concept, see
      :class:`MBH`).
    - single atoms which are allowed to vibrate individually

    The internal geometry of the blocks does not have to be optimized. The
    positions of the fixed atoms also do not have to be optimized.
    This make the PHVA_MBH an appropriate method to perform NMA in partially
    optimized structures.
    """
    def __init__(self, fixed, blocks, do_gradient_correction=True, svd_threshold=1e-5):
        """
           Two arguments:
            | ``fixed`` -- a list with fixed atoms, counting starts from zero.
            | ``blocks`` -- a list of blocks, each block is a list of atoms

           Optional arguments:
            | ``svd_threshold`` -- threshold for zero singular values in svd
            | ``do_gradient_correction`` -- boolean, whether gradient correction
                                            to MBH part should be added
                                            [default=True]
        """
        # QA:
        if len(fixed) == 0:
            raise ValueError("At least one fixed atom is required.")
        # check if no atoms both in fixedatoms AND in blocks
        for block in blocks:
            for at in block:
                if at in fixed:
                    raise ValueError("Atoms in blocks can not be part of fixed atom region: atom "+str(at)+" is in both regions.")
        # Rest of init:
        self.fixed = np.array(fixed)
        MBH.__init__(self, blocks, do_gradient_correction=do_gradient_correction, svd_threshold=svd_threshold)

    def compute_zeros(self, molecule, do_modes):
        """See :meth:`Treatment.compute_zeros`"""
        # [ See explanation PHVA ]
        U, W, Vt = np.linalg.svd(molecule.external_basis, full_matrices=False)
        rank = (abs(W) > abs(W[0])*self.svd_threshold).sum()
        external_basis = Vt[:rank]
        fixed3 = []
        for i in self.fixed:
            fixed3.append(3*i)
            fixed3.append(3*i+1)
            fixed3.append(3*i+2)
        system = external_basis[:,fixed3].transpose()
        U, W, Vt = np.linalg.svd(system, full_matrices=False)
        self.num_zeros = (abs(W) < abs(W[0])*self.svd_threshold).sum()
        if do_modes and self.num_zeros > 0:
            # TODO: fix this
            #self.exernal_basis = ...
            raise NotImplementedError

    def compute_hessian(self, molecule,do_modes):
        """See :meth:`Treatment.compute_hessian`.

        First, the non-fixed atoms are cut out of the molecular system, i.e.
        the single atoms and the atoms belonging to the blocks and stored in a
        *submolecule*.
        Next, the regular MBH concept is applied on this submolecule.
        """

        # Make submolecule
        selectedatoms = [at for at in xrange(molecule.size) if at not in self.fixed]
        selectedcoords = sum([[3*at,3*at+1,3*at+2] for at in selectedatoms],[])

        from tamkin.data import Molecule
        submolecule = Molecule(
            np.take(molecule.numbers, selectedatoms),
            np.take(molecule.coordinates, selectedatoms, 0),
            np.take(molecule.masses, selectedatoms),
            molecule.energy,
            np.take(molecule.gradient,selectedatoms,0),
            np.take(np.take(molecule.hessian,selectedcoords,0),selectedcoords,1),
            molecule.multiplicity,
            0, # undefined molecule.symmetry_number
            molecule.periodic
        )

        # adapt numbering in blocks
        shifts = np.zeros((molecule.size),int)
        for fixat in self.fixed:
            shifts[fixat:] = shifts[fixat:]+1
        for bl,block in enumerate(self.blocks):
            for at,atom in enumerate(block):
                self.blocks[bl][at] = atom - shifts[atom]

        MBH.compute_hessian(self, submolecule, do_modes)

        if do_modes:   # adapt self.transform to include the fixed atom rows/cols
            transf = np.zeros((3*molecule.size, self.transform.matrix.shape[1]),float)
            transf[selectedcoords,:] = self.transform.matrix
            self.transform = Transform(transf)


class Constrain(Treatment):
    """Perform a normal mode analysis where part of the internal coordinates are
       constrained to a fixed value.

       The gradient corrections are taken into account correctly. At present, only
       distance constraints are implemented. In principle, the routine can be
       adapted to angle and dihedral angle constraints.
    """
    def __init__(self, constraints, do_gradient_correction=True, svd_threshold=1e-5):
        """
           One argument:
            | ``constraints`` -- a list with constraints of internal coordinates:
                                 [at1,at2] to constrain a distance,
                                 [at1,at2,at3] to constrain an angle,
                                 [at1,at2,at3,at4] to constrain a dihedral angle.

           Optional:
            | ``do_gradient_correction`` -- whether gradient correction should
                                            be applied
            | ``svd_threshold`` -- threshold for singular value decomposition
        """
        # QA:
        if len(constraints) == 0:
            raise ValueError("At least one constraint is required.")
        # Rest of init:
        self.constraints = constraints
        self.do_gradient_correction = do_gradient_correction
        self.svd_threshold = svd_threshold
        Treatment.__init__(self)

    def compute_zeros(self, molecule, do_modes):
        """See :meth:`Treatment.compute_zeros`.

        The number of zeros should be:

        - 3 for a single atom, nonperiodic calculation
        - 5 for a linear molecule, nonperiodic calculation
        - 6 for a nonlinear molecule, nonperiodic calculation
        - 3 in periodic calculations
        """
        # determine nb of zeros
        U, W, Vt = np.linalg.svd(molecule.external_basis, full_matrices=False)
        rank = (abs(W) > abs(W[0])*self.svd_threshold).sum()
        self.num_zeros = rank

        # check
        if self.num_zeros not in [3,5,6] and not molecule.periodic :
            raise ValueError("Number of zeros is expected to be 3, 5 or 6, but found %i." % self.num_zeros)
        if self.num_zeros != 3 and molecule.periodic :
            raise ValueError("Number of zeros is expected to be 3 (periodic calculation), but found %i." % self.num_zeros)

        # return mass-weighted basis vectors for external degrees of freedom
        if do_modes:
            self.external_basis = Vt[:rank]


    def compute_hessian(self, molecule, do_modes):
        """See :meth:`Treatment.compute_hessian`"""

        # make constraint matrix
        constrmat = np.zeros((3*molecule.size,len(self.constraints)))
        count = 0
        for i,constraint in enumerate(self.constraints):
            # constrain a distance
            if len(constraint) == 2:
                at1 = constraint[0]
                at2 = constraint[1]
                dist = np.sqrt(np.sum((molecule.coordinates[at1] - molecule.coordinates[at2])**2))
                dist = 1.0
                constrmat[3*at1:3*at1+3, count] = (molecule.coordinates[at1] - molecule.coordinates[at2])/dist
                constrmat[3*at2:3*at2+3, count] = (molecule.coordinates[at2] - molecule.coordinates[at1])/dist
                count += 1
        # determine the orthogonal complement of the basis of small
        # displacements determined by the constraints.
        U, W, Vt = np.linalg.svd(constrmat.transpose(), full_matrices=True)
        rank = (W/W[0] > self.svd_threshold).sum()
        nullspace = Vt[rank:].transpose()

        # mass matrix small = nullspace^T . M . nullspace
        # hessian     small = nullspace^T . H . nullspace + gradient correction
        self.mass_matrix_small = MassMatrix(np.dot(nullspace.transpose(),nullspace * molecule.masses3.reshape((-1,1))) )
        self.hessian_small     = np.dot(nullspace.transpose(), np.dot(molecule.hessian, nullspace))

        # check if gradient is small enough in this complement: overlap with nullspace should be small enough
        # print  np.sum(np.dot(nullspace.transpose(), np.ravel(molecule.gradient))**2)

        if self.do_gradient_correction:
            self._do_the_gradient_correction(constrmat.transpose(), U,W,Vt, rank, nullspace, molecule.gradient)

        if do_modes:
            self.transform = Transform(nullspace)


    def _do_the_gradient_correction(self, K,u,s,vh, rank, nullspace, gradient):

        # GRADIENT CORRECTION
        # construct generalized inverse L of K
        # zero rows for free atoms not included in L

        # construct generalized inverse
        L = np.zeros((K.shape[1],K.shape[0]),float)
        for i,sigma in enumerate(s[:rank]):
           L[i,:] = 1/sigma*u[:,i]
        L = np.dot(vh.transpose(),L)

        # construct right hand side of the equation
        Y = np.zeros((K.shape[0],len(self.hessian_small)**2))
        y = np.zeros((K.shape[0],1))
        count = 0
        for i in xrange(len(self.hessian_small)):
            for j in xrange(len(self.hessian_small)):
                for k,constraint in enumerate(self.constraints):
                    if len(constraint) == 2:
                        at1 = constraint[0]
                        at2 = constraint[1]
                        dist = 1.0
                        for a in range(3*at1,3*at1+3):
                            y[k] += 2 * nullspace[a,i]*nullspace[a,j]/dist
                        for b in range(3*at2,3*at2+3):
                            y[k] += 2 * nullspace[b,i]*nullspace[b,j]/dist
                        for (a,b) in zip(range(3*at1,3*at1+3),range(3*at2,3*at2+3)):
                            y[k] -= 2 * nullspace[a,i]*nullspace[b,j]/dist

                Y[:,count] = y[:,0]
                count += 1

        # solve equation
        X = -np.dot(L, Y)

        # apply correction
        count = 0
        for i in xrange(len(self.hessian_small)):
            for j in xrange(len(self.hessian_small)):

                self.hessian_small[i,j] += np.sum( np.ravel(gradient)*X[:,count] )
                count += 1
