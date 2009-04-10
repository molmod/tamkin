# TAMkin is a post-processing toolkit for thermochemistry and kinetics analysis.
# Copyright (C) 2008-2009 Toon Verstraelen <Toon.Verstraelen@UGent.be>,
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


from tamkin.io import load_chk, dump_chk

import numpy


__all__ = [
    "NMA", "AtomDivision", "Transform", "MassMatrix", "Treatment",
    "Full", "ConstrainExt", "PHVA", "VSA","VSA_no_mass",
]


class NMA(object):
    """A generic normal mode analysis class.

       This class gathers the functionality that is common between all types of
       nma variations, i.e. computation of frequencies and modes, once the
       problem is transformed to reduced coordinates. The actual nature of the
       reduced coordinates is determined by the treatment argument.
    """

    def __init__(self, molecule, treatment=None, do_modes=True):
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

        if do_modes:
            evals, modes_small_mw = numpy.linalg.eigh(hessian_small_mw)
        else:
            evals = numpy.linalg.eigvalsh(hessian_small_mw)
            modes_small_mw = None

        # frequencies
        self.freqs = numpy.sqrt(abs(evals))/(2*numpy.pi)
        self.freqs *= (evals > 0)*2-1

        if do_modes:
            # At this point the transform object transforms unweighted reduced
            # coordinates into Cartesian coordinates. Now we will alter it so that
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
                to_try = abs(self.freqs).argsort()[:num_try]
                overlaps = numpy.zeros(num_try, float)
                for counter, i in enumerate(to_try):
                    components = numpy.dot(treatment.external_basis, self.modes[:,i])
                    overlaps[counter] = numpy.linalg.norm(components)
                self.zeros = to_try[overlaps.argsort()[-treatment.num_zeros:]]
            else:
                self.zeros = abs(self.freqs).argsort()[:treatment.num_zeros]

        # a few more attributes that are worth keeping
        self.mass = molecule.mass
        self.masses = molecule.masses
        self.numbers = molecule.numbers
        self.coordinates = molecule.coordinates
        self.inertia_tensor = molecule.inertia_tensor
        self.multiplicity = molecule.multiplicity
        self.symmetry_number = molecule.symmetry_number
        self.periodic = molecule.periodic
        self.energy = molecule.energy

    def write_to_file(self, filename, fields='all'):
        if fields == 'all':
            data = dict((key, val) for key, val in self.__dict__.iteritems())
        elif fields == 'modes':
            keys = ["freqs", "modes", "masses", "numbers", "coordinates", "zeros"]
            data = dict((key, self.__dict__[key]) for key in keys)
        elif fields == 'partf':
            keys = ["freqs", "mass", "inertia_tensor", "multiplicity", "symmetry_number", "periodic", "energy", "zeros"]
            data = dict((key, self.__dict__[key]) for key in keys)
        dump_chk(filename, data)

    @classmethod
    def read_from_file(cls, filename):
        # ugly way to bypass the default constructor
        result = cls.__new__(cls)
        # load the file
        data = load_chk(filename)
        # check the names of the fields:
        possible_fields = set([
            "freqs", "modes", "mass", "masses", "numbers", "coordinates",
            "inertia_tensor", "multiplicity", "symmetry_number", "periodic",
            "energy", "zeros",
        ])
        if not set(data.iterkeys()).issubset(possible_fields):
            raise IOError("The Checkpoint file does not contain the correct fields.")
        # assign the attributes
        result.__dict__.update(data)
        return result


class AtomDivision(object):
    """A division of atoms into transformed, free and fixed.

       transformed: transformed in the reduced coordinates
       free:        identical in the reduced coordinates
       fixed:       absent in the reduced coordinates
    """

    def __init__(self, transformed, free, fixed):
        self.transformed = numpy.array(transformed, int)
        self.free = numpy.array(free, int)
        self.fixed = numpy.array(fixed, int)

        self.num_cartesian = 3*(len(self.transformed)+len(self.free)+len(self.fixed))
        self.to_cartesian_order = numpy.zeros(self.num_cartesian, int)
        self.to_reduced_order = numpy.zeros(self.num_cartesian, int)
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
        """Intialize the transformation object

           Arguments:
              atom_division -- set AtomDivision class
              matrix -- the linear transformation from the transformed
                        displacements to Cartesian coordinates.

           Attributes:
              matrix  --  see above
              scalars  --  ...
        """
        if matrix is None:
            matrix = numpy.zeros((0,0), float)
        if atom_division is None:
            # internal usage only:
            self._num_reduced = matrix.shape[1]
        else:
            # Quality Assurance:
            if matrix.shape[0] != 3*len(atom_division.transformed):
                raise ValueError("The matrix must have %i columns (matching the number of transformed atoms), got %i." %
                    3*len(transformed_atoms), matrix.shape[0]
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

    weighted = property(lambda self: self._weighted)

    def __call__(self, modes):
        """A transform object behaves like a function that transforms small
           displacement vectors from new to Cartesian coordinates.

           The array modes is a matrix with columns corresponding to the mass
           weighted modes in the reduced coordinates.
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
            return numpy.dot(self.matrix, modes)
        else:
            result = numpy.zeros((self.atom_division.num_cartesian, modes.shape[1]), float)
            i1 = 3*len(self.atom_division.transformed)
            i2 = i1 + 3*len(self.atom_division.free)
            result[:i1] = numpy.dot(self.matrix, modes[:self.matrix.shape[1]])
            if self.weighted:
                result[i1:i2] = modes[self.matrix.shape[1]:]*self.scalars
            else:
                result[i1:i2] = modes[self.matrix.shape[1]:]
            #    result[:,i2:] remains zero because these atoms are fixed
            # Reorder the atoms and return the result
            tmp = result[self.atom_division.to_cartesian_order]
            return tmp

    def make_weighted(self, mass_matrix):
        if self.weighted:
            raise Exception("The transformation is already weighted.")
        self.matrix = numpy.dot(self.matrix, mass_matrix.mass_block_inv_sqrt)
        self.scalars = mass_matrix.mass_diag_inv_sqrt.reshape((-1,1))
        self._weighted = True


class MassMatrix(object):
    """A clever mass matrix object. It is sparse when atom coordinates remain
       Cartesian in the reduced coordinates.

    """

    def __init__(self, *args):
        """Initialize the mass matrix object

           Arguments, if one is given and it is a two-dimensional matrix:
              mass_block -- the mass matrix associated with the transformed
                            coordinates

           Arguments, if one is given and it is a one-dimensional matrix:
              mass_diag -- the diagonal of the mass matrix associated with the
                           free atoms (each mass appears three times)

           Arguments, if two are given:
              mass_block -- the mass matrix associated with the transformed
                            coordinates
              mass_diag -- the diagonal of the mass matrix associated with the
                           free atoms (each mass appears three times)

           The mass of the fixed atoms does not really matter here.
        """
        if len(args) == 1:
            if len(args[0].shape) == 1:
                self.mass_diag = args[0]
                self.mass_block = numpy.zeros((0,0), float)
            elif len(args[0].shape) == 2:
                self.mass_diag = numpy.zeros((0,), float)
                self.mass_block = args[0]
            else:
                raise TypeError("When MassMatrix.__init__ gets one argument, it must be a one- or two-dimensional array.")
        elif len(args) == 2:
            self.mass_block = mass_block
            self.mass_diag = mass_diag
        else:
            raise TypeError("MassMatrix.__init__ takes one or two arguments, %i given." % len(args))

        # the square root of the inverse
        if len(self.mass_block) == 0:
            self.mass_block_inv_sqrt = numpy.zeros((0,0), float)
        else:
            evals, evecs = numpy.linalg.eigh(self.mass_block)
            self.mass_block_inv_sqrt = numpy.dot(evecs/numpy.sqrt(evals), evecs.transpose())
        self.mass_diag_inv_sqrt = 1/numpy.sqrt(self.mass_diag)

    def get_weighted_hessian(self, hessian):
        hessian_mw = numpy.zeros(hessian.shape,float)
        n = len(self.mass_block)
        # transform block by block:
        hessian_mw[:n,:n] = numpy.dot(numpy.dot(self.mass_block_inv_sqrt, hessian[:n,:n]), self.mass_block_inv_sqrt)
        hessian_mw[:n,n:] = numpy.dot(self.mass_block_inv_sqrt, hessian[:n,n:])*self.mass_diag_inv_sqrt
        hessian_mw[n:,:n] = hessian[:n,n:].transpose()
        hessian_mw[n:,n:] = (hessian[n:,n:]*self.mass_diag_inv_sqrt).transpose()*self.mass_diag_inv_sqrt
        return hessian_mw


class Treatment(object):
    """An abstract base class for the treatments. Derived classes must
       override the __call__ function. Parameters specific for the treatment
       are passed to the constructor, see for example the PHVA implementation."""
    def __init__(self):
        self.hessian_small = None
        self.mass_matrix_small = None
        self.transform = None
        self.num_zeros = None
        self.external_basis = None

    def __call__(self, molecule, do_modes):
        self.compute_hessian(molecule, do_modes)
        self.compute_zeros(molecule, do_modes)

    def compute_zeros(self, molecule, do_modes):
        # to be computed in derived classes:
        # treatment.num_zeros:
        #    the number of zero eigenvalues to expect
        # treatment.external_basis: (None if do_modes=False)
        #    the basis of external degrees of freedom. number of basis vectors
        #    matches the number of zeros.
        raise NotImplementedError

    def compute_hessian(self, molecule, do_modes):
        # to be computed in derived classes:
        # treatment.hessian_small:
        #    the Hessian in reduced coordinates
        # treatment.mass_matrix_small:
        #    the mass matrix in reduced coordinates (see MassMatrix class)
        # treatment.transform: (None if do_modes=False)
        #    the transformation from small displacements in reduced coordinates
        #    to small displacements in Cartesian coordinates. (see Transform class)
        #
        # For the implementation of certain treatments, it is easier to produce
        # a mass-weighted small Hessian immediately. In such cases, the
        # transform is also readily mass-weighted and mass_matrix_small is
        # None.
        raise NotImplementedError


class Full(Treatment):
    """A full vibrational analysis, without transforming to a new set of
       coordinates.
    """
    def __init__(self, svd_threshold=1e-5):
        self.svd_threshold = svd_threshold
        Treatment.__init__(self)

    def compute_zeros(self, molecule, do_modes):
        # An unambigous way to define the 'external' degrees of freedom is as
        # follows: first construct an external basis of the entire systems,
        U, W, Vt = numpy.linalg.svd(molecule.external_basis, full_matrices=False)
        rank = (abs(W) > abs(W[0])*self.svd_threshold).sum()
        self.num_zeros = rank
        if do_modes:
            self.external_basis = Vt[:rank]

    def compute_hessian(self, molecule, do_modes):
        self.hessian_small = molecule.hessian
        self.mass_matrix_small = MassMatrix(molecule.masses3)
        if do_modes:
            atom_division = AtomDivision([], numpy.arange(molecule.size), [])
            self.transform = Transform(None, atom_division)


class ConstrainExt(Treatment):
    """Almost a full vibrational analysis, but with constrained external degrees
       of freedom.

       Note that the current implementation only works correctly when the
       gradient is zero.
    """

    def __init__(self, gradient_threshold=1e-4, svd_threshold=1e-5):
        """Initialize the GassPhase treatment.

           One optional argument:
             gradient_threshold  --  The maximum allowed value of the components
                                     of the Cartesian gradient in atomic units.
                                     When the threshold is exceeded, a
                                     ValueError is raised. [default=1-e4]
        """
        self.gradient_threshold = gradient_threshold
        self.svd_threshold = svd_threshold
        Treatment.__init__(self)

    def compute_zeros(self, molecule, do_modes):
        self.num_zeros = 0

    def compute_hessian(self, molecule, do_modes):
        if abs(molecule.gradient).max() > self.gradient_threshold:
            raise ValueError(
                "Some components of the gradient exceed the threshold "
                "(%.1e > %.1e). The current implementation of the ConstrainExt "
                "treatment only works on optimized geometries." % (
                    abs(molecule.gradient).max(), self.gradient_threshold
                )
            )
        # project the hessian on the orthogonal complement of the basis of small
        # displacements in the external degrees of freedom.
        U, W, Vt = numpy.linalg.svd(molecule.external_basis, full_matrices=True)
        rank = (W/W[0] > self.svd_threshold).sum()
        internal_basis_mw = (Vt[rank:]/numpy.sqrt(molecule.masses3)).transpose()
        # the following hessian is already mass-weighted;
        self.hessian_small = numpy.dot(internal_basis_mw.transpose(), numpy.dot(molecule.hessian, internal_basis_mw))
        # we do not define mass_matrix_small since it is useless when the hessian
        # is already mass-weighted
        if do_modes:
            # also mass-weighted because the hessian is mass-weighted too:
            self.transform = Transform(internal_basis_mw)


class PHVA(Treatment):
    def __init__(self, fixed, svd_threshold=1e-5):
        """Initialize the PHVA treatment.

           One argument:
             fixed  --  a list with fixed atoms, counting starts from zero.
        """
        # QA:
        if len(fixed) == 0:
            raise ValueError("At least one fixed atom is required.")
        # Rest of init:
        self.fixed = numpy.array(fixed)
        self.fixed.sort()
        self.svd_threshold = svd_threshold
        Treatment.__init__(self)

    def compute_zeros(self, molecule, do_modes):
        # This is a bit tricky. Most of the times the number of zero eigenvalues
        # is zero, but there are a few exceptions. When there is one fixed
        # point, there are in general three zeros. When there are two (or more
        # colinear fixed atoms), there is in general one zero. When both the
        # fixed and the free atoms are colinear, there are no zeros.
        #
        # An unambigous way to define the 'external' degrees of freedom is as
        # follows: first construct an external basis of the entire systems,
        U, W, Vt = numpy.linalg.svd(molecule.external_basis, full_matrices=False)
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
        # the nullity of the system, i.e. the number of zero singual values.
        U, W, Vt = numpy.linalg.svd(system, full_matrices=False)
        self.num_zeros = (abs(W) < abs(W[0])*self.svd_threshold).sum()
        if do_modes and self.num_zeros > 0:
            # TODO: fix this
            #self.exernal_basis = ...
            raise NotImplementedError

    def compute_hessian(self, molecule, do_modes):
        free = numpy.zeros(molecule.size - len(self.fixed), int)
        free3 = numpy.zeros(len(free)*3, int)
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
    def __init__(self, subs, svd_threshold=1e-5):
        """Initialize the VSA treatment.

           Frequencies and modes are calculated with the VSA approach:
           Vibrational Subsystem Analysis
           - Zheng and Brooks, ... (2006)
           - Woodcock, ... (2008)

           One argument:
             subs  --  a list with the subsystem atoms, counting starts from zero.
        """
        # QA:
        if len(subs) == 0:
            raise ValueError("At least one subsystem atom is required.")
        # Rest of init:
        self.subs = numpy.array(subs)
        #self.subs.sort()
        self.svd_threshold = svd_threshold
        Treatment.__init__(self)

    def compute_zeros(self, molecule, do_modes):
        # Number of zeros for VSA:
        # 6 zeros if atoms of subsystem are non-collinear
        # 5 zeros if atoms of subsystem are collinear, e.g. when the subsystem contains only 2 atoms
        # 3 zeros if subsystem contains 1 atom
        #
        # method:
        # Construct a kind of inertia matrix (without inertia) of the subsystem
        #           A = transrot_subs . transrot_sub**T
        # Diagonalize A. The rank is the number of expected zero freqs.

        # information subsystem/environment atoms/coordinates:
        subs  = self.subs.tolist()
        subs3 = sum([[3*at, 3*at+1, 3*at+2] for at in xrange(molecule.size) if at in subs],[])
        masses3_subs = molecule.masses3[subs3]

        # transrot_subs contains the 6 global translations and rotations of the subsystem
        # dimension trabsrot: 6 x 3*molecule.size  =>  select subs atoms (6 x 3*len(self.subs))
        transrot_subs = numpy.take(molecule.external_basis,subs3,1)
        A    = numpy.dot( transrot_subs, transrot_subs.transpose() )
        eigv = numpy.linalg.eigvalsh(A)
        rank = (abs(eigv) > abs(eigv[-1])*self.svd_threshold).sum()
        self.num_zeros = rank

        if do_modes and self.num_zeros > 0:
            if self.num_zeros == 3:
                self.external_basis = molecule.external_basis[:3,:]
            elif self.num_zeros == 5:
                # TODO: fix this
                # select 3 translations + 2 rotations
                #self.exernal_basis = external_basis
                raise NotImplementedError
            elif self.num_zeros == 6:
                self.external_basis = molecule.external_basis
            else:
                raise ValueError("Number of zeros is expected to be 3, 5 or 6, but found %i." % self.num_zeros)


    def compute_hessian(self, molecule, do_modes):

        # fill lists with subsystem/environment atoms/coordinates
        subs = self.subs.tolist()
        envi = sum([[at] for at in xrange(molecule.size) if at not in subs],[])
        subs3 = sum([[3*at, 3*at+1, 3*at+2] for at in xrange(molecule.size) if at in subs],[])
        envi3 = sum([[3*at, 3*at+1, 3*at+2] for at in xrange(molecule.size) if at not in subs],[])

        # 1. Construct Hessian (small: 3Nsubs x 3Nsubs)
        # construct H_ss, H_ee, H_es
        hessian_ss = numpy.take(numpy.take(molecule.hessian,subs3,0), subs3,1)
        hessian_ee = numpy.take(numpy.take(molecule.hessian,envi3,0), envi3,1)
        hessian_es = numpy.take(numpy.take(molecule.hessian,envi3,0), subs3,1)
        # construct H_ee**-1 and H_ee**-1 . H_es
        hessian_e1 = numpy.linalg.inv(hessian_ee)
        hessian_e1_es = numpy.dot(hessian_e1,hessian_es)
        # construct H_ss - H_se . H_ee**-1 . H_es
        self.hessian_small = hessian_ss - numpy.dot( hessian_es.transpose(), hessian_e1_es)

        # 2. Construct mass matrix (small: 3Nsubs x 3Nsubs)
        # with corrected mass matrix
        masses3_subs = molecule.masses3[subs3]               # masses subsystem
        masses3_envi = molecule.masses3[envi3]               # masses environment
        tempmat = numpy.zeros((len(envi3),len(subs3)),float) # temporary matrix

        # construct   M_e . H_ee**-1 . H_es
        tempmat = masses3_envi.reshape((-1,1))*hessian_e1_es
        # construct   H_se . H_ee**-1 . M_e . H_ee**-1 . H_es
        massmatrixsmall = numpy.dot(hessian_e1_es.transpose(), tempmat)
        # construct   M_s + H_se . H_ee**-1 . M_e . H_ee**-1 . H_es
        # by adding the diagonal contributions
        massmatrixsmall.ravel()[::len(massmatrixsmall)+1] += masses3_subs
        self.mass_matrix_small = MassMatrix( massmatrixsmall )

        if do_modes:
            atom_division = AtomDivision(envi+subs,[],[])
            self.transform = Transform( numpy.concatenate( (- hessian_e1_es, numpy.identity(len(subs3))),0), atom_division)


class VSA_no_mass(Treatment):
    def __init__(self, subs, svd_threshold=1e-5):
        """Initialize the VSA treatment.

           Frequencies and modes are calculated with the VSA approach:
           Vibrational Subsystem Analysis
           - Zheng and Brooks, ... (2006)
           - Woodcock, ... (2008)

           VSA is performed according to the original version of 2006:
           no mass correction for the environment is included.
           The version of VSA corresponds to the approximation
           of zero mass for all environment atoms.

           One argument:
             subs  --  a list with the subsystem atoms, counting starts from zero.
        """
        # QA:
        if len(subs) == 0:
            raise ValueError("At least one subsystem atom is required.")
        # Rest of init:
        self.subs = numpy.array(subs)
        #self.subs.sort()
        self.svd_threshold = svd_threshold
        Treatment.__init__(self)

    def compute_zeros(self, molecule, do_modes):
        # Number of zeros for VSA:
        # 6 zeros if atoms of subsystem are non-collinear
        # 5 zeros if atoms of subsystem are collinear, e.g. when the subsystem contains only 2 atoms
        # 3 zeros if subsystem contains 1 atom
        #
        # method:
        # Construct a kind of inertia matrix (without inertia) of the subsystem
        #           A = transrot_subs . transrot_sub**T
        # Diagonalize A. The rank is the number of expected zero freqs.

        # information subsystem/environment atoms/coordinates:
        subs  = self.subs.tolist()
        subs3 = sum([[3*at, 3*at+1, 3*at+2] for at in xrange(molecule.size) if at in subs],[])
        masses3_subs = molecule.masses3[subs3]

        # transrot_subs contains the 6 global translations and rotations of the subsystem
        # dimension trabsrot: 6 x 3*molecule.size  =>  select subs atoms (6 x 3*len(self.subs))
        transrot_subs = numpy.take(molecule.external_basis,subs3,1)
        A    = numpy.dot( transrot_subs, transrot_subs.transpose() )
        eigv = numpy.linalg.eigvalsh(A)
        rank = (abs(eigv) > abs(eigv[-1])*self.svd_threshold).sum()
        self.num_zeros = rank

        if do_modes and self.num_zeros > 0:
            if self.num_zeros == 3:
                self.external_basis = molecule.external_basis[:3,:]
            elif self.num_zeros == 5:
                # TODO: fix this
                # select 3 translations + 2 rotations
                #self.exernal_basis = external_basis
                raise NotImplementedError
            elif self.num_zeros == 6:
                self.external_basis = molecule.external_basis
            else:
                raise ValueError("Number of zeros is expected to be 3, 5 or 6, but found %i." % self.num_zeros)


    def compute_hessian(self, molecule, do_modes):

        # fill lists with subsystem/environment atoms/coordinates
        subs = self.subs.tolist()
        envi = sum([[at] for at in xrange(molecule.size) if at not in subs],[])
        subs3 = sum([[3*at, 3*at+1, 3*at+2] for at in xrange(molecule.size) if at in subs],[])
        envi3 = sum([[3*at, 3*at+1, 3*at+2] for at in xrange(molecule.size) if at not in subs],[])

        # 1. Construct Hessian (small: 3Nsubs x 3Nsubs)
        # construct H_ss, H_ee, H_es
        hessian_ss = numpy.take(numpy.take(molecule.hessian,subs3,0), subs3,1)
        hessian_ee = numpy.take(numpy.take(molecule.hessian,envi3,0), envi3,1)
        hessian_es = numpy.take(numpy.take(molecule.hessian,envi3,0), subs3,1)
        # construct H_ee**-1 and H_ee**-1 . H_es
        hessian_e1 = numpy.linalg.inv(hessian_ee)
        hessian_e1_es = numpy.dot(hessian_e1,hessian_es)
        # construct H_ss - H_se . H_ee**-1 . H_es
        self.hessian_small = hessian_ss - numpy.dot( hessian_es.transpose(), hessian_e1_es)

        # 2. Construct mass matrix (small: 3Nsubs x 3Nsubs)
        # with plain submatrix M_s
        self.mass_matrix_small = MassMatrix( numpy.diag(numpy.take(molecule.masses3,subs3)) )

        if do_modes:
            atom_division = AtomDivision(envi+subs,[],[])
            self.transform = Transform( numpy.concatenate( (- hessian_e1_es, numpy.identity(len(subs3))),0), atom_division)


