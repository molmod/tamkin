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


from tamkin.geom import transrot_basis, rank_linearity
from tamkin.io.internal import load_chk, dump_chk

import numpy


__all__ = [
    "NMA", "AtomDivision", "Transform", "MassMatrix", "Treatment",
    "Full", "ConstrainExt", "PHVA", "VSA","VSANoMass", "MBH",
    "PHVA_MBH", "Constrain",
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
        self.masses3 = molecule.masses3
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
            keys = [
                "freqs", "mass", "masses3", "inertia_tensor", "multiplicity",
                "symmetry_number", "periodic", "energy", "zeros"
            ]
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
            "freqs", "modes", "mass", "masses", "masses3", "numbers",
            "coordinates", "inertia_tensor", "multiplicity", "symmetry_number",
            "periodic", "energy", "zeros",
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
             | atom_division -- set AtomDivision class
             | matrix -- the linear transformation from the transformed
                         displacements to Cartesian coordinates.

           Attributes:
             |  matrix  --  see above
             |  scalars  --  ...
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
             | mass_block -- the mass matrix associated with the transformed
                             coordinates

           Arguments, if one is given and it is a one-dimensional matrix:
             | mass_diag -- the diagonal of the mass matrix associated with the
                            free atoms (each mass appears three times)

           Arguments, if two are given:  ! Attention for order of arguments.
             | mass_block -- the mass matrix associated with the transformed
                            coordinates
             | mass_diag -- the diagonal of the mass matrix associated with the
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
            self.mass_block = args[0]  #mass_block is first argument
            self.mass_diag  = args[1]  #mass_diag is second argument
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
             | gradient_threshold  --  The maximum allowed value of the components
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
             |  fixed  --  a list with fixed atoms, counting starts from zero.
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

           Frequencies and modes are computed with the VSA approach:
           Vibrational Subsystem Analysis
           - Zheng and Brooks, ... (2006) TODO add reference
           - Woodcock, ... (2008)

           One argument:
             | subs  --  a list with the subsystem atoms, counting starts from zero.
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
        #- If nonperiodic system:
        # 6 zeros if atoms of subsystem are non-collinear
        # 5 zeros if atoms of subsystem are collinear, e.g. when the subsystem contains only 2 atoms
        # 3 zeros if subsystem contains 1 atom
        #- If periodic system:
        # 3 zeros in all cases
        if not molecule.periodic:
            self.num_zeros = rank_linearity(numpy.take(molecule.coordinates,self.subs,0), svd_threshold = self.svd_threshold)
        else:
            self.num_zeros = 3
        if do_modes and self.num_zeros > 0:
            if self.num_zeros == 3:
                self.external_basis = molecule.external_basis[:3,:]  # three translations
            elif self.num_zeros == 5:
                # Compute direction of the linear SUBSystem (with two atoms) and check for highest alignment with one of the axes.
                diff = molecule.coordinates[self.subs[0]] - molecule.coordinates[self.subs[1]]
                axis = 3+abs(diff).tolist().index(max(abs(diff)))   # axis to be skipped
                alphas = [i for i in range(6) if i is not axis]
                self.external_basis = numpy.take(molecule.external_basis,alphas,0)
            elif self.num_zeros == 6:
                self.external_basis = molecule.external_basis
            else:
                raise ValueError("Number of zeros is expected to be 3, 5 or 6, but found %i." % self.num_zeros)

    def compute_hessian(self, molecule, do_modes):

        # fill lists with subsystem/environment atoms/coordinates
        subs = self.subs.tolist()
        envi = sum([[at] for at in xrange(molecule.size) if at not in subs],[])
        subs3 = sum([[3*at, 3*at+1, 3*at+2] for at in subs],[])
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


class VSANoMass(Treatment):
    def __init__(self, subs, svd_threshold=1e-5):
        """Initialize the VSA treatment.

           Frequencies and modes are computed with the VSA approach:
           Vibrational Subsystem Analysis
           - Zheng and Brooks, ... (2006)  TODO ADD REFERENCE
           - Woodcock, ... (2008)

           VSA is performed according to the original version of 2006:
           no mass correction for the environment is included.
           The version of VSA corresponds to the approximation
           of zero mass for all environment atoms.

           One argument:
             | subs  --  a list with the subsystem atoms, counting starts from zero.
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
        # Number of zeros for VSANoMass:
        #- If nonperiodic system:
        # 6 zeros if atoms of subsystem are non-collinear
        # 5 zeros if atoms of subsystem are collinear, e.g. when the subsystem contains only 2 atoms
        # 3 zeros if subsystem contains 1 atom
        #- If periodic system:
        # 3 zeros in all cases
        if not molecule.periodic:
            self.num_zeros = rank_linearity(numpy.take(molecule.coordinates,self.subs,0), svd_threshold = self.svd_threshold)
        else:
            self.num_zeros = 3
        if do_modes and self.num_zeros > 0:
            if self.num_zeros == 3:
                self.external_basis = molecule.external_basis[:3,:]  # three translations
            elif self.num_zeros == 5:
                # Compute direction of the linear SUBSystem (with two atoms) and check for highest alignment with one of the axes.
                diff = molecule.coordinates[self.subs[0]] - molecule.coordinates[self.subs[1]]
                axis = 3+abs(diff).tolist().index(max(abs(diff)))   # axis to be skipped
                alphas = [i for i in range(6) if i is not axis]
                self.external_basis = numpy.take(molecule.external_basis,alphas,0)
            elif self.num_zeros == 6:
                self.external_basis = molecule.external_basis
            else:
                raise ValueError("Number of zeros is expected to be 3, 5 or 6, but found %i." % self.num_zeros)

    def compute_hessian(self, molecule, do_modes):

        # fill lists with subsystem/environment atoms/coordinates
        subs = self.subs.tolist()
        envi = sum([[at] for at in xrange(molecule.size) if at not in subs],[])
        subs3 = sum([[3*at, 3*at+1, 3*at+2] for at in subs],[])
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


class MBH(Treatment):
    def __init__(self, blocks, do_gradient_correction=True, svd_threshold=1e-5):
        """Initialize the MBH treatment.

           Frequencies and modes are computed with the MBH approach:
           Mobile Block Hessian method
           -  J. Chem. Phys. 126 (22): Art. No. 224102, 2007
           -  J. Chem. Phys. 127 (16), Art. No. 164108, 2007
           Mobile Block Hessian method with linked blocks
           -  J. Chem. Theory Comput. 5 (5), 1203-1215, 2009
           -  J. Chem. Phys. 130 (18), Art. No. 084107, 2009

           MBH is ...  # TODO fill in

           One argument:
             | blocks  --  a list of blocks, each block is a list of atoms,
                           counting starts from zero.
           Optional arguments:
             | do_gradient_correction  --  logical, whether gradient correction
                                           to MBH should be added
             | svd_threshold  --  threshold for zero singular values in svd
        """
        # QA:
        if len(blocks) == 0:
            raise ValueError("At least one block is required.")
        # Rest of init:
        self.blocks = blocks
        self.svd_threshold = svd_threshold
        if type(do_gradient_correction).__name__ == 'bool':
            self.do_gradient_correction = do_gradient_correction
        else: raise TypeError("Optional argument do_grad_correction should be boolean.")
        Treatment.__init__(self)

    def compute_zeros(self, molecule, do_modes):
        # Number of zeros for MBH:
        #- If nonperiodic system:
        # 6 zeros if atoms of system are non-collinear
        # 5 zeros if atoms of system are collinear, e.g. when the subsystem contains only 2 atoms
        # 3 zeros if system contains just 1 atom
        #- If periodic system:
        # 3 zeros in all cases
        if not molecule.periodic:
            self.num_zeros = rank_linearity(molecule.coordinates, svd_threshold = self.svd_threshold)
        else:
            self.num_zeros = 3
        if do_modes and self.num_zeros > 0:
            if self.num_zeros == 3:
                self.external_basis = molecule.external_basis[:3,:]
            elif self.num_zeros == 5:
                # Compute direction of the linear SYSTEM (with two atoms) and check for highest alignment with one of the axes.
                diff = molecule.coordinates[0,:] - molecule.coordinates[1,:]
                axis = 3+abs(diff).tolist().index(max(abs(diff)))   # axis to be skipped
                alphas = [i for i in range(6) if i is not axis]
                self.external_basis = numpy.take(molecule.external_basis,alphas,0)
            elif self.num_zeros == 6:
                self.external_basis = molecule.external_basis
            else:
                raise ValueError("Number of zeros is expected to be 3, 5 or 6, but found %i." % self.num_zeros)

    def compute_hessian(self, molecule, do_modes):
        if do_modes:
            self.hessian_small,self.mass_matrix_small,self.transform = \
                     self.compute_matrices_small(molecule,do_modes)
        else:
            self.hessian_small,self.mass_matrix_small = \
                     self.compute_matrices_small(molecule,do_modes)

    def compute_matrices_small(self, molecule, do_modes):
        # Notation: b,b0,b1   --  a block index
        #           block  --  a list of atoms, e.g. [at1,at4,at6]
        #           alphas  --  the 6 block parameter indices (or 5 for linear block)

        # Block information
        blkinfo = Blocks(self.blocks, molecule, self.svd_threshold)
        mbhdim1 = 6*blkinfo.nb_nlin + 5*blkinfo.nb_lin + 3*len(blkinfo.free)

        # TRANSFORM from CARTESIAN to BLOCK PARAMETERS
        U = self.construct_U(molecule,mbhdim1,blkinfo)

        # Construct Hessian in block parameters: Hp = U**T . H . U + correction
        Hp = numpy.dot(numpy.dot( U.transpose(), molecule.hessian) , U)

        # gradient correction
        if self.do_gradient_correction:
         for b,block in enumerate(blkinfo.blocks_nlin_strict+blkinfo.blocks_lin_strict):   # Nonlinear AND linear blocks
            GP = numpy.zeros((3,3),float)
            G = numpy.take(molecule.gradient,block,0)
            P = numpy.take(molecule.coordinates,block,0)
            for i in range(len(block)):
                GP += numpy.dot(G[i,:].reshape((3,1)),P[i,:].reshape((1,3)))
            # note: GP is not symmetric
            p = GP[0,0] # Gx.x
            q = GP[1,1] # Gy.y
            r = GP[2,2] # Gz.z
            s = GP[1,0] # Gy.x
            t = GP[2,0] # Gz.x
            u = GP[2,1] # Gz.y

            corr    = numpy.zeros((6,6),float)
            corr[3:,3:] = numpy.diag([-q-r, -p-r, -p-q])
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
                Hp[col:(col+dim),col:(col+dim)] += numpy.take(numpy.take(corr,alphas,0),alphas,1)

        # Construct mass matrix in block parameters: Mp = U**T . M . U
        Mp = numpy.dot(U.transpose(),  U * molecule.masses3.reshape((-1,1)))

        if blkinfo.is_linked:
            # SECOND TRANSFORM: from BLOCK PARAMETERS to Y VARIABLES
            # Necessary if blocks are linked to each other.
            nullspace = self.construct_nullspace_K(molecule,mbhdim1,blkinfo)

            My = numpy.dot(nullspace.transpose(), numpy.dot( Mp,nullspace) )
            Hy = numpy.dot(nullspace.transpose(), numpy.dot( Hp,nullspace) )

            # TODO
            # gradient correction of the second transform...

        if do_modes:
            if not blkinfo.is_linked:
                return Hp, MassMatrix(Mp), Transform(U)
            else:
                return Hy, MassMatrix(My), Transform(numpy.dot(U, nullspace))
        else:
            if not blkinfo.is_linked:
                return Hp, MassMatrix(Mp)
            else:
                return Hy, MassMatrix(My)

    def construct_U(self,molecule,mbhdim1,blkinfo):
        # Construct first transformation matrix
        D = transrot_basis(molecule.coordinates)   # is NOT mass-weighted

        U = numpy.zeros((3*molecule.size, mbhdim1),float)

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

    def construct_nullspace_K(self,molecule,mbhdim1,blkinfo):
        # SECOND TRANSFORM: from BLOCK PARAMETERS to Y VARIABLES
        # Necessary if blocks are linked to each other.
            # Construct K matrix, with constraints
            D = transrot_basis(molecule.coordinates)   # is NOT mass-weighted
            nbrows = (numpy.sum(blkinfo.sharenbs)-molecule.size)*3
            K = numpy.zeros(( nbrows, mbhdim1-3*len(blkinfo.free)), float)
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
                        D0 = numpy.take(D[:,3*at:3*(at+1)],alphas,0)

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
                            D1 = numpy.take(D[:,3*at:3*(at+1)],alphas,0)

                        K[row:row+3,sta1:end1] = -D1.transpose()
                        row += 3

            # Do SVD of matrix K
            u,s,vh = numpy.linalg.svd(K)

            # construct nullspace of K
            rank = sum(s>max(s)*self.svd_threshold)
            nullspace = vh[rank:,:].transpose()
            [r_null,c_null] = nullspace.shape
            n = numpy.zeros((mbhdim1,c_null+3*len(blkinfo.free)),float)
            n[:r_null,:c_null] = nullspace
            n[r_null:,c_null:] = numpy.identity(3*len(blkinfo.free),float)
            return n


class Blocks(object):
    def __init__(self,blocks,molecule,svd_threshold):
        """
        initialization:  Blocks(blocks,N) with

        Arguments:
          | blocks --   a list of lists of atoms
                       [ [at1,at5,at3], [at4,at5], ...]
                       with a list of atoms for each block
          | molecule -- Molecule object, necessary for N (total nb
                        of atoms) and positions (linearity of blocks).
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
        #remove single atoms and empty blocks
        for i in range(len(to_remove)):
            del blocks[to_remove[len(to_remove)-i-1]]   #remove starting from largests b

        #checking
        fixed = set(sum(blocks,[]))
        free  = [atom for atom in range(N) if atom not in fixed] #list of N_E integers or empty list
        fixed = [atom for atom in range(N) if atom in fixed]

        # check for linearity and fill in dimensions
        #D = molecule.external_basis
        #dim_block=numpy.zeros((len(blocks)),int)
        indices_blocks_nlin = []    # nonlinear blocks
        indices_blocks_lin  = []    # linear blocks
        for b,block in enumerate(blocks):
            rank = rank_linearity(numpy.take(molecule.coordinates,block,0), svd_threshold=svd_threshold)
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
        bA1 = numpy.zeros((molecule.size),int)
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
        skip_axis_lin = numpy.zeros((nb_lin),int)
        for b,block in enumerate(blocks_lin):       # do not use strict partition here
            diff = molecule.coordinates[block[0]] - molecule.coordinates[block[1]]
            skip_axis_lin[b] = 3+abs(diff).tolist().index(max(abs(diff)))   # axis to be skipped


        # Check if there are linked blocks
        sharenbs = numpy.zeros((molecule.size),int)  # share number of each atom
        for (at,apps) in appearances.iteritems():
            sharenbs[at] = len(apps)
        is_linked = False
        for sharenb in numpy.ravel(sharenbs):
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
    def __init__(self, fixed, blocks, svd_threshold=1e-5):
        """Initialize the PHVA_MBH treatment.

           Two arguments:
             | fixed  --  a list with fixed atoms, counting starts from zero.
             | blocks  --  a list of blocks, each block is a list of atoms
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
        self.fixed = numpy.array(fixed)
        MBH.__init__(self, blocks, svd_threshold=svd_threshold)

    def compute_zeros(self, molecule, do_modes):
        # [ See explanation PHVA ]
        U, W, Vt = numpy.linalg.svd(molecule.external_basis, full_matrices=False)
        rank = (abs(W) > abs(W[0])*self.svd_threshold).sum()
        external_basis = Vt[:rank]
        fixed3 = []
        for i in self.fixed:
            fixed3.append(3*i)
            fixed3.append(3*i+1)
            fixed3.append(3*i+2)
        system = external_basis[:,fixed3].transpose()
        U, W, Vt = numpy.linalg.svd(system, full_matrices=False)
        self.num_zeros = (abs(W) < abs(W[0])*self.svd_threshold).sum()
        if do_modes and self.num_zeros > 0:
            # TODO: fix this
            #self.exernal_basis = ...
            raise NotImplementedError

    def compute_hessian(self,molecule,do_modes):
        # Make submolecule
        selectedatoms = [at for at in xrange(molecule.size) if at not in self.fixed]
        selectedcoords = sum([[3*at,3*at+1,3*at+2] for at in selectedatoms],[])

        from tamkin.data import Molecule
        submolecule = Molecule(
            numpy.take(molecule.numbers, selectedatoms),
            numpy.take(molecule.coordinates, selectedatoms, 0),
            numpy.take(molecule.masses, selectedatoms),
            molecule.energy,
            numpy.take(molecule.gradient,selectedatoms,0),
            numpy.take(numpy.take(molecule.hessian,selectedcoords,0),selectedcoords,1),
            molecule.multiplicity,
            0, # undefined molecule.symmetry_number
            False # molecule.is_periodic
        )

        # adapt numbering in blocks
        shifts = numpy.zeros((molecule.size),int)
        for fixat in self.fixed:
            shifts[fixat:] = shifts[fixat:]+1
        for bl,block in enumerate(self.blocks):
            for at,atom in enumerate(block):
                self.blocks[bl][at] = atom - shifts[atom]

        if do_modes:
            self.hessian_small, self.mass_matrix_small, transform = self.compute_matrices_small(submolecule, do_modes)
            transf = numpy.zeros((3*molecule.size, transform.matrix.shape[1]),float)
            transf[selectedcoords,:] = transform.matrix
            self.transform = Transform(transf)
        else:
            self.hessian_small, self.mass_matrix_small = self.compute_matrices_small(submolecule, do_modes)


class Constrain(Treatment):
    def __init__(self, constraints, do_grad_correction=True, svd_threshold=1e-5):
        """Initialize the Constrain treatment.

           One argument:
             | constraints  --  a list with constraints of internal coordinates:
                                [at1,at2] to constrain a distance,
                                [at1,at2,at3] to constrain an angle,
                                [at1,at2,at3,at4] to constrain a dihedral angle.
           Optional:
             | do_grad_correction  --  whether gradient correction should be applied
             | svd_threshold  --  threshold for singular value decomposition
        """
        # QA:
        if len(constraints) == 0:
            raise ValueError("At least one constraint is required.")
        # Rest of init:
        self.constraints = constraints
        self.do_grad_correction = do_grad_correction
        self.svd_threshold = svd_threshold
        Treatment.__init__(self)

    def compute_zeros(self, molecule, do_modes):
        # Number of zeros for Constrain:
        #- If nonperiodic system:
        # 6 zeros if atoms of subsystem are non-collinear
        # 5 zeros if atoms of subsystem are collinear, e.g. when the subsystem contains only 2 atoms
        #- If periodic system:
        # 3 zeros in all cases

        if not molecule.periodic:
            #self.num_zeros = rank_linearity(numpy.take(molecule.coordinates,self.subs,0), svd_threshold = self.svd_threshold)
            U, W, Vt = numpy.linalg.svd(molecule.external_basis, full_matrices=False)
            rank = (abs(W) > abs(W[0])*self.svd_threshold).sum()
            self.num_zeros = rank
            if do_modes:
                self.external_basis = Vt[:rank]
        else:
            self.num_zeros = 3
            if do_modes:
                self.external_basis = self.external_basis = molecule.external_basis[:3,:]  # three translations

    def compute_hessian(self, molecule, do_modes):

        # make constraint matrix
        constrmat = numpy.zeros((3*molecule.size,len(self.constraints)))
        count = 0
        for i,constraint in enumerate(self.constraints):
            # constrain a distance
            if len(constraint) == 2:
                at1 = constraint[0]
                at2 = constraint[1]
                dist = numpy.sqrt(numpy.sum((molecule.coordinates[at1] - molecule.coordinates[at2])**2))
                dist = 1.0
                constrmat[3*at1:3*at1+3, count] = (molecule.coordinates[at1] - molecule.coordinates[at2])/dist
                constrmat[3*at2:3*at2+3, count] = (molecule.coordinates[at2] - molecule.coordinates[at1])/dist
                count += 1
        # determine the orthogonal complement of the basis of small
        # displacements determined by the constraints.
        U, W, Vt = numpy.linalg.svd(constrmat.transpose(), full_matrices=True)
        rank = (W/W[0] > self.svd_threshold).sum()
        nullspace = Vt[rank:].transpose()

        # mass matrix small = nullspace^T . M . nullspace
        # hessian     small = nullspace^T . H . nullspace + gradient correction
        self.mass_matrix_small = MassMatrix(numpy.dot(nullspace.transpose(),nullspace * molecule.masses3.reshape((-1,1))) )
        self.hessian_small     = numpy.dot(nullspace.transpose(), numpy.dot(molecule.hessian, nullspace))

        # check if gradient is small enough in this complement: overlap with nullspace should be small enough
        # print  numpy.sum(numpy.dot(nullspace.transpose(), numpy.ravel(molecule.gradient))**2)

        if self.do_grad_correction:
            self.do_the_gradient_correction(constrmat.transpose(), U,W,Vt, rank, nullspace, molecule.gradient)

        if do_modes:
            self.transform = Transform(nullspace)


    def do_the_gradient_correction(self, K,u,s,vh, rank, nullspace, gradient):


#         print "Cart grad", self.gradient
#         print "W grad", G_W
#         print "orthogonal?",numpy.dot(nullspace.transpose(),G_W)
#         print "G_W, sum", sum(abs(G_W))
#         print "grad, sum", sum(abs(self.gradient.ravel()))


        # GRADIENT CORRECTION
        # construct generalized inverse L of K
        # zero rows for free atoms not included in L

        # construct generalized inverse
        L = numpy.zeros((K.shape[1],K.shape[0]),float)
        for i,sigma in enumerate(s[:rank]):
           L[i,:] = 1/sigma*u[:,i]
        L = numpy.dot(vh.transpose(),L)

        # construct right hand side of the equation
        Y = numpy.zeros((K.shape[0],len(self.hessian_small)**2))
        y = numpy.zeros((K.shape[0],1))
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
        X = -numpy.dot(L, Y)

        # apply correction
        count = 0
        for i in xrange(len(self.hessian_small)):
            for j in xrange(len(self.hessian_small)):

                self.hessian_small[i,j] += numpy.sum( numpy.ravel(gradient)*X[:,count] )
                count += 1


