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
"""Implementation of free and hindered rotors

   The current implementation supports one-dimensional free and hindered rotors.
   For practical applications it is apparently not necessary to consider higher-
   dimensional hindered rotors. [1]

   [1] Chemical Physics, Vol. 328 (1-3) 251 - 258, 2006
"""

from tamkin.partf import Info, StatFysTerms, helper_vibrations, \
    helpert_vibrations, helpertt_vibrations, helper_levels, helpert_levels, \
    helpertt_levels
from tamkin.nma import NMA, MBH
from tamkin.geom import transrot_basis

from molmod import deg, kjmol, angstrom, centimeter, amu, boltzmann, \
    lightspeed, cached

import numpy as np


__all__ = [
    "HarmonicBasis", "RotorError", "compute_cancel_frequency_mbh",
    "compute_moments", "Rotor"
]


class HarmonicBasis(object):
    """A harmonic basis set for periodic one-dimensional QM systems

       In addition to the definition of the basis set, this class also
       implements the kinetic and potential energy operators required for the
       solution of the schrodinger equation. The workflow is as follows:

       >>> a = 10.0                             # the size of the system
       >>> hb = HarmonicBasis(10, a)            # create basis object
       >>> grid = np.arange(0.0, 10.01, 1.0)    # define a grid
       >>> v = -np.exp(-((grid-5)/2)**2)        # define a potential on the grid
       >>> v_coeffs = hb.fit_fn(grid, v, 10)    # expand the potential in the basis
       >>> mass = 1.0
       >>> energies, wfns = hb.solve(mass, v_coeffs, evecs=True) # solve problem
    """
    def __init__(self, nmax, a):
        """
           Arguments:
            | ``nmax`` -- The maximum index used for the gaussian basis. Hence
                          the the number of basis functions is 2*nmax+1
            | ``a`` -- The length of the periodic system


           The basis functions are:
            | 1/sqrt(a/2)
            | cos(2*pi*x/a)/sqrt(a/2)
            | sin(2*pi*x/a)/sqrt(a/2)
            | ...
            | cos(nmax*2*pi*x/a)/sqrt(a/2)
            | sin(nmax*2*pi*x/a)/sqrt(a/2)

           This is an orthonormal basis. This object has methods to construct a
           Hamiltonian operator in this basis and to transform a function on a
           grid into an expansion in this basis (and back).
        """
        self.nmax = nmax
        self.a = a

    size = property(lambda self: 2*self.nmax+1,
        doc="The size of the basis set. (read-only attribute)")

    def get_empty_op(self):
        """Returns an empty operator (zero)"""
        return np.zeros((self.size, self.size), float)

    def get_hamiltonian_op(self, mass, potential):
        """Returns the Hamiltonian operator for the given mass and potential

           Arguments:
            | ``mass`` -- the mass of the particle
            | ``potential`` -- the expansion coefficients of the potential energy
        """
        op = self.get_empty_op()
        self._add_kinetic_op(op, mass)
        self._add_potential_op(op, potential)
        return op

    def _add_kinetic_op(self, op, mass):
        """Add the kinetic energy to the given operator"""
        factor = (0.5/mass*(2*np.pi/self.a)**2)
        tmp = factor*np.arange(1, self.nmax+1)**2
        diag = op.ravel()[::self.size+1]
        diag[1::2] += tmp
        diag[2::2] += tmp

    def _add_potential_op(self, op, potential):
        """Add the potential energy to the given operator"""
        c = np.zeros(self.nmax+1,float)
        s = np.zeros(self.nmax+1,float)
        c[0] = potential[0]/np.sqrt(self.a)
        c[1:] = potential[1::2]/np.sqrt(self.a)
        s[1:] = potential[2::2]/np.sqrt(self.a)
        op[0,0] += c[0]
        for i0 in xrange(1,self.nmax+1):
            op[2*i0-1,0] += c[i0]
            op[2*i0-0,0] += s[i0]
            op[0,2*i0-1] += c[i0]
            op[0,2*i0-0] += s[i0]
        c[1:] /= np.sqrt(2)
        s[1:] /= np.sqrt(2)
        for i0 in xrange(1,self.nmax+1):
            for i1 in xrange(1,self.nmax+1):
                k = i0+i1
                if k <= self.nmax:
                    op[2*i0-1,2*i1-1] += c[k] # cc
                    op[2*i0-0,2*i1-0] -= c[k] # ss
                    op[2*i0-1,2*i1-0] += s[k] # cs
                    op[2*i0-0,2*i1-1] += s[k] # sc
                k = abs(i0-i1)
                if k <= self.nmax:
                    sign = 2*(i0>i1)-1
                    op[2*i0-1,2*i1-1] += c[k] # cc
                    op[2*i0-0,2*i1-0] += c[k] # ss
                    op[2*i0-1,2*i1-0] -= sign*s[k] # cs
                    op[2*i0-0,2*i1-1] += sign*s[k] # sc

    def solve(self, mass, potential, evecs=False):
        """Return the energies and wavefunctions for the given mass and potential

           Arguments:
            | ``mass`` -- the mass of the particle
            | ``potential`` -- the expansion coefficients of the potential energy

           Optional argument:
            | ``evecs`` -- When True, also the eigenstates are returned.
                           [default=False]
        """
        H = self.get_hamiltonian_op(mass, potential)
        if evecs:
            return np.linalg.eigh(H)
        else:
            return np.linalg.eigvalsh(H)

    def eval_fn(self, grid, coeffs):
        """Evaluate the function represented by coeffs

           Arguments:
            | ``grid`` -- the values at which the function must be evaluated
            | ``coeffs`` -- the expansion coefficients
        """
        result = np.zeros(grid.shape, float) + coeffs[0]/np.sqrt(self.a)
        for i in xrange(self.nmax):
            arg = ((i+1)*2*np.pi/self.a)*grid
            result += (coeffs[2*i+1]/np.sqrt(self.a/2))*np.cos(arg)
            result += (coeffs[2*i+2]/np.sqrt(self.a/2))*np.sin(arg)
        return result

    def eval_deriv(self, grid, coeffs):
        """Evaluate the derivative of function represented by coeffs

           Arguments:
            | ``grid`` -- the values at which the derivative must ben evaluated
            | ``coeffs`` -- the expansion coefficients
        """
        result = np.zeros(grid.shape, float)
        for i in xrange(self.nmax):
            scale = ((i+1)*2*np.pi/self.a)
            arg = scale*grid
            result -= (coeffs[2*i+1]*scale/np.sqrt(self.a/2))*np.sin(arg)
            result += (coeffs[2*i+2]*scale/np.sqrt(self.a/2))*np.cos(arg)
        return result

    def eval_deriv2(self, grid, coeffs):
        """Evaluate the second derivative of function represented by coeffs

           Arguments:
            | ``grid`` -- the values at which the second derivative must be
                          evaluated
            | ``coeffs`` -- the expansion coefficients
        """
        result = np.zeros(grid.shape, float)
        for i in xrange(self.nmax):
            scale = ((i+1)*2*np.pi/self.a)
            arg = scale*grid
            result -= (coeffs[2*i+1]*scale**2/np.sqrt(self.a/2))*np.cos(arg)
            result -= (coeffs[2*i+2]*scale**2/np.sqrt(self.a/2))*np.sin(arg)
        return result

    def fit_fn(self, grid, v, dofmax, rotsym=1, even=False, rcond=0.0, v_threshold=0.01):
        """Fit the expansion coefficients that represent function f

           Arguments:
            | ``grid`` -- The x values on which the function f is known.
            | ``v`` -- The function to be represented by expansion coefficients.
            | ``dofmax`` -- The maximum number of cosines in the fit. When
                            even==False, the same number of sines is also
                            included.

           Optional arguments:
            | ``rotsym`` -- Impose this rotational symmetry [default=1]. Must be
                            an integer and is least 1..
            | ``even`` -- Only fit even functions, i.e. cosines. [default=False]
            | ``rcond`` -- The cutoff for the singular values in the least
                           squares fit. [default=1e-10]
            | ``v_threshold`` -- Tolerance on the relative error between the
                                 Fourier expansion and the data points of the
                                 scan. [default=0.01]. Absolute errors smaller
                                 than 1 kJ/mol are always ignored.

           In case the Fourier expansion represents a poor fit (determined by
           v_threshold), a ValueError is raised. It means that you have to
           check your torsional scan datapoints for errors.
        """
        if rotsym < 1:
            raise ValueError("rotsym must be at least 1.")

        # construct the design matrix
        ncos = min(dofmax, self.nmax/rotsym)
        if even:
            A = np.zeros((len(grid), ncos+1), float)
        else:
            A = np.zeros((len(grid), 2*ncos+1), float)
        A[:,0] = 1.0/np.sqrt(self.a)
        counter = 1
        for i in xrange(ncos):
            arg = ((i+1)*rotsym*2*np.pi/self.a)*grid
            A[:,counter] = np.cos(arg)/np.sqrt(self.a/2)
            counter += 1
            if not even:
                A[:,counter] = np.sin(arg)/np.sqrt(self.a/2)
                counter += 1

        coeffs, residuals, rank, S = np.linalg.lstsq(A, v, rcond)

        # check the error
        residual = np.dot(A, coeffs) - v
        rmsd = (residual**2).mean()**0.5
        rms = (v**2).mean()**0.5
        abs_threshold = max(1*kjmol, rms*v_threshold)
        if (rmsd > abs_threshold).any():
            raise ValueError("Residual is too large. (poor Fourier expansion.) rmsd [kJ/mol] = %f, rms [kJ/mol] = %f" % (rmsd/kjmol, rms/kjmol))

        # collect the parameters in a convenient array
        result = np.zeros(self.size)
        result[0] = coeffs[0]
        if even:
            tmp = result[2*rotsym-1::2*rotsym]
            tmp[:ncos] = coeffs[1:]
        else:
            tmp = result[2*rotsym-1::2*rotsym]
            tmp[:ncos] = coeffs[1::2]
            tmp = result[2*rotsym-0::2*rotsym]
            tmp[:ncos] = coeffs[2::2]
        return result


class RotorError(Exception):
    """This exception is raised when :class:`Rotor` specific errors are
       encountered.
    """
    pass


def compute_cancel_frequency_mbh(molecule, dihedral, top_indexes):
    """Compute the frequency of the rotor in the HO approximation

       This function is based on the MBH method and returns the frequency that
       has to be canceled when this mode is replaced by a free or hindered
       rotor.

       Arguments:
        | ``molecule`` -- A Molecule object. (see :mod:`tamkin.data`)
        | ``top_indexes`` -- The indexes of the rotor atoms.
    """
    axis = tuple(dihedral[1:3])
    other_top_indexes = tuple(set(xrange(molecule.size)) - set(top_indexes) - set(axis))
    blocks = [
        axis + tuple(top_indexes),
        axis + tuple(other_top_indexes),
    ]
    nma = NMA(molecule, MBH(blocks))
    if len(nma.freqs) != 7:
        raise RotorError("Expecting 7 frequencies, got %i" % len(nma.freqs))
    non_zero = [i for i in xrange(7) if i not in nma.zeros][0]
    return nma.freqs[non_zero]


def compute_moments(coordinates, masses3, center, axis, indexes):
    """Computes the absolute and the relative moment of an internal rotor

       Arguments:
        | ``coordinates`` -- The coordinates of all atoms, float numpy array
                             with shape (N,3).
        | ``masses3`` -- The diagonal of the mass matrix, each mass is repeated
                         three tines, float numpy array with shape 3N.
        | ``center`` -- A point on the rotation axis. Float numpy array with
                        shape 3.
        | ``axis`` -- A unit vector with the direction of the rotation axis.
                      Float numpy array with shape 3.
        | ``indexes`` -- The indexes of the atoms that belong to the rotor.

       The implementation is based on the transformation of the mass-matrix to
       the internal rotation coordinate. The derivative of the internal
       coordinate towards cartesian coordinates represents a small displacement
       of the atoms. This displacement is constrained to show no external
       linear or rotational moment.
    """
    # the derivative of the cartesian coordinates towards the rotation
    # angle of the top:
    rot_tangent = np.zeros((len(coordinates)*3), float)
    for i in indexes:
        rot_tangent[3*i:3*i+3] = np.cross(coordinates[i]-center, axis)
    # transform this to a derivative without global linear or angular
    # momentum. This is done by projecting on the basis of external degrees
    # of freedom in mass-weighted coordinates, and subsequently subtracting
    # that projection from the original tangent
    basis = transrot_basis(coordinates)
    A = np.dot(basis*masses3, basis.transpose())
    B = np.dot(basis*masses3, rot_tangent)
    alphas = np.linalg.solve(A,B)
    rot_tangent_relative = rot_tangent - np.dot(alphas, basis)
    return (
        (rot_tangent**2*masses3).sum(),
        (rot_tangent_relative**2*masses3).sum()
    )


class Rotor(Info, StatFysTerms):
    """Partition function term for a one-dimensional rotor

       The contribution from the free or hindered rotor to the partition
       function is based on the quantum mechanical solution of the rotational
       motion. To avoid double counting problems, one must also provide the
       frequency of this motion as if it was treated as a harmonic oscillator.
       The corresponding contribution to the partition function is subtracted.
       (Use compute_cancel_frequency to obtain this frequency.)
    """
    def __init__(self, rot_scan, molecule=None, cancel_freq='mbh', suffix=None,
                 rotsym=1, even=False, num_levels=50, dofmax=5,
                 v_threshold=0.01, large_fixed=False):
        """
           Arguments:
            | ``rot_scan`` -- A rotational scan object. (free or hindered rotor
                              information)

           Optional arguments:
            | ``molecule`` -- Molecule to which the rotor applies, is used to
                              compute the cancelation frequency. Not required
                              when the cancel_freq argument is present.
            | ``cancel_freq`` -- The frequency to cancel in the vibrational
                                 partition function. This can also be 'mbh' or
                                 'scan' to indicate that the cancel frequency
                                 should be computed using the MBH method or
                                 based on the second order derivative of the
                                 rotational potential. Note that the latter
                                 option is only possible in the case of the
                                 hindered rotor formalism. [default='mbh']
            | ``suffix`` -- A name suffix used to distinguish between different
                            rotors.
            | ``rotsym`` -- The rotational symmetry of the rotor. [default=1]
            | ``even`` -- True of the rotor is not chiral, i.e. when it has an
                          even potential
            | ``num_levels`` -- The number of energy levels considered in the
                                QM treatment of the rotor [default=50]
            | ``dofmax`` -- The maximum number of cosines used to represent the
                            torsional potential. if the potential is not even,
                            the same number of sines is also used. [default=5]
            | ``v_threshold`` -- Tolerance on the relative error between the
                                 Fourier expansion and the data points of the
                                 scan. [default=0.01]. Absolute errors smaller
                                 than 1 kJ/mol are always ignored.
            | ``large_fixed`` -- When True, assume that the large part of the
                                 system is fixed in space while the small part
                                 rotates. (this means that the absolute moment
                                 of the rotor is used instead of the relative
                                 moment)

           In case the Fourier expansion of the potential represents a poor fit
           (determined by v_threshold), a ValueError is raised. It means that
           you have to check your torsional scan datapoints for errors.
        """
        self.rot_scan = rot_scan
        self.cancel_freq = cancel_freq
        self.cancel_method = 'given by the user'
        # the cancelation frequency
        if self.cancel_freq == 'mbh':
            self.cancel_freq = compute_cancel_frequency_mbh(
                molecule, self.rot_scan.dihedral, self.rot_scan.top_indexes
            )
            self.cancel_method = 'computed with mbh'
        elif self.cancel_freq == 'scan':
            if self.rot_scan.potential is None:
                raise RotorError("The option cancel_freq='scan' can only be used with hindered rotors.")
            self.cancel_method = 'derived from the torsional potential'
            # the actual computation is performed later
        else:
            raise ValueError("Could not interpret cancel_freq=%s" % cancel_freq)
        #
        if suffix is None:
            self.suffix = "_".join(str(i) for i in rot_scan.top_indexes)
        else:
            self.suffix = suffix
        self.rotsym = rotsym
        self.even = even
        self.num_levels = num_levels
        self.dofmax = dofmax
        self.v_threshold = v_threshold
        self.large_fixed = large_fixed
        if rot_scan.potential is None:
            Info.__init__(self, "free_rotor_%s" % self.suffix)
        else:
            Info.__init__(self, "hindered_rotor_%s" % self.suffix)
        StatFysTerms.__init__(self, 2) # two terms

    def init_part_fun(self, nma, partf):
        """See :meth:`tamkin.partf.StatFys.init_part_fun`"""
        if nma.periodic:
            raise NotImplementedError("Rotors in periodic systems are not supported yet")

        self.center = nma.coordinates[self.rot_scan.dihedral[1]]
        self.axis = nma.coordinates[self.rot_scan.dihedral[2]] - self.center
        self.axis /= np.linalg.norm(self.axis)
        self.moment, self.reduced_moment = compute_moments(
            nma.coordinates, nma.masses3, self.center, self.axis, self.rot_scan.top_indexes
        )
        if self.large_fixed:
            moment = self.moment
        else:
            moment = self.reduced_moment
        from molmod.ic import dihed_angle
        self.nma_angle = dihed_angle(nma.coordinates[self.rot_scan.dihedral])[0]
        # the energy levels
        if self.rot_scan.potential is None:
            # free rotor
            self.energy_levels = np.zeros(self.num_levels, float)
            for i in xrange(self.num_levels-1):
                index = i/2+1
                self.energy_levels[i+1] = index**2/(2*moment)
            self.hb = None
            self.v_coeffs = None
            self.v_ref = 0.0
        else:
            # hindered rotor
            self.hb = HarmonicBasis(self.num_levels, 2*np.pi)
            angles, energies = self.potential

            self.v_coeffs = self.hb.fit_fn(angles, energies, self.dofmax,
                self.rotsym, self.even, self.v_threshold)
            self.energy_levels = self.hb.solve(moment, self.v_coeffs)
            self.energy_levels = self.energy_levels[:self.num_levels]

        # the cancelation frequency based on the scan
        if self.hb is None:
            if not isinstance(self.cancel_freq, float):
                raise ValueError('No cancelation frequency was computed for rotor "%s"' % self.name)
        else:
            force_constant = self.hb.eval_deriv2(np.array([self.nma_angle]), self.v_coeffs)[0]
            scan_cancel_freq = np.sqrt(force_constant/moment)/(2*np.pi)
            if self.cancel_freq == 'scan':
                self.cancel_freq = scan_cancel_freq
            elif abs(self.cancel_freq - scan_cancel_freq) > 0.3*abs(self.cancel_freq):
                print 'WARNING: The cancelation frequency of rotor "%s" obtained with MBH (%.1f cm^-1) deviates a lot from the one derived from the scan (%.1f cm^-1).' % (
                    self.name, self.cancel_freq/(lightspeed/centimeter), scan_cancel_freq/(lightspeed/centimeter)
                )

        # print a warning is the cancelation frequency is rather high.
        if self.cancel_freq > 500*(lightspeed/centimeter):
            print 'WARNING: the cancelation frequency of rotor "%s" is rather high: %.1f cm^-1.' % (
                self.name, self.cancel_freq/(lightspeed/centimeter)
            )
        elif self.cancel_freq <= 0:
            print 'WARNING: the cancelation frequency of rotor "%s" is negative: %.1f cm^-1.' % (
                self.name, self.cancel_freq/(lightspeed/centimeter)
            )


        # scaling factors
        self.freq_scaling = partf.vibrational.freq_scaling
        self.zp_scaling = partf.vibrational.zp_scaling
        self.classical = partf.vibrational.classical

    @cached
    def potential(self):
        """A tuple with angles and potential energies (hindered only)

           If the rotor is free, the result is None. The reference for the
           potential energy is the energy of the reference geometry used in the
           partition function.
        """
        if self.rot_scan.potential is None:
            return None
        else:
            a = 2*np.pi
            angles, energies = self.rot_scan.potential.copy()
            # apply periodic boundary conditions
            angles -= np.floor(angles/a)*a
            # set reference energy, which is take to be the energy of the geometry
            # with a dihedral angle closest to the dihedral angle of the reference
            # geometry in the nma object. We can not take the nma energy because
            # the rotational energy barriers may be computed at another level
            # that the nma energy.
            deltas = angles - self.nma_angle
            # apply periodic boundary conditions
            deltas -= np.floor(deltas/a)*a
            # get the right energy
            energies -= energies[deltas.argmin()]
            return angles, energies

    def dump(self, f):
        """Write all the information about the rotor to a file

           This method is part of the PartFun API and should never be called
           directly. It will only work properly once the init_part_fun method
           is called.

           Arguments:
            | ``f`` -- A file-like object.
        """
        Info.dump(self, f)
        # parameters
        print >> f, "    Dihedral indexes: %s" % " ".join(str(i) for i in self.rot_scan.dihedral)
        print >> f, "    Top indexes: %s" % " ".join(str(i) for i in self.rot_scan.top_indexes)
        print >> f, "    Rotational symmetry: %i" % self.rotsym
        print >> f, "    Even potential: %s" % self.even
        if self.rot_scan.potential is None:
            print >> f, "    This is a free rotor"
        else:
            angles, energies = self.rot_scan.potential
            print >> f, "    This is a hindered rotor"
            print >> f, "    Maximum number of cosines in the fit: %i" % self.dofmax
            angles, energies = self.potential
            fit_energies = self.hb.eval_fn(angles, self.v_coeffs)
            rmsd = ((fit_energies - energies)**2).mean()**0.5
            rms = (energies**2).mean()**0.5
            rrmsd = rmsd/rms
            print >> f, "    RMSD of the fit [kJ/mol]: %.2f" % (rmsd/kjmol)
            print >> f, "    Relative RMSD of the fit the fit [%%]: %.1f" % (rrmsd*100)
            print >> f, "    Pearson R^2 of the fit [%%]: %.3f" % ((1-rrmsd**2)*100)

            print >> f, "    Potential: Angle [deg]    Energy [kJ/mol]"
            for i in xrange(len(angles)):
                print >> f, "              % 7.2f         %6.1f" % (angles[i]/deg, energies[i]/kjmol)
        print >> f, "    Number of QM energy levels: %i" % self.num_levels
        # derived quantities
        print >> f, "    Center [A]: % 8.2f % 8.2f % 8.2f" % tuple(self.center/angstrom)
        print >> f, "    Axis [1]: % 8.2f % 8.2f % 8.2f" % tuple(self.axis)
        print >> f, "    Moment [amu*bohr**2]: %f" % (self.moment/amu)
        print >> f, "    Reduced moment [amu*bohr**2]: %f" % (self.reduced_moment/amu)
        print >> f, "    Cancel wavenumber [1/cm]: %.1f" % (self.cancel_freq/(lightspeed/centimeter))
        print >> f, "    The cancelation wavenumber is %s." % self.cancel_method
        self.dump_values(f, "Energy levels [kJ/mol]", self.energy_levels/kjmol, "% 8.2f", 8)
        if self.hb is not None:
            print >> f, "    Number of basis functions: %i" % (self.hb.size)
        print >> f, "    Zero-point contribution [kJ/mol]: %.7f" % (self.free_energy(0.0)/kjmol)

    def plot_levels(self, prefix, temp, do_levels=True, do_data=True):
        """Plots the potential with the energy levels

           Arguments:
            | ``prefix`` -- A filename prefix for the png files.
            | ``temp`` -- A temperature that is used to indicate the statistical
                          weight of each level in the plots

           Optional argument:
            | ``do_levels`` -- When True, the energy levels are plotted.
                               [default=True]
            | ``do_data`` -- When True, the data points are plotted.
                             [default=True]

           One image will be generated:
            | ``${prefix}.png`` -- The potential and the energy levels.
        """
        import matplotlib.pyplot as pt
        pt.clf()
        if do_data:
            # plot the original potential data
            if self.rot_scan.potential is not None:
                angles, energies = self.potential
                pt.plot(angles/deg, energies/kjmol, "rx", mew=2)
            pt.axvline(self.nma_angle/deg, color="silver")
        # plot the fitted potential
        if self.hb is not None:
            step = 0.001
            x = np.arange(0.0, 2*np.pi*(1+0.5*step), 2*np.pi*step)
            v = self.hb.eval_fn(x, self.v_coeffs)
            pt.plot(x/deg, v/kjmol, "k-", linewidth=2)
        if do_levels:
            # plot the energy levels
            eks = self.energy_levels/(temp*boltzmann)
            bfs = np.exp(-eks)
            bfs /= bfs.sum()
            for i in xrange(self.num_levels):
                e = (self.energy_levels[i])/kjmol
                pt.axhline(e, color="b", linewidth=0.5)
                pt.axhline(e, xmax=bfs[i], color="b", linewidth=2)
        if self.rot_scan.potential is not None:
            angles, energies = self.potential
            pt.ylim(
                energies.min()/kjmol,
                1.5*energies.max()/kjmol
            )
            fit_energies = self.hb.eval_fn(angles, self.v_coeffs)
            rmsd = ((fit_energies - energies)**2).mean()**0.5
            rms = (energies**2).mean()**0.5
            rrmsd = rmsd/rms
            title = 'RMSD [kJ/mol] = %.1f    RRMSD [%%] = %.0f    dofmax = %i   rotsym = %i' % (
                rmsd, rrmsd*100, self.dofmax, self.rotsym
            )
            if self.even:
                title += '   even'
            pt.title(title)
        pt.xlim(0, 360)
        pt.ylabel("Energy [kJ/mol]")
        pt.xlabel("Dihedral angle [deg]")
        pt.savefig(prefix)

    def helper_terms(self, temp, n):
        """See :meth:`tamkin.partf.StatFysTerms.helper_terms`"""
        return np.array([
            -helper_vibrations(temp, n, self.cancel_freq, self.classical,
                                 self.freq_scaling, self.zp_scaling),
            helper_levels(temp, n, self.energy_levels, check=True) - temp**n*np.log(self.rotsym),
        ])

    def helpert_terms(self, temp, n):
        """See :meth:`tamkin.partf.StatFysTerms.helpert_terms`"""
        return np.array([
            -helpert_vibrations(temp, n, self.cancel_freq, self.classical,
                                  self.freq_scaling, self.zp_scaling),
            helpert_levels(temp, n, self.energy_levels, check=True),
        ])

    def helpertt_terms(self, temp, n):
        """See :meth:`tamkin.partf.StatFysTerms.helpertt_terms`"""
        return np.array([
            -helpertt_vibrations(temp, n, self.cancel_freq, self.classical,
                                   self.freq_scaling, self.zp_scaling),
            helpertt_levels(temp, n, self.energy_levels, check=True),
        ])
