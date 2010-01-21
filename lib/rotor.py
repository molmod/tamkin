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
"""Implementation of free and hindered rotors

   The current implementation supports one-dimensional free and hindered rotors.
   For practical applications it is apparently not necessary to consider higher-
   dimensional hindered rotors. [1]

   [1] Chemical Physics, Vol. 328 (1-3) 251 - 258, 2006
"""

from tamkin.partf import Info, StatFysTerms, helper0_vibrations, \
    helper1_vibrations, helper2_vibrations, helper0_levels, helper1_levels, \
    helper2_levels
from tamkin.nma import NMA, MBH
from tamkin.geom import transrot_basis

from molmod.units import deg, kjmol, angstrom, cm, amu
from molmod.constants import boltzmann, lightspeed

import numpy


__all__ = [
    "HarmonicBasis", "compute_cancel_frequency", "compute_moments", "Rotor"
]


class HarmonicBasis(object):
    """A harmonic basis set for periodic one-dimensional QM systems

       In addition to the definition of the basis set, this class also
       implements the kinetic and potential energy operators required for the
       solution of the schrodinger equation. The workflow is as follows:

       >>> a = 10.0                             # the size of the system
       >>> hb = HarmonicBasis(10, a)            # create basis object
       >>> grid = numpy.arange(0.0, 10.01, 1.0) # define a grid
       >>> v = -numpy.exp(-((grid-5)/2)**2)     # define a potential on the grid
       >>> v_coeffs = hb.fit_fn(grid, v, 10)    # expand the potential in the basis
       >>> mass = 1.0
       >>> energies, wfns = hb.solve(mass, v_coeffs, evecs=True) # solve problem
    """
    def __init__(self, nmax, a):
        """Initialize a harmonic basis

           Basis functions:
             1/sqrt(a/2)
             cos(2*pi*x/a)/sqrt(a/2)
             sin(2*pi*x/a)/sqrt(a/2)
             ...
             cos(nmax*2*pi*x/a)/sqrt(a/2)
             sin(nmax*2*pi*x/a)/sqrt(a/2)

           This is an orthonormal basis. This object has methods to construct a
           Hamiltonian operator in this basis and to transform a function on a
           grid into an expansion in this basis (and back).
        """
        self.nmax = nmax
        self.a = a

    size = property(lambda self: 2*self.nmax+1)

    def get_empty_op(self):
        """Returns an empty operator (zero)"""
        return numpy.zeros((self.size, self.size), float)

    def get_hamiltonian_op(self, mass, potential):
        """Returns the Hamiltonian operator for the given mass and potential

           Arguments:
             mass  --  the mass of the particle
             potential  --  the expansion coefficients of the potential energy
        """
        op = self.get_empty_op()
        self._add_kinetic_op(op, mass)
        self._add_potential_op(op, potential)
        return op

    def _add_kinetic_op(self, op, mass):
        """Add the kinetic energy to the given operator"""
        factor = (0.5/mass*(2*numpy.pi/self.a)**2)
        tmp = factor*numpy.arange(1, self.nmax+1)**2
        diag = op.ravel()[::self.size+1]
        diag[1::2] += tmp
        diag[2::2] += tmp

    def _add_potential_op(self, op, potential):
        """Add the potential energy to the given operator"""
        c = numpy.zeros(self.nmax+1,float)
        s = numpy.zeros(self.nmax+1,float)
        c[0] = potential[0]/numpy.sqrt(self.a)
        c[1:] = potential[1::2]/numpy.sqrt(self.a)
        s[1:] = potential[2::2]/numpy.sqrt(self.a)
        op[0,0] += c[0]
        for i0 in xrange(1,self.nmax+1):
            op[2*i0-1,0] += c[i0]
            op[2*i0-0,0] += s[i0]
            op[0,2*i0-1] += c[i0]
            op[0,2*i0-0] += s[i0]
        c[1:] /= numpy.sqrt(2)
        s[1:] /= numpy.sqrt(2)
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
             mass  --  the mass of the particle
             potential  --  the expansion coefficients of the potential energy
        """
        H = self.get_hamiltonian_op(mass, potential)
        if evecs:
            return numpy.linalg.eigh(H)
        else:
            return numpy.linalg.eigvalsh(H)

    def eval_fn(self, grid, coeffs):
        """Evaluate the function represented by coeffs on the given grid

           Arguments:
             grid  --  the values at which the function must be evaluated
             coeffs  --  the expansion coefficients
        """
        result = numpy.zeros(grid.shape, float) + coeffs[0]/numpy.sqrt(self.a)
        for i in xrange(self.nmax):
            arg = ((i+1)*2*numpy.pi/self.a)*grid
            result += (coeffs[2*i+1]/numpy.sqrt(self.a/2))*numpy.cos(arg)
            result += (coeffs[2*i+2]/numpy.sqrt(self.a/2))*numpy.sin(arg)
        return result

    def fit_fn(self, grid, f, dofmax, rotsym=1, even=False, rcond=0.0):
        """Fit the expansion coefficients that represent function f

           Arguments:
             grid  --  the x values on which the function f is known
             f  --  the function to be represented by expansion coefficients
             dofmax  --  the maximum number of cosines in the fit
                         when even==False, the same number of sines is also
                         included

           Optional arguments:
             rotsym  --  impose this rotational symmetry (default=1)
             even  --  only fit even functions, i.e. cosines (default=False)
             rcond  --  the cutoff for the singular values in the least squares
                        fit (default=1e-10)
        """
        if rotsym < 1:
            raise ValueError("rotym must be at least 1")
        ncos = min(dofmax, self.nmax/rotsym)
        if even:
            A = numpy.zeros((len(grid), ncos+1), float)
        else:
            A = numpy.zeros((len(grid), 2*ncos+1), float)
        A[:,0] = 1.0/numpy.sqrt(self.a)
        counter = 1
        for i in xrange(ncos):
            arg = ((i+1)*rotsym*2*numpy.pi/self.a)*grid
            A[:,counter] = numpy.cos(arg)/numpy.sqrt(self.a/2)
            counter += 1
            if not even:
                A[:,counter] = numpy.sin(arg)/numpy.sqrt(self.a/2)
                counter += 1

        coeffs, residuals, rank, S = numpy.linalg.lstsq(A, f, rcond)

        result = numpy.zeros(self.size)
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


def compute_cancel_frequency(molecule, dihedral, top_indexes):
    """Compute the frequency of the rotor in the HO approximation

       This function is based on the MNH method and returns the frequency that
       has to be canceled when this mode is replaced by a free or hindered
       rotor.

       Arguments:
         molecule  --  A Molecule object (see data.py)
         top_indexes  --  the indexes of the rotor atoms
    """
    axis = list(dihedral[1:3])
    other_top_indexes = list(set(xrange(molecule.size)) - set(top_indexes) - set(axis))
    blocks = [
        axis + top_indexes,
        axis + other_top_indexes,
    ]
    nma = NMA(molecule, MBH(blocks))
    if len(nma.freqs) != 7:
        raise RuntimeError("Expecting 7 frequencies, got %i" % len(nma.freqs))
    non_zero = [i for i in xrange(7) if i not in nma.zeros][0]
    return nma.freqs[non_zero]


def compute_moments(coordinates, masses3, center, axis, indexes):
    # the derivative of the cartesian coordinates towards the rotation
    # angle of the top:
    rot_tangent = numpy.zeros((len(coordinates)*3), float)
    for i in indexes:
        rot_tangent[3*i:3*i+3] = numpy.cross(coordinates[i]-center, axis)
    # transform this to a derivative without global linear or angular
    # momentum. This is done by projecting on the basis of external degrees
    # of freedom in mass-weighted coordinates, and subsequently subtracting
    # that projection from the original tangent
    basis = transrot_basis(coordinates)
    A = numpy.dot(basis*masses3, basis.transpose())
    B = numpy.dot(basis*masses3, rot_tangent)
    alphas = numpy.linalg.solve(A,B)
    rot_tangent_relative = rot_tangent - numpy.dot(alphas, basis)
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
    def __init__(self, dihedral, indexes, cancel_freq, suffix=None, rotsym=1, even=False, potential=None, num_levels=50):
        """Initialize a Rotor term for the partition function

           Arguments:
             dihedral  --  the index of the atoms that define the dihedral angle
             indexes  --  a list of atom indexes involved in the rotor.
             cancel_freq  --  the frequency to cancel in the vibrational
                              partition function

           Optional arguments:
             suffix  --  a name suffix used to distinguish between different
                         rotors
             rotsym  --  the rotational symmetry of the rotor (default=1)
             even    --  True of the rotor is not chiral, i.e. when it has a
                         even potential
             potential  --  a tuple with two arrays, the first containing the
                            angles and the second containing the corresponding
                            energies
             num_levels  --  the number of energy levels considered in the
                             QM treatment of the rotor (default=50)
        """
        if len(dihedral) != 4:
            raise ValueError("The first argument must be a list of 4 integers")
        if len(indexes) < 1:
            raise ValueError("A rotor must have at least one atoms")
        self.dihedral = dihedral
        self.indexes = indexes
        self.cancel_freq = cancel_freq
        if suffix is None:
            self.suffix = "_".join(str(i) for i in indexes)
        else:
            self.suffix = suffix
        self.rotsym = rotsym
        self.even = even
        self.potential = potential
        self.num_levels = num_levels
        if self.potential is None:
            Info.__init__(self, "free_rotor_%s" % self.suffix)
        else:
            Info.__init__(self, "hindered_rotor_%s" % self.suffix)
        StatFysTerms.__init__(self, 2) # two terms

    def init_part_fun(self, nma, partf):
        """Compute all the partition function parameters based on the nma

           This method is part of the PartFun API. It should never be called
           directly.

           Arguments:
             nma  --  An NMA object (see nma.py)
        """
        if nma.periodic:
            raise NotImplementedError("Rotors in periodic systems are not supported yet")
        self.center = nma.coordinates[self.dihedral[1]]
        self.axis = nma.coordinates[self.dihedral[2]] - self.center
        self.axis /= numpy.linalg.norm(self.axis)
        self.moment, self.reduced_moment = compute_moments(
            nma.coordinates, nma.masses3, self.center, self.axis, self.indexes
        )
        # the energy levels
        if self.potential is None:
            # free rotor
            self.energy_levels = numpy.zeros(self.num_levels, float)
            for i in xrange(self.num_levels-1):
                index = i/2+1
                self.energy_levels[i+1] = index**2/(2*self.reduced_moment)
            self.hb = None
            self.v_coeffs = None
            self.v_ref = 0.0
        else:
            # hindered rotor
            a = 2*numpy.pi
            self.hb = HarmonicBasis(self.num_levels, a)
            angles, energies, dofmax = self.potential
            # apply periodic boundary conditions
            angles -= numpy.floor(angles/a)*a
            # set reference energy, which is take to be the energy of the geometry
            # with a dihedral angle closest to the dihedral angle of the refrence
            # geometry in the nma object. We can not take the nma energy because
            # the rotational energy barriers may be computed at another level
            # that the nma energy.
            from molmod.ic import dihed_angle
            nma_angle = dihed_angle(
                nma.coordinates[self.dihedral[0]], nma.coordinates[self.dihedral[1]],
                nma.coordinates[self.dihedral[2]], nma.coordinates[self.dihedral[3]],
            )[0]
            deltas = angles - nma_angle
            # apply periodic boundary conditions
            deltas -= numpy.floor(deltas/a)*a
            # get the right energy
            energies -= energies[deltas.argmin()]

            self.v_coeffs = self.hb.fit_fn(angles, energies, dofmax, self.rotsym, self.even)
            self.energy_levels = self.hb.solve(self.reduced_moment, self.v_coeffs)
            self.energy_levels = self.energy_levels[:self.num_levels]
        # scaling factors
        self.freq_scaling = partf.vibrational.freq_scaling
        self.zp_scaling = partf.vibrational.zp_scaling
        self.classical = partf.vibrational.classical

    def dump(self, f):
        """Write all the information about the rotor to a file

           This method is part of the PartFun API and should never be called
           directly. It will only work properly once the init_part_fun method
           is called.

           Arguments:
             f  --  a file-like object
        """
        Info.dump(self, f)
        # parameters
        print >> f, "    Indexes: %s" % " ".join(str(i) for i in self.indexes)
        print >> f, "    Rotational symmetry: %i" % self.rotsym
        print >> f, "    Even potential: %s" % self.even
        if self.potential is None:
            print >> f, "    This is a free rotor"
        else:
            angles, energies, dofmax = self.potential
            print >> f, "    This is a hindered rotor"
            print >> f, "    Maximum number of cosines in the fit: %i" % dofmax
            print >> f, "    Potential: Angle [deg]    Energy [kJ/mol]"
            for i in xrange(len(angles)):
                print >> f, "              % 7.2f         %6.1f" % (angles[i]/deg, energies[i]/kjmol)
        print >> f, "    Number of QM energy levels: %i" % self.num_levels
        # derived quantities
        print >> f, "    Center [A]: % 8.2f % 8.2f % 8.2f" % tuple(self.center/angstrom)
        print >> f, "    Axis [1]: % 8.2f % 8.2f % 8.2f" % tuple(self.axis)
        print >> f, "    Moment [amu*bohr**2]: %f" % (self.moment/amu)
        print >> f, "    Reduced moment [amu*bohr**2]: %f" % (self.reduced_moment/amu)
        print >> f, "    Cancel wavenumber [1/cm]: %.1f" % (self.cancel_freq/(lightspeed/cm))
        self.dump_values(f, "Energy levels [kJ/mol]", self.energy_levels/kjmol, "% 8.2f", 8)
        if self.hb is not None:
            print >> f, "    Number of basis functions: %i" % (self.hb.size)
        print >> f, "    Free energy contribution at T=0 [kJ/mol]: %.7f" % (self.free_energy(0.0)/kjmol)

    def plot_levels(self, filename, temp, num=20):
        """Plots the potential with the energy levels

           Arguments:
             prefix  --  a filename prefix for the png files
             temp  --  a temperature that is used to indicate the statistical
                       weight of each level in the plots

           Optional argument:
             num  --  the number of energy levels and wavefunctions to be
                      plotted (default=10)

           One image will be generated:
             ${prefix}.png  --  the potential and the energy levels
        """
        import pylab
        pylab.clf()
        # plot the original potential data
        if self.potential is not None:
            angles, energies, dofmax = self.potential
            pylab.plot(angles/deg, energies/kjmol, "rx", mew=2)
        # plot the fitted potential
        if self.hb is not None:
            step = 0.001
            x = numpy.arange(0.0, 2*numpy.pi*(1+0.5*step), 2*numpy.pi*step)
            v = self.hb.eval_fn(x, self.v_coeffs)
            pylab.plot(x/deg, v/kjmol, "k-", linewidth=2)
        # plot the energy levels
        eks = self.energy_levels/(temp*boltzmann)
        bfs = numpy.exp(-eks)
        bfs /= bfs.sum()
        for i in xrange(min(num, self.num_levels)):
            e = (self.energy_levels[i])/kjmol
            pylab.axhline(e, color="b", linewidth=0.5)
            pylab.axhline(e, xmax=bfs[i], color="b", linewidth=2)
        pylab.xlim(0,360)
        if self.potential is not None:
            pylab.ylim(self.potential[1].min()/kjmol, 1.5*self.potential[1].max()/kjmol)
        pylab.ylabel("Energy [kjmol]")
        pylab.xlabel("Dihedral angle [deg]")
        pylab.savefig(filename)

    def helper0_terms(self, temp, n):
        """Compute the conbtributions to the logarithm of the partition function

           Argument:
             temp  --  the temperature

           Returns: 1D numpy arrays with contributions
        """
        return numpy.array([
            -helper0_vibrations(temp, n, self.cancel_freq, self.classical,
                                 self.freq_scaling, self.zp_scaling),
            helper0_levels(temp, n, self.energy_levels) - temp**n*numpy.log(self.rotsym),
        ])

    def helper1_terms(self, temp, n):
        """Compute the conbtributions to the derivative of the logarithm of the
           partition function

           Argument:
             temp  --  the temperature

           Returns: 1D numpy arrays with contributions
        """
        return numpy.array([
            -helper1_vibrations(temp, n, self.cancel_freq, self.classical,
                                  self.freq_scaling, self.zp_scaling),
            helper1_levels(temp, n, self.energy_levels),
        ])

    def helper2_terms(self, temp, n):
        """Compute the conbtributions to the second derivative of the logarithm
           of the partition function

           Argument:
             temp  --  the temperature

           Returns: 1D numpy arrays with contributions
        """
        return numpy.array([
            -helper2_vibrations(temp, n, self.cancel_freq, self.classical,
                                   self.freq_scaling, self.zp_scaling),
            helper2_levels(temp, n, self.energy_levels),
        ])


