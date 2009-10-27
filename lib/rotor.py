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


import numpy


__all__ = ["HarmonicBasis"]


class HarmonicBasis(object):
    def __init__(self, nmax, a):
        """Initialize a harmonic basis

           Basis functions:
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

    size = property(lambda self: 2*self.nmax)

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
        diag[::2] = tmp
        diag[1::2] = tmp

    def _add_potential_op(self, op, potential):
        """Add the potential energy to the given operator"""
        # blocks
        c = potential[::2]/(numpy.sqrt(2*self.a))
        s = potential[1::2]/(numpy.sqrt(2*self.a))
        for i0 in xrange(self.nmax):
            for i1 in xrange(self.nmax):
                k = (i0+1)+(i1+1)-1
                if k < self.nmax:
                    op[2*i0  ,2*i1  ] += c[k] # cc
                    op[2*i0+1,2*i1+1] -= c[k] # ss
                    op[2*i0  ,2*i1+1] += s[k] # cs
                    op[2*i0+1,2*i1  ] += s[k] # sc
                k = abs(i0-i1)-1
                if k>= 0 and k < self.nmax:
                    op[2*i0  ,2*i1  ] += c[k] # cc
                    op[2*i0+1,2*i1+1] += c[k] # ss
                    op[2*i0  ,2*i1+1] -= s[k] # cs
                    op[2*i0+1,2*i1  ] -= s[k] # sc

    def solve(self, mass, potential):
        """Return the energies and wavefunctions for the given mass and potential

           Arguments:
             mass  --  the mass of the particle
             potential  --  the expansion coefficients of the potential energy
        """
        H = self.get_hamiltonian_op(mass, potential)
        return numpy.linalg.eigh(H)

    def eval_fn(self, grid, coeffs):
        """Evaluate the function represented by coeffs on the given grid

           Arguments:
             grid  --  the values at which the function must be evaluated
             coeffs  --  the expansion coefficients
        """
        result = numpy.zeros(grid.shape, float)
        for i in xrange(0,self.nmax):
            arg = 2*(i+1)*numpy.pi*grid/self.a
            result += (coeffs[2*i  ]/numpy.sqrt(self.a/2))*numpy.cos(arg)
            result += (coeffs[2*i+1]/numpy.sqrt(self.a/2))*numpy.sin(arg)
        return result

    def fit_fn(self, grid, f, rotsym=1, even=False, rcond=1e-10):
        """Fit the expansion coefficients that represent function f

           Arguments:
             grid  --  the x values on which the function f is known
             f  --  the function to be represented by expansion coefficients

           Optional arguments:
             rotsym  --  impose this rotational symmetry (default=1)
             even  --  only fit even functions, i.e. cosines (default=False)
             rcond  --  the cutoff for the singular values in the least squares
                        fit (default=1e-10)
        """
        if rotsym < 1:
            raise ValueError("rotym must be at least 1")
        if even:
            A = numpy.zeros((len(grid), self.nmax/rotsym), float)
        else:
            A = numpy.zeros((len(grid), self.size/rotsym), float)
        counter = 0
        for i in xrange(self.nmax/rotsym):
            arg = 2*(i+1)*rotsym*numpy.pi*grid/self.a
            A[:,counter] = numpy.cos(arg)/numpy.sqrt(self.a/2)
            counter += 1
            if not even:
                A[:,counter] = numpy.sin(arg)/numpy.sqrt(self.a/2)
                counter += 1
        coeffs, residuals, rank, S = numpy.linalg.lstsq(A, f, rcond)
        result = numpy.zeros(self.size)
        if even:
            result[2*rotsym-2::2*rotsym] = coeffs
        else:
            result[2*rotsym-2::2*rotsym] = coeffs[::2]
            result[2*rotsym-1::2*rotsym] = coeffs[1::2]
        return result

