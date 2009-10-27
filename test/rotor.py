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


from tamkin import *

import unittest
import numpy


__all__ = ["RotorTestCase"]


class RotorTestCase(unittest.TestCase):
    def assertArraysAlmostEqual(self, a, b, eps=1e-5, relative=False):
        self.assert_(isinstance(b, numpy.ndarray))
        self.assertEqual(a.shape, b.shape)
        if relative:
            self.assert_(abs(2*(a-b)/(a+b)).max() <= eps)
        else:
            self.assert_(abs(a-b).max() <= eps)

    def test_potential_op(self):
        a = 10.0
        hb = HarmonicBasis(3, a)
        coeffs = numpy.arange(6, dtype=float)*2
        op = hb.get_empty_op()
        hb._add_potential_op(op, coeffs)
        self.assertAlmostEqual(op[0,0], coeffs[2]/numpy.sqrt(2*a))
        self.assertAlmostEqual(op[0,1], coeffs[3]/numpy.sqrt(2*a))
        self.assertAlmostEqual(op[1,0], coeffs[3]/numpy.sqrt(2*a))
        self.assertAlmostEqual(op[1,1], -coeffs[2]/numpy.sqrt(2*a))

    def test_eval_fn(self):
        a = 10.0
        hb = HarmonicBasis(3, a)
        grid = numpy.arange(0.0, 10.01, 0.1)
        fn = hb.eval_fn(grid, [0,0,0,0,3,0])
        expected = 3*numpy.cos(grid*6.0*numpy.pi/a)/numpy.sqrt(a/2)
        self.assertArraysAlmostEqual(fn, expected)
        fn = hb.eval_fn(grid, [0,0,0,2,0,0])
        expected = 2*numpy.sin(grid*4.0*numpy.pi/a)/numpy.sqrt(a/2)
        self.assertArraysAlmostEqual(fn, expected)

    def test_fit_fn(self):
        a = 10.0
        hb = HarmonicBasis(10, a)
        grid = numpy.arange(0.0, 10.01, 1.0)
        f = numpy.exp(-((grid-5)/2)**2)
        coeffs = hb.fit_fn(grid, f)
        g = hb.eval_fn(grid, coeffs)
        self.assertArraysAlmostEqual(f, g)

    def test_fit_fn_sym(self):
        a = 9.0
        hb = HarmonicBasis(90, a)
        grid = numpy.arange(0.0, 1.501, 0.1)
        f = numpy.exp(-(grid/2)**2)
        f -= f.mean()
        coeffs = hb.fit_fn(grid, f, rotsym=3, even=True)
        g = hb.eval_fn(grid, coeffs)
        self.assertArraysAlmostEqual(f, g)
        grid = numpy.arange(0.0, 9.001, 0.1)
        g = hb.eval_fn(grid, coeffs)
        if False:
            import pylab
            pylab.clf()
            pylab.plot(grid, g)
            pylab.plot(grid[:16], f)
            pylab.savefig("output/test.png")
        self.assertArraysAlmostEqual(f, g[0:16])
        self.assertArraysAlmostEqual(f, g[30:46])
        self.assertArraysAlmostEqual(f, g[60:76])
        self.assertArraysAlmostEqual(f[::-1], g[15:31])
        self.assertArraysAlmostEqual(f[::-1], g[45:61])
        self.assertArraysAlmostEqual(f[::-1], g[75:91])

    def test_flat(self):
        a = 10.0
        mass = 1.0
        hb = HarmonicBasis(10, a)
        energies, orbitals = hb.solve(mass, numpy.zeros(hb.size))
        indexes = numpy.array([1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10])
        expected = 0.5/mass*(2*indexes*numpy.pi/a)**2
        self.assertArraysAlmostEqual(energies, expected, 1e-4)

    def test_harmonic_oscillator(self):
        a = 20.0
        hb = HarmonicBasis(20, a)
        x = numpy.arange(0.0, a, 0.1)
        v = 0.5*(x-a/2)**2
        #v = 5*(1+numpy.cos(2*numpy.pi*x/a))**2
        v_coeffs = hb.fit_fn(x, v)
        v_ref = -v.mean()
        energies, orbitals = hb.solve(1, v_coeffs)
        expected = numpy.arange(10) + 0.5

        if False:
            import pylab
            x = numpy.arange(0.0, a, 0.001)
            pylab.clf()
            for i in xrange(10):
                f = hb.eval_fn(x, orbitals[:,i])
                pylab.plot(x, f)
            pylab.savefig("output/wavefunctions.png")
            pylab.clf()
            v = hb.eval_fn(x, v_coeffs)
            pylab.plot(x, v-v_ref)
            for energy in energies[:10]:
                pylab.axhline(energy-v_ref)
            pylab.xlim(0,a)
            pylab.savefig("output/energy_levels.png")

        self.assertAlmostEqual(energies[0]-v_ref, 1.5, 1)
        self.assertAlmostEqual(energies[2]-v_ref, 3.5, 1)
        self.assertAlmostEqual(energies[4]-v_ref, 5.5, 1)
        self.assertAlmostEqual(energies[6]-v_ref, 7.5, 1)

