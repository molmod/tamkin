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

from molmod.units import cm, amu, kjmol, joule, mol, kelvin
from molmod.constants import lightspeed, boltzmann

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
        fn = hb.eval_fn(grid, [0,0,0,0,0,3,0])
        expected = 3*numpy.cos(grid*6.0*numpy.pi/a)/numpy.sqrt(a/2)
        self.assertArraysAlmostEqual(fn, expected)
        fn = hb.eval_fn(grid, [0,0,0,0,2,0,0])
        expected = 2*numpy.sin(grid*4.0*numpy.pi/a)/numpy.sqrt(a/2)
        self.assertArraysAlmostEqual(fn, expected)

    def test_fit_fn(self):
        a = 10.0
        hb = HarmonicBasis(10, a)
        grid = numpy.arange(0.0, 10.01, 1.0)
        f = numpy.exp(-((grid-5)/2)**2)
        coeffs = hb.fit_fn(grid, f, 10)
        g = hb.eval_fn(grid, coeffs)
        self.assertArraysAlmostEqual(f, g)

    def test_fit_fn_sym(self):
        a = 9.0
        hb = HarmonicBasis(90, a)
        grid = numpy.arange(0.0, 1.501, 0.1)
        f = numpy.exp(-(grid/2)**2)
        f -= f.mean()
        coeffs = hb.fit_fn(grid, f, 30, rotsym=3, even=True)
        g = hb.eval_fn(grid, coeffs)
        self.assertArraysAlmostEqual(f, g)
        grid = numpy.arange(0.0, 9.001, 0.1)
        g = hb.eval_fn(grid, coeffs)
        self.assertArraysAlmostEqual(f, g[0:16])
        self.assertArraysAlmostEqual(f, g[30:46])
        self.assertArraysAlmostEqual(f, g[60:76])
        self.assertArraysAlmostEqual(f[::-1], g[15:31])
        self.assertArraysAlmostEqual(f[::-1], g[45:61])
        self.assertArraysAlmostEqual(f[::-1], g[75:91])

        import pylab
        pylab.clf()
        pylab.plot(grid, g, "k-", lw=2)
        pylab.plot(grid[:16], f, "rx", mew=2)
        pylab.savefig("output/test_fit_fn_sym.png")

    def test_potential_op(self):
        a = 10.0
        mass = 1.0
        nmax = 10
        v_exp = numpy.zeros(nmax+1, complex)
        v_exp[0] = numpy.random.normal(0,1)
        v_exp[1:] += numpy.random.normal(0,1,nmax)
        v_exp[1:] += 1j*numpy.random.normal(0,1,nmax)
        #v_exp[3] = 1.0
        def get_v(index):
            if index>nmax or -index>nmax:
                return 0
            elif index>=0:
                return v_exp[index]
            else:
                return numpy.conjugate(v_exp[-index])
        v_op_exp = numpy.zeros((2*nmax+1,2*nmax+1), complex)
        for i0 in xrange(2*nmax+1):
            k0 = ((i0-1)/2+1)*(2*(i0%2)-1)
            for i1 in xrange(2*nmax+1):
                k1 = ((i1-1)/2+1)*(2*(i1%2)-1)
                #print (i0,i1), (k0,k1), k0-k1
                v_op_exp[i0,i1] = get_v(k0-k1)/numpy.sqrt(a)
        #for row in v_op_exp:
        #    print "".join({True: " ", False: "X"}[v==0] for v in row)
        hb = HarmonicBasis(nmax, a)
        v_cs = numpy.zeros(2*nmax+1, float)
        v_cs[0] = v_exp.real[0]
        v_cs[1::2] = numpy.sqrt(2.0)*v_exp.real[1:]
        v_cs[2::2] = -numpy.sqrt(2.0)*v_exp.imag[1:]
        v_op_cs = hb.get_empty_op()
        hb._add_potential_op(v_op_cs, v_cs)

        lc = numpy.array([
            [1.0, -1.0j],
            [1.0, 1.0j],
        ])/numpy.sqrt(2)

        lc_dagger = lc.transpose().conjugate()
        for i0 in xrange(nmax):
            for i1 in xrange(nmax):
                check = numpy.dot(lc_dagger, numpy.dot(v_op_exp[2*i0+1:2*i0+3,2*i1+1:2*i1+3], lc))
                self.assert_(abs(check.imag).max() < 1e-3)
                check = check.real
                self.assertArraysAlmostEqual(
                    v_op_cs[2*i0+1:2*i0+3,2*i1+1:2*i1+3],
                    check,
                    1e-3
                )

    def test_flat(self):
        a = 10.0
        mass = 1.0
        hb = HarmonicBasis(10, a)
        energies, orbitals = hb.solve(mass, numpy.zeros(hb.size), evecs=True)

        import pylab
        x = numpy.arange(0.0, a, 0.001)
        pylab.clf()
        for i in xrange(10):
            f = hb.eval_fn(x, orbitals[:,i])
            pylab.plot(x, f+i)
        pylab.savefig("output/flat_wavefunctions.png")

        indexes = numpy.array([0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10])
        expected = 0.5/mass*(2*indexes*numpy.pi/a)**2
        self.assertArraysAlmostEqual(energies, expected, 1e-4)

    def test_flat2(self):
        molecule = load_molecule_g03fchk("input/ethane/gaussian.fchk")
        nma = NMA(molecule)
        dihedral, angles, energies, geometries, top_indexes = load_rotscan_g03(
            "input/rotor/gaussian.com", "input/rotor/gaussian.fchk"
        )
        energies[:] = 0.0
        cancel_freq = compute_cancel_frequency(molecule, top_indexes)
        self.assertAlmostEqual(cancel_freq/lightspeed*cm, 314, 0)
        rotor1 = Rotor(
            top_indexes, cancel_freq, rotsym=3, even=True,
            potential=(angles, energies, 5), num_levels=50
        )
        pf1 = PartFun(nma, [
            ExternalTranslation(),
            ExternalRotation(6),
            rotor1,
        ])
        rotor2 = Rotor(
            top_indexes, cancel_freq, rotsym=3, even=True,
            potential=None, num_levels=50
        )
        pf2 = PartFun(nma, [
            ExternalTranslation(),
            ExternalRotation(6),
            rotor2,
        ])
        self.assertArraysAlmostEqual(rotor1.energy_levels, rotor2.energy_levels)


    def test_harmonic(self):
        a = 20.0
        hb = HarmonicBasis(20, a)
        x = numpy.arange(0.0, a, 0.1)
        v = 0.5*(x-a/2)**2
        #v = 5*(1+numpy.cos(2*numpy.pi*x/a))**2
        v_coeffs = hb.fit_fn(x, v, 20, even=True)
        energies, orbitals = hb.solve(1, v_coeffs, evecs=True)
        expected = numpy.arange(10) + 0.5
        self.assertAlmostEqual(energies[0], 0.5, 1)
        self.assertAlmostEqual(energies[1], 1.5, 1)
        self.assertAlmostEqual(energies[2], 2.5, 1)
        self.assertAlmostEqual(energies[3], 3.5, 1)
        self.assertAlmostEqual(energies[4], 4.5, 1)
        self.assertAlmostEqual(energies[5], 5.5, 1)
        self.assertAlmostEqual(energies[6], 6.5, 1)

        import pylab
        x = numpy.arange(0.0, a, 0.001)
        pylab.clf()
        for i in xrange(10):
            f = hb.eval_fn(x, orbitals[:,i])
            pylab.plot(x, f)
        pylab.savefig("output/harmonic_wavefunctions.png")
        pylab.clf()
        v = hb.eval_fn(x, v_coeffs)
        pylab.plot(x, v)
        for energy in energies[:10]:
            pylab.axhline(energy)
        pylab.xlim(0,a)
        pylab.savefig("output/harmonic_levels.png")

    def test_ethane_hindered(self):
        molecule = load_molecule_g03fchk("input/ethane/gaussian.fchk")
        nma = NMA(molecule)
        dihedral, angles, energies, geometries, top_indexes = load_rotscan_g03(
            "input/rotor/gaussian.com", "input/rotor/gaussian.fchk"
        )
        cancel_freq = compute_cancel_frequency(molecule, top_indexes)
        self.assertAlmostEqual(cancel_freq/lightspeed*cm, 314, 0)
        rotor = Rotor(
            top_indexes, cancel_freq, rotsym=3, even=True,
            potential=(angles, energies, 5), num_levels=50
        )
        pf = PartFun(nma, [
            ExternalTranslation(),
            ExternalRotation(6),
            rotor,
        ])
        self.assertArraysAlmostEqual(
            rotor.hb.eval_fn(angles, rotor.v_coeffs),
            energies
        )
        # reference data from legacy code (Veronique & co)
        self.assertAlmostEqual(rotor.absolute_moment/amu, 11.092362911176032, 2)
        self.assertAlmostEqual(rotor.relative_moment/amu, 5.5461814555880098, 2)
        self.assertAlmostEqual(numpy.exp(rotor.log_eval_terms(100.0)[1]), 0.12208E+00, 1)
        self.assertAlmostEqual(rotor.heat_capacity_terms(100.0)[1]/(joule/mol/kelvin), 2.567, 0)
        self.assertAlmostEqual(rotor.entropy_terms(100.0)[1]/(joule/mol), 0.766, 0)
        self.assertAlmostEqual(numpy.exp(rotor.log_eval_terms(800.0)[1]), 0.21108E+01, 1)
        self.assertAlmostEqual(rotor.heat_capacity_terms(800.0)[1]/(joule/mol/kelvin), 6.346, 1)
        self.assertAlmostEqual(rotor.entropy_terms(800.0)[1]/(joule/mol), 14.824, 1)

        rotor.plot_levels("output/ethane_hindered_levels.png", 300)
        pf.write_to_file("output/ethane_hindered.txt")
        ta = ThermoAnalysis(pf, [200,300,400,500,600,700,800,900])
        ta.write_to_file("output/ethane_hindered_thermo.csv")

    def test_legacy1(self):
        a = 2*numpy.pi
        mass = 5.5*amu
        hb = HarmonicBasis(100, a)
        v_coeffs = numpy.zeros(hb.size, float)
        v_coeffs[0] = 0.5*11.5*kjmol*numpy.sqrt(a)
        v_coeffs[5] = 0.5*11.5*kjmol*numpy.sqrt(a/2)
        self.assertArraysAlmostEqual(
            hb.eval_fn(numpy.array([0, a/6]), v_coeffs),
            numpy.array([11.5*kjmol, 0.0]),
        )
        energies, orbitals = hb.solve(mass, v_coeffs, evecs=True)

        import pylab
        x = numpy.arange(0.0, a, 0.001)
        pylab.clf()
        for i in xrange(10):
            f = hb.eval_fn(x, orbitals[:,i])
            pylab.plot(x, f+i)
        pylab.savefig("output/legacy_wavefunctions.png")

        # check energy levels
        expected = numpy.array([
            1.7635118, 1.76361979, 5.11795465, 8.04553104, 8.1095722,
            10.3876796, 11.8999683, 12.9078395, 14.6739639, 16.7836847,
            19.1722507,
            1.76361979,  5.11795465, 5.12218335, 8.1095722, 10.3876796,
            10.8661504, 12.9078395, 14.6739639, 16.7544718, 19.1722507,
        ])
        expected.sort()
        self.assertArraysAlmostEqual(energies[:10]/kjmol, expected[:10], 1e-3)
        self.assertAlmostEqual(numpy.exp(-energies/(100*boltzmann)).sum()/3.0, 0.12208E+00, 5)

