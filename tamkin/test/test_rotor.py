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


import os
import pkg_resources
import numpy as np
import unittest

from molmod import centimeter, amu, kjmol, joule, mol, kelvin, lightspeed, boltzmann
from molmod.periodic import periodic
from molmod.io import XYZFile
from molmod.test.common import tmpdir

from tamkin import *


__all__ = ["RotorTestCase"]


class RotorTestCase(unittest.TestCase):
    def assertArraysAlmostEqual(self, a, b, eps=1e-5, relative=False):
        self.assert_(isinstance(b, np.ndarray))
        self.assertEqual(a.shape, b.shape)
        if relative:
            self.assert_(abs(2*(a-b)/(a+b)).max() <= eps)
        else:
            self.assert_(abs(a-b).max() <= eps)

    def test_potential_op(self):
        a = 10.0
        hb = HarmonicBasis(3, a)
        coeffs = np.arange(6, dtype=float)*2
        op = hb.get_empty_op()
        hb._add_potential_op(op, coeffs)
        self.assertAlmostEqual(op[0,0], coeffs[2]/np.sqrt(2*a))
        self.assertAlmostEqual(op[0,1], coeffs[3]/np.sqrt(2*a))
        self.assertAlmostEqual(op[1,0], coeffs[3]/np.sqrt(2*a))
        self.assertAlmostEqual(op[1,1], -coeffs[2]/np.sqrt(2*a))

    def test_eval_fn(self):
        a = 10.0
        hb = HarmonicBasis(3, a)
        grid = np.arange(0.0, 10.01, 0.1)
        fn = hb.eval_fn(grid, [0,0,0,0,0,3,0])
        expected = 3*np.cos(grid*6.0*np.pi/a)/np.sqrt(a/2)
        self.assertArraysAlmostEqual(fn, expected)
        fn = hb.eval_fn(grid, [0,0,0,0,2,0,0])
        expected = 2*np.sin(grid*4.0*np.pi/a)/np.sqrt(a/2)
        self.assertArraysAlmostEqual(fn, expected)

    def test_eval_deriv(self):
        a = 10.0
        hb = HarmonicBasis(3, a)
        grid = np.arange(0.0, 10.01, 1.0)
        coeffs = [-0.5,1.2,2.3,-0.7,0.1,0.3,-1.0]
        eps = 1e-6
        aderiv = hb.eval_deriv(grid, coeffs)
        nderiv = (hb.eval_fn(grid+eps, coeffs) - hb.eval_fn(grid-eps, coeffs))/(2*eps)
        self.assertArraysAlmostEqual(aderiv, nderiv)

    def test_eval_deriv2(self):
        a = 10.0
        hb = HarmonicBasis(3, a)
        grid = np.arange(0.0, 10.01, 1.0)
        coeffs = [-0.5,1.2,2.3,-0.7,0.1,0.3,-1.0]
        eps = 1e-6
        aderiv2 = hb.eval_deriv2(grid, coeffs)
        nderiv2 = (hb.eval_deriv(grid+eps, coeffs) - hb.eval_deriv(grid-eps, coeffs))/(2*eps)
        self.assertArraysAlmostEqual(aderiv2, nderiv2)

    def test_fit_fn(self):
        a = 10.0
        hb = HarmonicBasis(10, a)
        grid = np.arange(0.0, 10.01, 1.0)
        f = np.exp(-((grid-5)/2)**2)
        coeffs = hb.fit_fn(grid, f, 10)
        g = hb.eval_fn(grid, coeffs)
        self.assertArraysAlmostEqual(f, g)

    def test_fit_fn_sym(self):
        a = 9.0
        hb = HarmonicBasis(90, a)
        grid = np.arange(0.0, 1.501, 0.1)
        f = np.exp(-(grid/2)**2)
        f -= f.mean()
        coeffs = hb.fit_fn(grid, f, 30, rotsym=3, even=True)
        g = hb.eval_fn(grid, coeffs)
        self.assertArraysAlmostEqual(f, g)
        grid = np.arange(0.0, 9.001, 0.1)
        g = hb.eval_fn(grid, coeffs)
        self.assertArraysAlmostEqual(f, g[0:16])
        self.assertArraysAlmostEqual(f, g[30:46])
        self.assertArraysAlmostEqual(f, g[60:76])
        self.assertArraysAlmostEqual(f[::-1], g[15:31])
        self.assertArraysAlmostEqual(f[::-1], g[45:61])
        self.assertArraysAlmostEqual(f[::-1], g[75:91])

        import matplotlib.pyplot as pt
        pt.clf()
        pt.plot(grid, g, "k-", lw=2)
        pt.plot(grid[:16], f, "rx", mew=2)
        with tmpdir(__name__, 'test_fit_fn_sym') as dn:
            pt.savefig(os.path.join(dn, "test_fit_fn_sym.png"))

    def test_potential_op(self):
        a = 10.0
        mass = 1.0
        nmax = 10
        v_exp = np.zeros(nmax+1, complex)
        v_exp[0] = np.random.normal(0,1)
        v_exp[1:] += np.random.normal(0,1,nmax)
        v_exp[1:] += 1j*np.random.normal(0,1,nmax)
        #v_exp[3] = 1.0
        def get_v(index):
            if index>nmax or -index>nmax:
                return 0
            elif index>=0:
                return v_exp[index]
            else:
                return np.conjugate(v_exp[-index])
        v_op_exp = np.zeros((2*nmax+1,2*nmax+1), complex)
        for i0 in xrange(2*nmax+1):
            k0 = ((i0-1)/2+1)*(2*(i0%2)-1)
            for i1 in xrange(2*nmax+1):
                k1 = ((i1-1)/2+1)*(2*(i1%2)-1)
                #print (i0,i1), (k0,k1), k0-k1
                v_op_exp[i0,i1] = get_v(k0-k1)/np.sqrt(a)
        #for row in v_op_exp:
        #    print "".join({True: " ", False: "X"}[v==0] for v in row)
        hb = HarmonicBasis(nmax, a)
        v_cs = np.zeros(2*nmax+1, float)
        v_cs[0] = v_exp.real[0]
        v_cs[1::2] = np.sqrt(2.0)*v_exp.real[1:]
        v_cs[2::2] = -np.sqrt(2.0)*v_exp.imag[1:]
        v_op_cs = hb.get_empty_op()
        hb._add_potential_op(v_op_cs, v_cs)

        lc = np.array([
            [1.0, -1.0j],
            [1.0, 1.0j],
        ])/np.sqrt(2)

        lc_dagger = lc.transpose().conjugate()
        for i0 in xrange(nmax):
            for i1 in xrange(nmax):
                check = np.dot(lc_dagger, np.dot(v_op_exp[2*i0+1:2*i0+3,2*i1+1:2*i1+3], lc))
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
        energies, orbitals = hb.solve(mass, np.zeros(hb.size), evecs=True)

        import matplotlib.pyplot as pt
        x = np.arange(0.0, a, 0.001)
        pt.clf()
        for i in xrange(10):
            f = hb.eval_fn(x, orbitals[:,i])
            pt.plot(x, f+i)
        with tmpdir(__name__, 'test_flat') as dn:
            pt.savefig(os.path.join(dn, "flat_wavefunctions.png"))

        indexes = np.array([0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10])
        expected = 0.5/mass*(2*indexes*np.pi/a)**2
        self.assertArraysAlmostEqual(energies, expected, 1e-4)

    def test_flat2(self):
        molecule = load_molecule_g03fchk(
            pkg_resources.resource_filename(__name__, "../data/test/ethane/gaussian.fchk"))
        nma = NMA(molecule)
        rotscan1 = load_rotscan_g03log(
            pkg_resources.resource_filename(__name__, "../data/test/rotor/gaussian.log"))
        my_potential = rotscan1.potential.copy()
        my_potential[1][:] = nma.energy
        rotscan1 = rotscan1.copy_with(potential=my_potential)
        rotor1 = Rotor(rotscan1, molecule, rotsym=3, even=True)
        self.assertAlmostEqual(rotor1.cancel_freq/lightspeed*centimeter, 314, 0)
        pf1 = PartFun(nma, [ExtTrans(), ExtRot(6), rotor1])
        rotscan2 = RotScan(rotscan1.dihedral, top_indexes=rotscan1.top_indexes)
        rotor2 = Rotor(rotscan2, molecule, rotsym=3, even=True)
        pf2 = PartFun(nma, [ExtTrans(), ExtRot(6), rotor2])
        self.assertArraysAlmostEqual(rotor1.energy_levels, rotor2.energy_levels)

    def test_harmonic(self):
        a = 20.0
        hb = HarmonicBasis(20, a)
        x = np.arange(0.0, a, 0.1)
        v = 0.5*(x-a/2)**2
        #v = 5*(1+np.cos(2*np.pi*x/a))**2
        v_coeffs = hb.fit_fn(x, v, 20, even=True, v_threshold=0.1)
        energies, orbitals = hb.solve(1, v_coeffs, evecs=True)
        expected = np.arange(10) + 0.5
        self.assertAlmostEqual(energies[0], 0.5, 1)
        self.assertAlmostEqual(energies[1], 1.5, 1)
        self.assertAlmostEqual(energies[2], 2.5, 1)
        self.assertAlmostEqual(energies[3], 3.5, 1)
        self.assertAlmostEqual(energies[4], 4.5, 1)
        self.assertAlmostEqual(energies[5], 5.5, 1)
        self.assertAlmostEqual(energies[6], 6.5, 1)

        import matplotlib.pyplot as pt
        x = np.arange(0.0, a, 0.001)
        pt.clf()
        for i in xrange(10):
            f = hb.eval_fn(x, orbitals[:,i])
            pt.plot(x, f)
        with tmpdir(__name__, 'test_harmonic1') as dn:
            pt.savefig(os.path.join(dn, "harmonic_wavefunctions.png"))
        pt.clf()
        v = hb.eval_fn(x, v_coeffs)
        pt.plot(x, v)
        for energy in energies[:10]:
            pt.axhline(energy)
        pt.xlim(0,a)
        with tmpdir(__name__, 'test_harmonic2') as dn:
            pt.savefig(os.path.join(dn, "harmonic_levels.png"))

    def test_ethane_hindered(self):
        molecule = load_molecule_g03fchk(
            pkg_resources.resource_filename(__name__, "../data/test/ethane/gaussian.fchk"))
        nma = NMA(molecule)
        rot_scan = load_rotscan_g03log(
            pkg_resources.resource_filename(__name__, "../data/test/rotor/gaussian.log"))
        rotor = Rotor(rot_scan, molecule, rotsym=3, even=True, cancel_freq='scan')
        pf = PartFun(nma, [ExtTrans(), ExtRot(6), rotor])
        self.assertAlmostEqual(rotor.cancel_freq/lightspeed*centimeter, 298, 0)
        rotor = Rotor(rot_scan, molecule, rotsym=3, even=True)
        self.assertAlmostEqual(rotor.cancel_freq/lightspeed*centimeter, 314, 0)
        pf = PartFun(nma, [ExtTrans(), ExtRot(6), rotor])
        self.assertArraysAlmostEqual(
            rotor.hb.eval_fn(rot_scan.potential[0], rotor.v_coeffs),
            rot_scan.potential[1] - rot_scan.potential[1].min(), 1e-4
        )
        # reference data from legacy code (Veronique & co)
        self.assertAlmostEqual(rotor.moment/amu, 11.092362911176032, 2)
        self.assertAlmostEqual(rotor.reduced_moment/amu, 5.5461814555880098, 2)
        self.assertAlmostEqual(np.exp(rotor.log_terms(100.0)[1]), 0.12208E+00, 1)
        self.assertAlmostEqual(rotor.heat_capacity_terms(100.0)[1]/(joule/mol/kelvin), 2.567, 0)
        self.assertAlmostEqual(rotor.entropy_terms(100.0)[1]/(joule/mol), 0.766, 0)
        self.assertAlmostEqual(np.exp(rotor.log_terms(800.0)[1]), 0.21108E+01, 1)
        self.assertAlmostEqual(rotor.heat_capacity_terms(800.0)[1]/(joule/mol/kelvin), 6.346, 1)
        self.assertAlmostEqual(rotor.entropy_terms(800.0)[1]/(joule/mol), 14.824, 1)

        with tmpdir(__name__, 'test_ethane_hindered') as dn:
            rotor.plot_levels(os.path.join(dn, "ethane_hindered_levels.png"), 300)
            pf.write_to_file(os.path.join(dn, "ethane_hindered.txt"))
            ta = ThermoAnalysis(pf, [200,300,400,500,600,700,800,900])
            ta.write_to_file(os.path.join(dn, "ethane_hindered_thermo.csv"))

    def test_ethyl_free(self):
        molecule = load_molecule_g03fchk(
            pkg_resources.resource_filename(__name__, "../data/test/ethyl/gaussian.fchk"))
        nma = NMA(molecule)
        dihedral = [5, 1, 0, 2]
        rot_scan = RotScan(dihedral, molecule)
        rotor = Rotor(rot_scan, molecule, rotsym=6, even=True)
        self.assertAlmostEqual(rotor.cancel_freq/lightspeed*centimeter, 141.2, 0)
        pf = PartFun(nma, [ExtTrans(), ExtRot(1), rotor])
        # reference data from legacy code (Veronique & co)
        self.assertAlmostEqual(rotor.reduced_moment/amu, 4.007, 1)
        self.assertAlmostEqual(np.exp(rotor.log_terms(100.0)[1]), 0.6386, 1)
        self.assertAlmostEqual(np.exp(-rotor.log_terms(100.0)[0]), 0.4168, 1)
        self.assertAlmostEqual(np.exp(rotor.log_terms(800.0)[1]), 1.8062, 1)
        self.assertAlmostEqual(np.exp(-rotor.log_terms(800.0)[0]), 3.9273, 1)

        with tmpdir(__name__, 'test_ethyl_free') as dn:
            rotor.plot_levels(os.path.join(dn, "ethyl_free_levels.png"), 300)
            pf.write_to_file(os.path.join(dn, "ethyl_free.txt"))
            ta = ThermoAnalysis(pf, [200,300,400,500,600,700,800,900])
            ta.write_to_file(os.path.join(dn, "ethyl_free_thermo.csv"))

    def test_imoms(self):
        cases = [
            ("caffeine.xyz", 2, 11, [0, 1, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 18, 19, 20, 21, 22, 23], (2598.9923066760343, 11.373318286792710)),
            ("caffeine.xyz", 2, 11, [16, 17, 15], (11.427609412414192, 11.373318286796099)),
            ("caffeine.xyz", 3, 12, [0, 1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14, 15, 16, 17, 21, 22, 23], (3874.9262249281255, 11.489403874005802)),
            ("caffeine.xyz", 3, 12, [20, 18, 19], (11.554047706686680, 11.489402424437612)),
            ("caffeine.xyz", 4, 13, [0, 1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15, 16, 17, 18, 19, 20], (2298.5324430380965, 11.334929532798469)),
            ("caffeine.xyz", 4, 13, [23, 22, 21], (11.394908129049933, 11.334928102722181)),
            ("ethane.xyz", 0, 1, [2, 3, 4], (11.330123438245337, 5.6648361869614678)),
            ("ethane.xyz", 0, 1, [5, 6, 7], (11.330123438245337, 5.6648361869614661)),
            ("glycerol.xyz", 0, 3, [11], (3.0794510843017311, 3.0113070937447430)),
            ("glycerol.xyz", 0, 3, [1, 2, 4, 5, 6, 7, 8, 9, 10, 12, 13], (951.85671473731713, 3.0113074736677845)),
            ("glycerol.xyz", 1, 4, [0, 2, 3, 5, 6, 7, 8, 9, 10, 11, 13], (1072.5177006846639, 2.9828627310014326)),
            ("glycerol.xyz", 1, 4, [12], (3.0988467514954592, 2.9828627310015023)),
            ("glycerol.xyz", 2, 5, [0, 1, 3, 4, 6, 7, 8, 9, 10, 11, 12], (1071.1143603815583, 2.9517497493009159)),
            ("glycerol.xyz", 2, 5, [13], (3.1115917762553726, 2.9517493768918146)),
            ("glycerol.xyz", 3, 4, [0, 2, 5, 6, 9, 10, 11, 13], (370.75539459124985, 61.588976994367783)),
            ("glycerol.xyz", 3, 4, [8, 1, 12, 7], (124.71612061985820, 61.588969223953136)),
            ("glycerol.xyz", 3, 5, [0, 1, 4, 6, 7, 8, 11, 12], (352.35483251604194, 55.690249341790206)),
            ("glycerol.xyz", 3, 5, [9, 2, 10, 13], (116.03080804859955, 55.690242315592506)),
            ("nicotine.xyz", 0, 7, [1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 22, 23, 24, 25], (4653.1199884792477, 11.230510638276883)),
            ("nicotine.xyz", 0, 7, [19, 20, 21], (11.272810801992463, 11.230510638275103)),
            ("nicotine.xyz", 2, 6, [0, 3, 4, 5, 7, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21], (846.76088427036848, 208.26235559926022)),
            ("nicotine.xyz", 2, 6, [1, 8, 9, 10, 11, 22, 23, 24, 25], (307.24370887361914, 208.26235559925976)),
            ("peroxide.xyz", 0, 1, [2], (3.3679370303334704, 1.5747785415198767)),
            ("peroxide.xyz", 0, 1, [3], (3.3679004879866139, 1.5747785415198774)),
        ]
        from molmod.io.xyz import XYZFile
        for fn_xyz, i0, i1, top, expected in cases:
            # preparation
            mol = XYZFile(pkg_resources.resource_filename(
                __name__, os.path.join("../data/test/imom", fn_xyz))).get_molecule()
            masses = np.array([periodic[n].mass for n in mol.numbers])
            masses3 = np.array([masses, masses, masses]).transpose().ravel()
            center = mol.coordinates[i0]
            axis = mol.coordinates[i1] - mol.coordinates[i0]
            axis /= np.linalg.norm(axis)
            # trivial computation of absolute moment
            mom = 0.0
            for i in top:
                delta = mol.coordinates[i] - center
                delta -= axis*np.dot(axis, delta)
                mom += masses[i]*np.linalg.norm(delta)**2
            self.assertAlmostEqual(mom/amu, expected[0], 2)
            # check tamkin routine
            mom, redmom = compute_moments(mol.coordinates, masses3, center, axis, top)
            self.assertAlmostEqual(mom/amu, expected[0], 2)
            self.assertAlmostEqual(redmom/amu, expected[1], 2)

    def test_legacy1(self):
        a = 2*np.pi
        mass = 5.5*amu
        hb = HarmonicBasis(100, a)
        v_coeffs = np.zeros(hb.size, float)
        v_coeffs[0] = 0.5*11.5*kjmol*np.sqrt(a)
        v_coeffs[5] = 0.5*11.5*kjmol*np.sqrt(a/2)
        self.assertArraysAlmostEqual(
            hb.eval_fn(np.array([0, a/6]), v_coeffs),
            np.array([11.5*kjmol, 0.0]),
        )
        energies, orbitals = hb.solve(mass, v_coeffs, evecs=True)

        import matplotlib.pyplot as pt
        x = np.arange(0.0, a, 0.001)
        pt.clf()
        for i in xrange(10):
            f = hb.eval_fn(x, orbitals[:,i])
            pt.plot(x, f+i)
        with tmpdir(__name__, 'test_legacy1') as dn:
            pt.savefig(os.path.join(dn, "legacy_wavefunctions.png"))

        # check energy levels
        expected = np.array([
            1.7635118, 1.76361979, 5.11795465, 8.04553104, 8.1095722,
            10.3876796, 11.8999683, 12.9078395, 14.6739639, 16.7836847,
            19.1722507,
            1.76361979,  5.11795465, 5.12218335, 8.1095722, 10.3876796,
            10.8661504, 12.9078395, 14.6739639, 16.7544718, 19.1722507,
        ])
        expected.sort()
        self.assertArraysAlmostEqual(energies[:10]/kjmol, expected[:10], 1e-3)
        self.assertAlmostEqual(np.exp(-energies/(100*boltzmann)).sum()/3.0, 0.12208E+00, 5)

    def test_load_rotor_margot(self):
        rot_scan = load_rotscan_g03log(
            pkg_resources.resource_filename(__name__, "../data/test/rotor/margot.log"))
        assert rot_scan.potential.shape == (2, 1)
        assert (rot_scan.dihedral == [2, 3, 4, 5]).all()
