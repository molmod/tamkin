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


import numpy as np
import os
import pkg_resources
import unittest

from molmod.periodic import periodic
from molmod.units import angstrom, amu, calorie, avogadro, electronvolt
from molmod.constants import lightspeed
from molmod.test.common import tmpdir

from tamkin import *


__all__ = ["IOTestCase"]


class IOTestCase(unittest.TestCase):
    def test_load_fixed_g03com(self):
        fixed_atoms = load_fixed_g03com(
            pkg_resources.resource_filename(__name__, "../data/test/mat/Zp_p_prod.18aug.com"))
        self.assertEqual(len(fixed_atoms), 48)
        self.assertEqual(fixed_atoms, range(114,114+48))

    def test_load_molecule_g03fchk(self):
        atoms = 181
        molecule = load_molecule_g03fchk(
            pkg_resources.resource_filename(__name__, "../data/test/mat/Zp_p_react.28aug.fchk"))
        self.assertEqual(molecule.hessian.shape,(atoms*3,atoms*3))
        self.assertAlmostEqual(molecule.energy, -3053.805846445570, 7)
        molecule = load_molecule_g03fchk(
            pkg_resources.resource_filename(__name__, "../data/test/mat/Zp_p_react.28aug.fchk"),
            energy=-123.0)
        self.assertAlmostEqual(molecule.energy, -123.0, 7)
        molecule = load_molecule_g03fchk(
            pkg_resources.resource_filename(__name__, "../data/test/mat/Zp_p_react.28aug.fchk"),
            pkg_resources.resource_filename(__name__, "../data/test/mat/Zp_p_react.14mei.fchk"))
        self.assertAlmostEqual(molecule.energy, -18613.135744186180, 7)

    def test_load_molecule_g98fchk(self):
        atoms = 6
        molecule = load_molecule_g98fchk(
            pkg_resources.resource_filename(__name__, "../data/test/g98/freqs.fchk"))
        self.assertEqual(molecule.hessian.shape,(atoms*3,atoms*3))
        self.assertAlmostEqual(molecule.masses[0]/amu, 12.011)
        self.assertAlmostEqual(molecule.masses[2]/amu, 1.0079)
        self.assertAlmostEqual(molecule.energy, -78.58745828877478, 7)
        molecule = load_molecule_g98fchk(
            pkg_resources.resource_filename(__name__, "../data/test/g98/freqs.fchk"),
            energy=-123.0)
        self.assertAlmostEqual(molecule.energy, -123.0, 7)

    def test_load_molecule_g03fchkvdw(self):
        atoms = 179
        molecule = load_molecule_g03fchk(
            pkg_resources.resource_filename(__name__, "../data/test/matvdw/R.fchk"),
            pkg_resources.resource_filename(__name__, "../data/test/matvdw/R_SCF.fchk"),
            pkg_resources.resource_filename(__name__, "../data/test/matvdw/R_b3lyp-d.out"))

        self.assertEqual(molecule.hessian.shape,(atoms*3,atoms*3))
        self.assertAlmostEqual(molecule.energy,-18612.352569964281 , 7)

    def test_load_fixed_fchk(self):
        molecule = load_molecule_g03fchk(
            pkg_resources.resource_filename(__name__, "../data/test/ethane/gaussian.fchk"))
        assert (molecule.fixed == [0, 6, 7]).all()

    def test_load_molecule_cp2k(self):
        molecule = load_molecule_cp2k(
            pkg_resources.resource_filename(__name__, "../data/test/cp2k/pentane/sp.out"),
            pkg_resources.resource_filename(__name__, "../data/test/cp2k/pentane/freq.out"))
        self.assertAlmostEqual(molecule.energy, 0.012255059530862)
        self.assertEqual(molecule.multiplicity, 1)
        self.assertEqual(molecule.numbers[0], 6)
        self.assertEqual(molecule.numbers[4], 1)
        self.assertAlmostEqual(molecule.masses[0], periodic[6].mass, 5)
        self.assertAlmostEqual(molecule.masses[4], periodic[1].mass, 0)
        self.assertAlmostEqual(molecule.coordinates[5,1]/angstrom, 13.928520)
        self.assertAlmostEqual(molecule.gradient[0,2], 0.0000000038, 9)
        self.assertAlmostEqual(molecule.gradient[11,0], 0.0000000177, 9)
        self.assertAlmostEqual(molecule.hessian[0,0], 49.629809*1e-6*molecule.masses[0], 6)
        self.assertAlmostEqual(molecule.hessian[-1,-1], 287.884198*1e-6*molecule.masses[-1], 6)
        self.assertAlmostEqual(molecule.unit_cell.matrix[0,0]/angstrom, 30.000,3)
        self.assertAlmostEqual(molecule.unit_cell.matrix[1,2]/angstrom, 0.000,3)

    def test_load_molecule_cp2k_23(self):
        molecule = load_molecule_cp2k(
            pkg_resources.resource_filename(__name__, "../data/test/cp2k/john/scf.out"),
            pkg_resources.resource_filename(__name__, "../data/test/cp2k/john/quickie.out"))
        self.assertAlmostEqual(molecule.energy, -141.189006820072791)
        self.assertEqual(molecule.multiplicity, 1)
        self.assertEqual(molecule.numbers[0], 14)
        self.assertEqual(molecule.numbers[4], 8)
        self.assertAlmostEqual(molecule.masses[0], periodic[14].mass)
        self.assertAlmostEqual(molecule.masses[4], periodic[8].mass)
        self.assertAlmostEqual(molecule.coordinates[5,1]/angstrom, 12.150867)
        self.assertAlmostEqual(molecule.gradient[0,2], -0.00008196, 9)
        self.assertAlmostEqual(molecule.gradient[11,0], -0.00041872, 9)
        self.assertAlmostEqual(molecule.hessian[0,0], 8.297378*1e-6*molecule.masses[0], 6)
        self.assertAlmostEqual(molecule.hessian[-1,-1], 71.133554*1e-6*molecule.masses[-1], 6)
        self.assertAlmostEqual(molecule.unit_cell.matrix[0,0]/angstrom, 20.000,3)
        self.assertAlmostEqual(molecule.unit_cell.matrix[1,2]/angstrom, 0.000,3)

    def test_load_molecule_cp2k_dan_phva(self):
        molecule = load_molecule_cp2k(
            pkg_resources.resource_filename(__name__, "../data/test/cp2k/dan/cp2k.out"),
            pkg_resources.resource_filename(__name__, "../data/test/cp2k/dan/freq.out"))
        self.assertAlmostEqual(molecule.energy, -290.409097595743333)
        self.assertEqual(molecule.multiplicity, 1)
        self.assertEqual(molecule.numbers[0], 6)
        self.assertEqual(molecule.numbers[4], 6)
        self.assertEqual(molecule.numbers[-1], 1)
        self.assertAlmostEqual(molecule.masses[0], periodic[6].mass, 0)
        self.assertAlmostEqual(molecule.masses[-1], periodic[1].mass, 0)
        self.assertAlmostEqual(molecule.coordinates[5,1], -2.458266*angstrom)
        self.assertAlmostEqual(molecule.gradient[-1,0], 0.00000044, 9)
        self.assertAlmostEqual(molecule.gradient[-3,2], 0.00000080, 9)
        self.assertAlmostEqual(molecule.hessian[0,0], 0.0, 6)
        self.assertAlmostEqual(molecule.hessian[0,-1], 0.0, 6)
        self.assertAlmostEqual(molecule.hessian[-1,0], 0.0, 6)
        self.assertAlmostEqual(molecule.hessian[-18,-1],
            0.5*(-0.418336-0.418196)*1e-6*(molecule.masses[-1]*molecule.masses[-18])**0.5, 6)
        self.assertAlmostEqual(molecule.hessian[-18,-18], 19.137757*1e-6*molecule.masses[-18], 6)
        self.assertAlmostEqual(molecule.hessian[-1,-1], 51.483068*1e-6*molecule.masses[-1], 6)
        self.assertAlmostEqual(molecule.unit_cell.matrix[0,0], 10.657*angstrom, 3)
        self.assertAlmostEqual(molecule.unit_cell.matrix[1,0], -6.153*angstrom, 3)

    def test_load_fixed_cp2k(self):
        fixed = load_fixed_cp2k(
            pkg_resources.resource_filename(__name__, "../data/test/cp2k/dan/freq.out"))
        np.testing.assert_equal(fixed, np.arange(52 - 6))

    def test_load_molecule_cpmd(self):
        molecule = load_molecule_cpmd(
            pkg_resources.resource_filename(__name__, "../data/test/cpmd/damp.out"),
            pkg_resources.resource_filename(__name__, "../data/test/cpmd/GEOMETRY.xyz"),
            pkg_resources.resource_filename(__name__, "../data/test/cpmd/MOLVIB"))
        self.assertAlmostEqual(molecule.energy, -17.14142079)
        self.assertEqual(molecule.multiplicity, 1)
        self.assertEqual(molecule.numbers[0], 8)
        self.assertEqual(molecule.numbers[2], 1)
        self.assertAlmostEqual(molecule.masses[0]/amu, 15.999400)
        self.assertAlmostEqual(molecule.masses[2]/amu, 1.007970)
        self.assertAlmostEqual(molecule.coordinates[2,0]/angstrom, .907207019799)
        self.assertAlmostEqual(molecule.gradient[0,2], 0.0)
        self.assertAlmostEqual(molecule.gradient[2,0], 0.0)
        self.assertAlmostEqual(molecule.hessian[0,0], 0.530679165332463953, 6)
        self.assertAlmostEqual(molecule.hessian[-1,-1], 0.921045226686428159E-01, 6)

    def test_load_molecule_charmm(self):
        molecule = load_molecule_charmm(
            pkg_resources.resource_filename(__name__, "../data/test/an/ethanol.cor"),
            pkg_resources.resource_filename(__name__, "../data/test/an/ethanol.hess.full"))
        self.assertAlmostEqual(molecule.energy/(1000*calorie/avogadro), -2.1303308955)
        self.assertEqual(molecule.multiplicity, 1)
        self.assertEqual(molecule.numbers[0], 6)
        self.assertEqual(molecule.numbers[4], 1)
        self.assertAlmostEqual(molecule.masses[0]/amu, 12.01100)
        self.assertAlmostEqual(molecule.masses[4]/amu, 1.00800)
        self.assertAlmostEqual(molecule.coordinates[5,1]/angstrom, 1.3582528196)
        self.assertAlmostEqual(molecule.gradient[0,2]/(1000*calorie/avogadro/angstrom), -0.0000000007, 9)
        self.assertAlmostEqual(molecule.gradient[8,0]/(1000*calorie/avogadro/angstrom), -0.0000001462, 9)
        self.assertAlmostEqual(molecule.hessian[0,0]/(1000*calorie/avogadro /angstrom**2), 1409.7091337384, 6)
        self.assertAlmostEqual(molecule.hessian[-1,-1]/(1000*calorie/avogadro /angstrom**2), 474.7950312957, 6)

    def test_load_molecule_qchem(self):
        molecule = load_molecule_qchem(
            pkg_resources.resource_filename(__name__, "../data/test/qchem/h2o2.hf.sto-3g.freq.out"),
            hessfile=pkg_resources.resource_filename(__name__, "../data/test/qchem/hessian.dat"))
        self.assertAlmostEqual(molecule.energy, -148.7649966058)
        self.assertEqual(molecule.multiplicity, 1)
        self.assertEqual(molecule.numbers[0], 1)
        self.assertEqual(molecule.numbers[3], 8)
        self.assertAlmostEqual(molecule.masses[0]/amu,  1.00783, 5)
        self.assertAlmostEqual(molecule.masses[3]/amu, 15.99491, 5)
        self.assertAlmostEqual(molecule.coordinates[2,1]/angstrom, -0.688720)
        self.assertAlmostEqual(molecule.gradient[0,2]/(1000*calorie/avogadro/angstrom), -0.00, 5)
        self.assertAlmostEqual(molecule.gradient[3,0]/(1000*calorie/avogadro/angstrom), -0.00, 5)
        self.assertAlmostEqual(molecule.hessian[0,0]/(1000*calorie/avogadro/angstrom**2), 364.769480916757800060, 6)
        self.assertAlmostEqual(molecule.hessian[-1,-1]/(1000*calorie/avogadro/angstrom**2), 338.870127396983150447, 6)

    def test_load_molecule_vasp_53(self):
        molecule = load_molecule_vasp(
            pkg_resources.resource_filename(__name__, "../data/test/lucas/vasp_5_3_5_complex/CONTCAR_opt"),
            pkg_resources.resource_filename(__name__, "../data/test/lucas/vasp_5_3_5_complex/OUTCAR_freq"))
        # contcar
        assert molecule.numbers[0] == 6
        assert (molecule.numbers[1:] == 1).all()
        assert molecule.size == 5
        assert molecule.unit_cell.matrix[0,0] == 15.0*angstrom
        assert molecule.unit_cell.matrix[1,2] == 0.0
        self.assertAlmostEqual(molecule.coordinates[0,0]/angstrom, 7.15840, 3)
        self.assertAlmostEqual(molecule.coordinates[1,2]/angstrom, 8.44640, 2) #?
        self.assertAlmostEqual(molecule.coordinates[-1,-1]/angstrom, 6.95131, 2) #?
        # outcar_freq
        assert molecule.masses[0] == 12.011*amu
        assert (molecule.masses[1:] == 1.000*amu).all()
        hunit = electronvolt/angstrom**2
        assert molecule.hessian[0,0] == 53.624756*hunit
        assert molecule.hessian[-1,-1] == 31.299419*hunit
        self.assertAlmostEqual(molecule.hessian[2,5], 0.5*(-7.551817 + 3.319877)*hunit)
        assert molecule.energy == -24.11901936*electronvolt
        gunit = electronvolt/angstrom
        assert molecule.gradient[0, 0] == 0.096977*gunit
        assert molecule.gradient[2, 1] == 0.100275*gunit
        assert molecule.gradient[-1, -1] == -0.212810*gunit

    def test_load_molecule_vasp_5_3_5_gamma(self):
        molecule = load_molecule_vasp(
            pkg_resources.resource_filename(__name__, "../data/test/julianna/vasp_5_3_5_gamma/CONTCAR_opt"),
            pkg_resources.resource_filename(__name__, "../data/test/julianna/vasp_5_3_5_gamma/OUTCAR_freq"))
        # contcar
        assert molecule.numbers[0] == 6
        assert (molecule.numbers[1:] == 1).all()
        assert molecule.size == 5
        assert molecule.unit_cell.matrix[0,0] == 15.0*angstrom
        assert molecule.unit_cell.matrix[1,2] == 0.0
        self.assertAlmostEqual(molecule.coordinates[0,0]/angstrom, 7.15782, 3)
        self.assertAlmostEqual(molecule.coordinates[1,2]/angstrom, 8.44278, 1) #??
        self.assertAlmostEqual(molecule.coordinates[-1,-1]/angstrom, 6.95393, 2)  #?
        # outcar_freq
        assert molecule.masses[0] == 12.011*amu
        assert (molecule.masses[1:] == 1.000*amu).all()
        hunit = electronvolt/angstrom**2
        assert molecule.hessian[0,0] == 47.756815*hunit
        assert molecule.hessian[-1,-1] == 31.561376*hunit
        self.assertAlmostEqual(molecule.hessian[2,5], 0.5*(-2.265871 + -3.645039)*hunit)
        assert molecule.energy == -24.12364199*electronvolt
        gunit = electronvolt/angstrom
        assert molecule.gradient[0, 0] == -0.005459*gunit
        assert molecule.gradient[2, 1] == -0.008215*gunit
        assert molecule.gradient[-1, -1] == 0.003424*gunit

    def test_load_molecule_vasp_5_3_5_gamma_part(self):
        molecule = load_molecule_vasp(
            pkg_resources.resource_filename(__name__, "../data/test/julianna/vasp_5_3_5_gamma/CONTCAR_opt"),
            pkg_resources.resource_filename(__name__, "../data/test/julianna/vasp_5_3_5_gamma/OUTCAR_freq_part"))
        # outcar_freq
        hunit = electronvolt/angstrom**2
        assert molecule.hessian[0,0] == 0.0
        assert molecule.hessian[2,5] == 0.0
        assert molecule.hessian[-1,-1] == 31.561374*hunit
        self.assertAlmostEqual(molecule.hessian[6,9], 0.5*(-2.601094 + -2.794160)*hunit)

    def test_load_molecule_vasp_5_2_11_complex(self):
        molecule = load_molecule_vasp(
            pkg_resources.resource_filename(__name__, "../data/test/julianna/vasp_5_2_11_complex/CONTCAR_opt"),
            pkg_resources.resource_filename(__name__, "../data/test/julianna/vasp_5_2_11_complex/OUTCAR_freq"))
        # outcar_freq
        assert molecule.masses[0] == 12.011*amu
        assert (molecule.masses[1:] == 1.000*amu).all()
        hunit = electronvolt/angstrom**2
        assert molecule.hessian[0,0] == 47.762604*hunit
        assert molecule.hessian[-1,-1] == 31.565279*hunit
        self.assertAlmostEqual(molecule.hessian[2,5], 0.5*(-3.648356 + -2.264335)*hunit)
        assert molecule.energy == -24.123642*electronvolt
        gunit = electronvolt/angstrom
        assert molecule.gradient[0, 0] == -0.005432*gunit
        assert molecule.gradient[2, 1] == -0.008204*gunit
        assert molecule.gradient[-1, -1] == 0.003425*gunit

    def test_load_molecule_vasp_5_2_11_complex_part(self):
        molecule = load_molecule_vasp(
            pkg_resources.resource_filename(__name__, "../data/test/julianna/vasp_5_2_11_complex/CONTCAR_opt"),
            pkg_resources.resource_filename(__name__, "../data/test/julianna/vasp_5_2_11_complex/OUTCAR_freq_part"))
        # outcar_freq
        hunit = electronvolt/angstrom**2
        assert molecule.hessian[0,0] == 0.0
        assert molecule.hessian[2,5] == 0.0
        assert molecule.hessian[-1,-1] == 31.565272*hunit
        self.assertAlmostEqual(molecule.hessian[6,9], 0.5*(-2.600373 + -2.836454)*hunit)
        assert molecule.energy == -24.123642*electronvolt

    def test_load_molecule_vasp_5_3_3_complex(self):
        molecule = load_molecule_vasp(
            pkg_resources.resource_filename(__name__, "../data/test/julianna/vasp_5_3_3_complex/CONTCAR_opt"),
            pkg_resources.resource_filename(__name__, "../data/test/julianna/vasp_5_3_3_complex/OUTCAR_freq"))
        # outcar_freq
        assert molecule.masses[0] == 12.011*amu
        assert (molecule.masses[1:] == 1.000*amu).all()
        hunit = electronvolt/angstrom**2
        assert molecule.hessian[0,0] == 47.757096*hunit
        assert molecule.hessian[-1,-1] == 31.561341*hunit
        self.assertAlmostEqual(molecule.hessian[2,5], 0.5*(-3.645047 + -2.265767)*hunit)
        assert molecule.energy == -24.12364179*electronvolt
        gunit = electronvolt/angstrom
        assert molecule.gradient[0, 0] == -0.005431*gunit
        assert molecule.gradient[2, 1] == -0.008279*gunit
        assert molecule.gradient[-1, -1] == 0.003335*gunit

    def test_load_molecule_vasp_5_3_3_complex_part(self):
        molecule = load_molecule_vasp(
            pkg_resources.resource_filename(__name__, "../data/test/julianna/vasp_5_3_3_complex/CONTCAR_opt"),
            pkg_resources.resource_filename(__name__, "../data/test/julianna/vasp_5_3_3_complex/OUTCAR_freq_part"))
        # outcar_freq
        hunit = electronvolt/angstrom**2
        assert molecule.hessian[0,0] == 0.0
        assert molecule.hessian[2,5] == 0.0
        assert molecule.hessian[-1,-1] == 31.561353*hunit
        self.assertAlmostEqual(molecule.hessian[6,9], 0.5*(-2.601060 + -2.794022)*hunit)

    def test_load_molecule_vasp_5_3_5_gamma_part_energy(self):
        molecule = load_molecule_vasp(
            pkg_resources.resource_filename(__name__, "../data/test/julianna/vasp_5_3_5_gamma/CONTCAR_opt"),
            pkg_resources.resource_filename(__name__, "../data/test/julianna/vasp_5_3_5_gamma/OUTCAR_freq_part"),
            energy=1.476)
        assert molecule.energy == 1.476

    def test_load_molecule_vasp_5_3_5_gamma_part_outcar_energy(self):
        molecule = load_molecule_vasp(
            pkg_resources.resource_filename(__name__, "../data/test/julianna/vasp_5_3_5_gamma/CONTCAR_opt"),
            pkg_resources.resource_filename(__name__, "../data/test/julianna/vasp_5_3_5_gamma/OUTCAR_freq_part"),
            pkg_resources.resource_filename(__name__, "../data/test/lucas/vasp_5_3_5_complex/OUTCAR_freq"))
        assert molecule.energy == -24.11901936*electronvolt

    def test_checkpoint(self):
        molecule = load_molecule_cp2k(
            pkg_resources.resource_filename(__name__, "../data/test/cp2k/pentane/sp.out"),
            pkg_resources.resource_filename(__name__, "../data/test/cp2k/pentane/freq.out"))
        nma1 = NMA(molecule)
        with tmpdir(__name__, 'test_checkpoint') as dn:
            fn_out = os.path.join(dn, 'test.chk')
            nma1.write_to_file(fn_out)
            nma2 = NMA.read_from_file(fn_out)

        self.assertEqual(nma1.freqs.shape, nma2.freqs.shape)
        self.assertEqual(nma1.modes.shape, nma2.modes.shape)
        self.assertEqual(nma1.masses.shape, nma2.masses.shape)
        self.assertEqual(nma1.numbers.shape, nma2.numbers.shape)
        self.assertEqual(nma1.coordinates.shape, nma2.coordinates.shape)
        self.assertEqual(nma1.inertia_tensor.shape, nma2.inertia_tensor.shape)

        self.assert_(abs(nma1.freqs - nma2.freqs).max()/abs(nma1.freqs).max() < 1e-15)
        self.assert_(abs(nma1.modes - nma2.modes).max()/abs(nma1.modes).max() < 1e-15)
        self.assert_(abs(nma1.masses - nma2.masses).max()/abs(nma1.masses).max() < 1e-15)
        self.assert_(abs(nma1.coordinates - nma2.coordinates).max()/abs(nma1.coordinates).max() < 1e-15)
        self.assert_(abs(nma1.inertia_tensor - nma2.inertia_tensor).max()/abs(nma1.inertia_tensor).max() < 1e-15)
        self.assert_((nma1.numbers==nma2.numbers).all())

        self.assertAlmostEqual(nma1.mass, nma2.mass)
        self.assertAlmostEqual(nma1.energy, nma2.energy)
        self.assertEqual(nma1.multiplicity, nma2.multiplicity)
        self.assertEqual(nma1.symmetry_number, nma2.symmetry_number)

    def test_load_indices(self):
        blocks = load_indices(pkg_resources.resource_filename(__name__, "../data/test/an/fixed.07.txt"), groups=True)
        self.assertEqual(blocks, [[3,2,6]])
        blocks = load_indices(pkg_resources.resource_filename(__name__, "../data/test/an/fixed.07.txt"))
        self.assertEqual(blocks, [3,2,6])
        blocks = load_indices(pkg_resources.resource_filename(__name__, "../data/test/an/fixed.08.txt"))
        self.assertEqual(blocks, [5,4,8])

    def test_load_indices_ranges(self):
        blocks = load_indices(pkg_resources.resource_filename(__name__, "../data/test/an/fixed_ranges.txt"), groups=True)
        self.assertEqual(blocks, [[0, 2, 3, 4], [9, 10, 11, 12, 13, 19], [21]])
        blocks = load_indices(pkg_resources.resource_filename(__name__, "../data/test/an/fixed_ranges.txt"), shift=0)
        self.assertEqual(blocks, [1, 3, 4, 5, 10, 11, 12, 13, 14, 20, 22])
        blocks = load_indices(pkg_resources.resource_filename(__name__, "../data/test/an/fixed_ranges_shift.txt"), groups=True)
        self.assertEqual(blocks, [[1, 3, 4, 5], [10, 11, 12, 13, 14, 20], [22]])
        blocks = load_indices(pkg_resources.resource_filename(__name__, "../data/test/an/fixed_ranges_shift.txt"), shift=0)
        self.assertEqual(blocks, [1, 3, 4, 5, 10, 11, 12, 13, 14, 20, 22])

    def test_dump_modes_xyz(self):
        molecule = load_molecule_charmm(
            pkg_resources.resource_filename(__name__, "../data/test/an/ethanol.cor"),
            pkg_resources.resource_filename(__name__, "../data/test/an/ethanol.hess.full"))
        nma = NMA(molecule)
        with tmpdir(__name__, 'test_dump_modes_xyz') as dn:
            prefix = os.path.join(dn, 'mode')
            dump_modes_xyz(nma, 6, prefix=prefix, amplitude=50.0)
            with open(prefix + ".6.xyz") as f:
                # 1st line
                line = f.readline().strip()
                self.assertEqual(line,"9")
                # 2nd line
                line = f.readline().strip()
                self.assertEqual(line,"frame 0")
                # 3rd line
                line = f.readline()
                words = line.split()
                self.assertEqual(len(words),4)
                self.assertEqual(words[0],"C")
                self.assertEqual(words[2],"0.081608346")
                for i in range(9):
                    line = f.readline().strip()
                self.assertEqual(line,"9")

    def test_dump_modes_gaussian(self):
        molecule = load_molecule_charmm(
            pkg_resources.resource_filename(__name__, "../data/test/an/ethanol.cor"),
            pkg_resources.resource_filename(__name__, "../data/test/an/ethanol.hess.full"))
        nma = NMA(molecule)
        with tmpdir(__name__, 'test_dump_modes_gaussian') as dn:
            # blind test, no double checking of the output file
            fn_log = os.path.join(dn, 'modes_gaussian.log')
            dump_modes_gaussian(fn_log, nma)

    def test_load_dump_indices1(self):
        subs = range(10)
        with tmpdir(__name__, 'test_load_dump_indices1') as dn:
            dump_indices(os.path.join(dn, "subs-atoms.1.txt"), subs, shift=0)
            dump_indices(os.path.join(dn, "subs-atoms.2.txt"), subs, shift=1)
            dump_indices(os.path.join(dn, "subs-atoms.3.txt"), subs)

            subs1 = load_indices(os.path.join(dn, "subs-atoms.1.txt"), shift=0)
            subs2 = load_indices(os.path.join(dn, "subs-atoms.2.txt"), shift=-1)
            subs22 = load_indices(os.path.join(dn, "subs-atoms.2.txt"))  # should not matter
            subs3 = load_indices(os.path.join(dn, "subs-atoms.3.txt"))

            blocks = [range(10), range(10,20)]
            dump_indices(os.path.join(dn, "blocks.1.txt"), blocks, shift=0)
            dump_indices(os.path.join(dn, "blocks.2.txt"), blocks, shift=1)
            dump_indices(os.path.join(dn, "blocks.3.txt"), blocks)

            blocks1 = load_indices(os.path.join(dn, "blocks.1.txt"), shift=0, groups=True)
            blocks2 = load_indices(os.path.join(dn, "blocks.2.txt"), shift=-1, groups=True)
            blocks22 = load_indices(os.path.join(dn, "blocks.2.txt"), groups=True)  # should not matter
            blocks3 = load_indices(os.path.join(dn, "blocks.3.txt"), groups=True)

        self.assertEqual(len(subs),len(subs1))
        for (i,j) in zip(subs,subs1):
            self.assertEqual(i,j)
        self.assertEqual(len(subs),len(subs2))
        for i,j in zip(subs,subs2):
            self.assertEqual(i,j)
        self.assertEqual(len(subs),len(subs22))
        for i,j in zip(subs,subs22):
            self.assertEqual(i,j)
        self.assertEqual(len(subs),len(subs3))
        for i,j in zip(subs,subs3):
            self.assertEqual(i,j)

        self.assertEqual(len(blocks),len(blocks1))
        for bl,bl1 in zip(blocks,blocks1):
            for i,j in zip(bl,bl1):
                self.assertEqual(i,j)
        self.assertEqual(len(blocks),len(blocks2))
        for bl,bl1 in zip(blocks,blocks2):
            for i,j in zip(bl,bl1):
                self.assertEqual(i,j)
        self.assertEqual(len(blocks),len(blocks22))
        for bl,bl1 in zip(blocks,blocks22):
            for i,j in zip(bl,bl1):
                self.assertEqual(i,j)
        self.assertEqual(len(blocks),len(blocks3))
        for bl,bl1 in zip(blocks,blocks3):
            for i,j in zip(bl,bl1):
                self.assertEqual(i,j)

    def test_load_dump_indices2(self):
        randint = np.random.randint
        for counter in xrange(20):
            for compact in True, False:
                shift = randint(-5,5)
                indices = []
                for i in xrange(randint(10,20)):
                    l = list(set(randint(0,10,randint(20))))
                    if len(l) > 0:
                        indices.append(l)
                indices_flat = sum(indices, [])

                with tmpdir(__name__, '{}_{}'.format('test_load_dump_indices2', counter)) as dn:
                    dump_indices(os.path.join(dn, "indices_blocks.txt"), indices, compact=compact, shift=shift)
                    dump_indices(os.path.join(dn, "indices_flat.txt"), indices_flat, compact=compact, shift=shift)

                    check = load_indices(os.path.join(dn, "indices_blocks.txt"), shift=-shift, groups=True)
                    self.assertEqual(indices, check)
                    check = load_indices(os.path.join(dn, "indices_blocks.txt"), shift=-shift)
                    self.assertEqual(indices_flat, check)
                    check = load_indices(os.path.join(dn, "indices_flat.txt"), shift=-shift, groups=True)
                    self.assertEqual([indices_flat], check)
                    check = load_indices(os.path.join(dn, "indices_flat.txt"), shift=-shift)
                    self.assertEqual(indices_flat, check)

    def test_punch(self):
        mol0 = load_molecule_g03fchk(
            pkg_resources.resource_filename(__name__, "../data/test/punch/gaussian.fchk"))
        mol1 = load_molecule_g03fchk(
            pkg_resources.resource_filename(__name__, "../data/test/punch/gaussian.fchk"),
            fn_punch=pkg_resources.resource_filename(__name__, "../data/test/punch/fort.7"))
        assert abs(mol0.gradient - mol1.gradient).max() < 1e-8
        assert abs(mol0.hessian - mol1.hessian).max() < 1e-8

    def test_dftd3(self):
        assert load_dftd3(
            pkg_resources.resource_filename(__name__, "../data/test/dftd3/dftd3.out")) == -0.00330057

    def test_dftd_orca(self):
        assert load_dftd_orca(
            pkg_resources.resource_filename(__name__, "../data/test/matvdw/R_b3lyp-d.out")) == -0.404083275
