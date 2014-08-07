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


from tamkin import *

from molmod.periodic import periodic
from molmod.units import angstrom, amu, calorie, avogadro, electronvolt
from molmod.constants import lightspeed

import unittest, numpy


__all__ = ["IOTestCase"]


class IOTestCase(unittest.TestCase):
    def test_load_fixed_g03com(self):
        fixed_atoms = load_fixed_g03com("test/input/mat/Zp_p_prod.18aug.com")
        self.assertEqual(len(fixed_atoms), 48)
        self.assertEqual(fixed_atoms, range(114,114+48))

    def test_load_molecule_g03fchk(self):
        atoms = 181
        molecule = load_molecule_g03fchk("test/input/mat/Zp_p_react.28aug.fchk")
        self.assertEqual(molecule.hessian.shape,(atoms*3,atoms*3))
        self.assertAlmostEqual(molecule.energy, -3053.805846445570, 7)
        molecule = load_molecule_g03fchk("test/input/mat/Zp_p_react.28aug.fchk", energy=-123.0)
        self.assertAlmostEqual(molecule.energy, -123.0, 7)
        molecule = load_molecule_g03fchk("test/input/mat/Zp_p_react.28aug.fchk", "test/input/mat/Zp_p_react.14mei.fchk")
        self.assertAlmostEqual(molecule.energy, -18613.135744186180, 7)

    def test_load_molecule_g98fchk(self):
        atoms = 6
        molecule = load_molecule_g98fchk("test/input/g98/freqs.fchk")
        self.assertEqual(molecule.hessian.shape,(atoms*3,atoms*3))
        self.assertAlmostEqual(molecule.masses[0]/amu, 12.011)
        self.assertAlmostEqual(molecule.masses[2]/amu, 1.0079)
        self.assertAlmostEqual(molecule.energy, -78.58745828877478, 7)
        molecule = load_molecule_g98fchk("test/input/g98/freqs.fchk", energy=-123.0)
        self.assertAlmostEqual(molecule.energy, -123.0, 7)

    def test_load_molecule_g03fchkvdw(self):
        atoms = 179
        molecule = load_molecule_g03fchk("test/input/matvdw/R.fchk","test/input/matvdw/R_SCF.fchk","test/input/matvdw/R_b3lyp-d.out")

        self.assertEqual(molecule.hessian.shape,(atoms*3,atoms*3))
        self.assertAlmostEqual(molecule.energy,-18612.352569964281 , 7)

    def test_load_molecule_cp2k(self):
        molecule = load_molecule_cp2k("test/input/cp2k/pentane/sp.out", "test/input/cp2k/pentane/freq.out")
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
        molecule = load_molecule_cp2k("test/input/cp2k/john/scf.out", "test/input/cp2k/john/quickie.out")
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

    def test_load_molecule_cpmd(self):
        molecule = load_molecule_cpmd("test/input/cpmd/damp.out", "test/input/cpmd/GEOMETRY.xyz", "test/input/cpmd/MOLVIB")
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
        molecule = load_molecule_charmm("test/input/an/ethanol.cor","test/input/an/ethanol.hess.full")
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
        molecule = load_molecule_qchem("test/input/qchem/h2o2.hf.sto-3g.freq.out", hessfile="test/input/qchem/hessian.dat")
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

    def test_load_molecule_vasp(self):
        molecule = load_molecule_vasp("test/input/vasp/xyz-structure","test/input/vasp/OUTCAR")
        self.assertEqual(molecule.numbers[0],14)
        self.assertEqual(molecule.numbers[107],13)
        self.assertAlmostEqual(molecule.masses[0]/amu, 28.085)
        self.assertAlmostEqual(molecule.masses[120]/amu, 1.000)
        self.assertAlmostEqual(molecule.coordinates[5,1]/angstrom, 2.93027 )
        self.assertAlmostEqual(molecule.gradient[0,2]/(electronvolt/angstrom), 0.016534, 5)
        self.assertAlmostEqual(molecule.gradient[8,0]/(electronvolt/angstrom), -0.035937, 5)
        self.assertAlmostEqual( - molecule.hessian[0,0]/(electronvolt/angstrom**2), -46.644731, 6)
        self.assertAlmostEqual( - molecule.hessian[-1,-1]/(electronvolt/angstrom**2), -5.524062, 6)
        self.assertAlmostEqual(molecule.unit_cell.matrix[0,0]/angstrom, 11.329193060, 5)
        self.assertAlmostEqual(molecule.unit_cell.matrix[1,2]/angstrom, -0.017392342, 5)
        # if VASP contains only a partial Hessian
        molecule = load_molecule_vasp("test/input/vasp/xyz-structure-part","test/input/vasp/OUTCAR-part")
        self.assertEqual(molecule.numbers[0],14)
        self.assertEqual(molecule.numbers[107],13)
        self.assertAlmostEqual(molecule.masses[0]/amu, 28.085)
        self.assertAlmostEqual(molecule.masses[120]/amu, 1.000)
        self.assertAlmostEqual(molecule.coordinates[5,1]/angstrom, 2.93027 )
        self.assertAlmostEqual(molecule.gradient[0,2]/(electronvolt/angstrom), 0.016923, 5)
        self.assertAlmostEqual(molecule.gradient[8,0]/(electronvolt/angstrom), -0.036518, 5)
        self.assertAlmostEqual( - molecule.hessian[0,0]/(electronvolt/angstrom**2), -46.646216, 6)
        self.assertAlmostEqual( - molecule.hessian[-1,-1]/(electronvolt/angstrom**2), -5.524077, 6)
        fixed = load_fixed_vasp("test/input/vasp/OUTCAR-part")
        self.assertEqual(fixed[0],2)
        self.assertEqual(fixed[30],53)
        self.assertAlmostEqual(molecule.unit_cell.matrix[0,0]/angstrom, 11.329193060, 5)
        self.assertAlmostEqual(molecule.unit_cell.matrix[1,2]/angstrom, -0.017392342, 5)

    def test_checkpoint(self):
        molecule = load_molecule_cp2k("test/input/cp2k/pentane/sp.out", "test/input/cp2k/pentane/freq.out")
        nma1 = NMA(molecule)
        nma1.write_to_file("test/output/test.chk")
        nma2 = NMA.read_from_file("test/output/test.chk")

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
        blocks = load_indices("test/input/an/fixed.07.txt", groups=True)
        self.assertEqual(blocks, [[3,2,6]])
        blocks = load_indices("test/input/an/fixed.07.txt")
        self.assertEqual(blocks, [3,2,6])

    def test_dump_modes_xyz(self):
        molecule = load_molecule_charmm("test/input/an/ethanol.cor","test/input/an/ethanol.hess.full")
        nma = NMA(molecule)
        dump_modes_xyz(nma, 6, prefix="test/output/mode", amplitude=50.0)
        f = file("test/output/mode.6.xyz")
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
        f.close()

    def test_dump_modes_gaussian(self):
        molecule = load_molecule_charmm("test/input/an/ethanol.cor","test/input/an/ethanol.hess.full")
        nma = NMA(molecule)
        dump_modes_gaussian("test/output/modes_gaussian.log", nma)
        # blind test, no double checking of the output file

    def test_load_dump_indices1(self):
        subs = range(10)
        dump_indices("test/output/subs-atoms.1.txt", subs, shift=0)
        dump_indices("test/output/subs-atoms.2.txt", subs, shift=1)
        dump_indices("test/output/subs-atoms.3.txt", subs)

        subs1 = load_indices("test/output/subs-atoms.1.txt", shift=0)
        self.assertEqual(len(subs),len(subs1))
        for (i,j) in zip(subs,subs1):
            self.assertEqual(i,j)
        subs2 = load_indices("test/output/subs-atoms.2.txt", shift=-1)
        self.assertEqual(len(subs),len(subs2))
        for i,j in zip(subs,subs2):
            self.assertEqual(i,j)
        subs22 = load_indices("test/output/subs-atoms.2.txt")  # should not matter
        self.assertEqual(len(subs),len(subs22))
        for i,j in zip(subs,subs22):
            self.assertEqual(i,j)
        subs3 = load_indices("test/output/subs-atoms.3.txt")
        self.assertEqual(len(subs),len(subs3))
        for i,j in zip(subs,subs3):
            self.assertEqual(i,j)

        blocks = [range(10), range(10,20)]
        dump_indices("test/output/blocks.1.txt", blocks, shift=0)
        dump_indices("test/output/blocks.2.txt", blocks, shift=1)
        dump_indices("test/output/blocks.3.txt", blocks)

        blocks1 = load_indices("test/output/blocks.1.txt", shift=0, groups=True)
        self.assertEqual(len(blocks),len(blocks1))
        for bl,bl1 in zip(blocks,blocks1):
            for i,j in zip(bl,bl1):
                self.assertEqual(i,j)
        blocks2 = load_indices("test/output/blocks.2.txt", shift=-1, groups=True)
        self.assertEqual(len(blocks),len(blocks2))
        for bl,bl1 in zip(blocks,blocks2):
            for i,j in zip(bl,bl1):
                self.assertEqual(i,j)
        blocks22 = load_indices("test/output/blocks.2.txt", groups=True)  # should not matter
        self.assertEqual(len(blocks),len(blocks2))
        for bl,bl1 in zip(blocks,blocks2):
            for i,j in zip(bl,bl1):
                self.assertEqual(i,j)
        blocks3 = load_indices("test/output/blocks.3.txt", groups=True)
        self.assertEqual(len(blocks),len(blocks3))
        for bl,bl1 in zip(blocks,blocks3):
            for i,j in zip(bl,bl1):
                self.assertEqual(i,j)

    def test_load_dump_indices2(self):
        randint = numpy.random.randint
        for counter in xrange(20):
            for compact in True, False:
                shift = randint(-5,5)
                indices = []
                for i in xrange(randint(10,20)):
                    l = list(set(randint(0,10,randint(20))))
                    if len(l) > 0:
                        indices.append(l)
                indices_flat = sum(indices, [])

                dump_indices("test/output/indices_blocks.txt", indices, compact=compact, shift=shift)
                dump_indices("test/output/indices_flat.txt", indices_flat, compact=compact, shift=shift)

                check = load_indices("test/output/indices_blocks.txt", shift=-shift, groups=True)
                self.assertEqual(indices, check)
                check = load_indices("test/output/indices_blocks.txt", shift=-shift)
                self.assertEqual(indices_flat, check)
                check = load_indices("test/output/indices_flat.txt", shift=-shift, groups=True)
                self.assertEqual([indices_flat], check)
                check = load_indices("test/output/indices_flat.txt", shift=-shift)
                self.assertEqual(indices_flat, check)

    def test_punch(self):
        mol0 = load_molecule_g03fchk('test/input/punch/gaussian.fchk')
        mol1 = load_molecule_g03fchk('test/input/punch/gaussian.fchk', fn_punch='test/input/punch/fort.7')
        assert abs(mol0.gradient - mol1.gradient).max() < 1e-8
        assert abs(mol0.hessian - mol1.hessian).max() < 1e-8

    def test_dftd3(self):
        assert load_dftd3('test/input/dftd3/dftd3.out') == -0.00330057
