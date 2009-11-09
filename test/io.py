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

from molmod.data.periodic import periodic
from molmod.units import angstrom, cm, amu, calorie, avogadro, eV
from molmod.constants import lightspeed

import unittest, numpy


__all__ = ["IOTestCase"]


class IOTestCase(unittest.TestCase):
    def test_load_fixed_g03com(self):
        fixed_atoms = load_fixed_g03com("input/mat/Zp_p_prod.18aug.com")
        self.assertEqual(len(fixed_atoms), 48)
        self.assertEqual(fixed_atoms, range(114,114+48))

    def test_load_molecule_g03fchk(self):
        atoms = 181
        molecule = load_molecule_g03fchk("input/mat/Zp_p_react.28aug.fchk","input/mat/Zp_p_react.14mei.fchk")

        self.assertEqual(molecule.hessian.shape,(atoms*3,atoms*3))
        self.assertAlmostEqual(molecule.energy, -18613.135744186180, 7)

    def test_load_molecule_g03fchkvdw(self):
        atoms = 179
        molecule = load_molecule_g03fchk("input/matvdw/R.fchk","input/matvdw/R_SCF.fchk","input/matvdw/R_b3lyp-d.out")

        self.assertEqual(molecule.hessian.shape,(atoms*3,atoms*3))
        self.assertAlmostEqual(molecule.energy,-18612.352569964281 , 7)

    def test_load_molecule_cp2k(self):
        molecule = load_molecule_cp2k("input/cp2k/pentane/opt.xyz", "input/cp2k/pentane/sp.out", "input/cp2k/pentane/freq.out")
        self.assertAlmostEqual(molecule.energy, 0.012255059530862)
        self.assertEqual(molecule.multiplicity, 1)
        self.assertEqual(molecule.numbers[0], 6)
        self.assertEqual(molecule.numbers[4], 1)
        self.assertAlmostEqual(molecule.masses[0], periodic[6].mass)
        self.assertAlmostEqual(molecule.masses[4], periodic[1].mass)
        self.assertAlmostEqual(molecule.coordinates[5,1]/angstrom, 13.9457396458)
        self.assertAlmostEqual(molecule.gradient[0,2], 0.0000000038, 9)
        self.assertAlmostEqual(molecule.gradient[11,0], 0.0000000177, 9)
        self.assertAlmostEqual(molecule.hessian[0,0], 1.08660340, 6)
        self.assertAlmostEqual(molecule.hessian[-1,-1], 0.528947590, 6)

    def test_load_molecule_cpmd(self):
        molecule = load_molecule_cpmd("input/cpmd/damp.out", "input/cpmd/GEOMETRY.xyz", "input/cpmd/MOLVIB")
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
        molecule = load_molecule_charmm("input/an/ethanol.cor","input/an/ethanol.hess.full")
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
        molecule = load_molecule_qchem("input/an/methanol.b3lypd.cc-pvtz.freq.out")
        self.assertAlmostEqual(molecule.energy, -115.7740000429)
        self.assertEqual(molecule.multiplicity, 1)
        self.assertEqual(molecule.numbers[0], 6)
        self.assertEqual(molecule.numbers[3], 1)
        self.assertAlmostEqual(molecule.masses[0]/amu, 12.01100, 2)
        self.assertAlmostEqual(molecule.masses[3]/amu, 1.00800, 2)
        self.assertAlmostEqual(molecule.coordinates[5,1]/angstrom, -0.756304)
        self.assertAlmostEqual(molecule.gradient[0,2]/(1000*calorie/avogadro/angstrom), -0.00, 5)
        self.assertAlmostEqual(molecule.gradient[5,0]/(1000*calorie/avogadro/angstrom), -0.00, 5)
        self.assertAlmostEqual(molecule.hessian[0,0], 0.4591905, 6)
        self.assertAlmostEqual(molecule.hessian[-1,-1], 0.0028270, 6)

    def test_load_molecule_vasp(self):
        molecule = load_molecule_vasp("input/vasp/xyz-structure","input/vasp/OUTCAR")
        self.assertEqual(molecule.numbers[0],14)
        self.assertEqual(molecule.numbers[107],13)
        self.assertAlmostEqual(molecule.masses[0]/amu, 28.085)
        self.assertAlmostEqual(molecule.masses[120]/amu, 1.000)
        self.assertAlmostEqual(molecule.coordinates[5,1]/angstrom, 2.93027 )
        self.assertAlmostEqual(molecule.gradient[0,2]/(eV/angstrom), -0.016534, 5)
        self.assertAlmostEqual(molecule.gradient[8,0]/(eV/angstrom), 0.035937, 5)
        self.assertAlmostEqual( - molecule.hessian[0,0]/(eV/angstrom**2), -46.644731, 6)
        self.assertAlmostEqual( - molecule.hessian[-1,-1]/(eV/angstrom**2), -5.524062, 6)
        # if VASP contains only a partial Hessian
        molecule = load_molecule_vasp("input/vasp/xyz-structure-part","input/vasp/OUTCAR-part")
        self.assertEqual(molecule.numbers[0],14)
        self.assertEqual(molecule.numbers[107],13)
        self.assertAlmostEqual(molecule.masses[0]/amu, 28.085)
        self.assertAlmostEqual(molecule.masses[120]/amu, 1.000)
        self.assertAlmostEqual(molecule.coordinates[5,1]/angstrom, 2.93027 )
        self.assertAlmostEqual(molecule.gradient[0,2]/(eV/angstrom), -0.016923, 5)
        self.assertAlmostEqual(molecule.gradient[8,0]/(eV/angstrom), 0.036518, 5)
        self.assertAlmostEqual( - molecule.hessian[0,0]/(eV/angstrom**2), -46.646216, 6)
        self.assertAlmostEqual( - molecule.hessian[-1,-1]/(eV/angstrom**2), -5.524077, 6)
        fixed = load_fixed_vasp("input/vasp/OUTCAR-part")
        self.assertEqual(fixed[0],2)
        self.assertEqual(fixed[30],53)

    def test_checkpoint(self):
        molecule = load_molecule_cp2k("input/cp2k/pentane/opt.xyz", "input/cp2k/pentane/sp.out", "input/cp2k/pentane/freq.out")
        nma1 = NMA(molecule)
        nma1.write_to_file("output/test.chk")
        nma2 = NMA.read_from_file("output/test.chk")

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


    def test_load_blocks_txt(self):
        blocks = load_blocks_txt("input/an/fixed.07.txt")

    def test_write_modes_for_VMD(self):
        molecule = load_molecule_charmm("input/an/ethanol.cor","input/an/ethanol.hess.full")
        nma = NMA(molecule)
        write_modes_for_VMD(molecule, nma, 6, filename="output/mode6.xyz")
        f = file("output/mode6.xyz")
        line = f.readline()    # first line
        words = line.split()
        self.assertEqual(len(words),1)
        self.assertEqual(words[0],"9")
        line = f.readline()    # 2nd line
        words = line.split()
        self.assertEqual(len(words),1)
        self.assertEqual(words[0],"i=0")
        line = f.readline()    # 3rd line
        words = line.split()
        self.assertEqual(len(words),4)
        self.assertEqual(words[0],"6")
        self.assertEqual(words[2],"0.0816")
        for i in range(9):
            line = f.readline()
        words = line.split()
        self.assertEqual(len(words),1)
        self.assertEqual(words[0],"9")
        f.close

