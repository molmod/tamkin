# -*- coding: utf-8 -*-
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
# "TAMkin: A Versatile Package for Vibrational Analysis and Chemical Kinetics",
# An Ghysels, Toon Verstraelen, Karen Hemelsoet, Michel Waroquier and Veronique
# Van Speybroeck, Journal of Chemical Information and Modeling, Articles ASAP
# (As Soon As Publishable)
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
# --


from tamkin import *

from molmod.periodic import periodic
from molmod import angstrom, amu, calorie, avogadro, electronvolt, lightspeed, \
    UnitCell

import unittest, numpy


__all__ = ["DataTestCase"]


class DataTestCase(unittest.TestCase):

    def test_get_submolecule(self):
        molecule = load_molecule_charmm("input/an/ethanol.cor","input/an/ethanol.hess.full")
        select = [3,4,2]
        molecule2 = molecule.get_submolecule(select)
        for i,at in enumerate(select):
            self.assertAlmostEqual(molecule.numbers[at], molecule2.numbers[i])
            self.assertAlmostEqual(molecule.masses[at], molecule2.masses[i])
            self.assertEqual(molecule.symbols[at], molecule2.symbols[i])
            for mu in range(3):
                self.assertAlmostEqual(molecule.coordinates[at,mu],molecule2.coordinates[i,mu])
                self.assertAlmostEqual(molecule.gradient[at,mu],molecule2.gradient[i,mu])
                for j,at2 in enumerate(select):
                    for nu in range(3):
                        self.assertAlmostEqual(molecule.hessian[3*at+mu,3*at2+nu],molecule2.hessian[3*i+mu,3*j+nu])
        self.assertEqual(molecule.multiplicity, molecule2.multiplicity)
        self.assertEqual(molecule.symmetry_number, molecule2.symmetry_number)
        self.assertEqual(molecule.periodic, molecule2.periodic)
        self.assertEqual(molecule.energy, molecule2.energy)

    def test_get_submolecule_cp2k(self):
        molecule = load_molecule_cp2k("input/cp2k/pentane/opt.xyz", "input/cp2k/pentane/sp.out", "input/cp2k/pentane/freq.out")
        select = range(5)+[9,11,14]
        molecule2 = molecule.get_submolecule(select, title="this is submol", energy=5., periodic=False, symmetry_number=6)  # just trying out something
        for i,at in enumerate(select):
            self.assertAlmostEqual(molecule.numbers[at], molecule2.numbers[i])
            self.assertAlmostEqual(molecule.masses[at], molecule2.masses[i])
            for mu in range(3):
                self.assertAlmostEqual(molecule.coordinates[at,mu],molecule2.coordinates[i,mu])
                self.assertAlmostEqual(molecule.gradient[at,mu],molecule2.gradient[i,mu])
                for j,at2 in enumerate(select):
                    for nu in range(3):
                        self.assertAlmostEqual(molecule.hessian[3*at+mu,3*at2+nu],molecule2.hessian[3*i+mu,3*j+nu])
        self.assertEqual(molecule.multiplicity, molecule2.multiplicity)
        self.assertEqual(molecule2.symmetry_number, 6)
        self.assertEqual(molecule2.periodic, False)
        self.assertEqual(molecule2.energy, 5.)

    def test_translate_pbc(self):
        molecule = load_molecule_cp2k("input/cp2k/pentane/opt.xyz", "input/cp2k/pentane/sp.out", "input/cp2k/pentane/freq.out")
        self.assertAlmostEqual(molecule.unit_cell.matrix[1,1]/angstrom, 30.000,3)
        self.assertAlmostEqual(molecule.unit_cell.matrix[0,1]/angstrom, 0.000,3)
        self.assertAlmostEqual(molecule.coordinates[5,1]/angstrom, 13.9457396458)
        selected = range(6)+[11,14]
        molecule2 = translate_pbc(molecule, selected, [1,-1,0])
        self.assertAlmostEqual(molecule2.coordinates[5,1]/angstrom, 13.9457396458-30.)
        for i in range(molecule.size):
            self.assertAlmostEqual(molecule.numbers[i], molecule2.numbers[i])
            self.assertAlmostEqual(molecule.masses[i], molecule2.masses[i])
            for mu in range(3):
                if i not in selected:
                    self.assertAlmostEqual(molecule.coordinates[i,mu],molecule2.coordinates[i,mu])
                else:
                    if mu is 0:
                        self.assertAlmostEqual(molecule.coordinates[i,mu]+30.*angstrom,molecule2.coordinates[i,mu])
                    if mu is 1:
                        self.assertAlmostEqual(molecule.coordinates[i,mu]-30.*angstrom,molecule2.coordinates[i,mu])
                    if mu is 2:
                        self.assertAlmostEqual(molecule.coordinates[i,mu],molecule2.coordinates[i,mu])
                self.assertAlmostEqual(molecule.gradient[i,mu],molecule2.gradient[i,mu])
                for j in selected:
                    for nu in range(3):
                        self.assertAlmostEqual(molecule.hessian[3*i+mu,3*j+nu],molecule2.hessian[3*i+mu,3*j+nu])
        self.assertEqual(molecule.multiplicity, molecule2.multiplicity)
        self.assertEqual(molecule.symmetry_number, molecule2.symmetry_number)
        self.assertAlmostEqual(molecule2.unit_cell.matrix[0,0]/angstrom, 30.000,3)
        self.assertAlmostEqual(molecule2.unit_cell.matrix[1,2]/angstrom, 0.000,3)

    def test_molecule_checkpoint_basic(self):
        mol1 = load_molecule_g03fchk("input/sterck/aa.fchk")
        mol1.write_to_file("output/molecule_checkpoint_basic.chk")
        mol2 = Molecule.read_from_file("output/molecule_checkpoint_basic.chk")

        self.assertEqual(mol1.numbers.shape, mol2.numbers.shape)
        self.assert_(abs(mol1.numbers - mol2.numbers).max() == 0)
        self.assertEqual(mol1.coordinates.shape, mol2.coordinates.shape)
        self.assert_(abs(mol1.coordinates - mol2.coordinates).max() < 1e-10)
        self.assertEqual(mol1.masses.shape, mol2.masses.shape)
        self.assert_(abs(mol1.masses - mol2.masses).max() < 1e-10)
        self.assertAlmostEqual(mol1.energy, mol2.energy)
        self.assertEqual(mol1.gradient.shape, mol2.gradient.shape)
        self.assert_(abs(mol1.gradient - mol2.gradient).max() < 1e-10)
        self.assertEqual(mol1.hessian.shape, mol2.hessian.shape)
        self.assert_(abs(mol1.hessian - mol2.hessian).max() < 1e-10)
        self.assertEqual(mol1.multiplicity, mol2.multiplicity)
        self.assertEqual(mol1.symmetry_number, mol2.symmetry_number)
        self.assertEqual(mol1.periodic, mol2.periodic)

    def test_molecule_checkpoint_full(self):
        mol1 = load_molecule_g03fchk("input/sterck/aa.fchk")
        mol1.set_default_graph()
        mol1.unit_cell = UnitCell(numpy.identity(3, float)*25)
        mol1.symbols = [periodic[n].symbol for n in mol1.numbers]
        mol1.write_to_file("output/molecule_checkpoint_full.chk")
        mol2 = Molecule.read_from_file("output/molecule_checkpoint_full.chk")

        self.assertEqual(mol1.numbers.shape, mol2.numbers.shape)
        self.assert_(abs(mol1.numbers - mol2.numbers).max() == 0)
        self.assertEqual(mol1.coordinates.shape, mol2.coordinates.shape)
        self.assert_(abs(mol1.coordinates - mol2.coordinates).max() < 1e-10)
        self.assertEqual(mol1.masses.shape, mol2.masses.shape)
        self.assert_(abs(mol1.masses - mol2.masses).max() < 1e-10)
        self.assertAlmostEqual(mol1.energy, mol2.energy)
        self.assertEqual(mol1.gradient.shape, mol2.gradient.shape)
        self.assert_(abs(mol1.gradient - mol2.gradient).max() < 1e-10)
        self.assertEqual(mol1.hessian.shape, mol2.hessian.shape)
        self.assert_(abs(mol1.hessian - mol2.hessian).max() < 1e-10)
        self.assertEqual(mol1.multiplicity, mol2.multiplicity)
        self.assertEqual(mol1.symmetry_number, mol2.symmetry_number)
        self.assertEqual(mol1.periodic, mol2.periodic)

        self.assertEqual(mol1.graph.edges, mol2.graph.edges)
        self.assertEqual(mol1.title, mol2.title)
        self.assertEqual(mol1.unit_cell.matrix.shape, mol2.unit_cell.matrix.shape)
        self.assert_(abs(mol1.unit_cell.matrix - mol2.unit_cell.matrix).max() < 1e-10)
        self.assertEqual(mol1.unit_cell.active.shape, mol2.unit_cell.active.shape)
        self.assert_((mol1.unit_cell.active == mol2.unit_cell.active).all())
        self.assert_(all((s1==s2 for s1, s2 in zip(mol1.symbols, mol2.symbols))))

    def test_copy_with(self):
        mol1 = load_molecule_g03fchk("input/sterck/aa.fchk")
        mol2 = mol1.copy_with(title="foo")
        self.assertEqual(mol2.title, "foo")

    def test_get_external_basis_new1(self):
        mol = load_molecule_g03fchk("input/linear/gaussian.fchk")
        eb = mol.get_external_basis_new()
        self.assertEqual(eb.shape[0], 5)
        # orthogonality test: the external basis should be orthogonal
        # to the basis of the four internal coordinates: stretch1, stertch2,
        # bend1 and bend2
        ib = numpy.array([
            [0,0,1,0,0,-1,0,0,0],
            [0,0,1,0,0,0,0,0,-1],
            [0,-1,0,0,0.5,0,0,0.5,0],
            [-1,0,0,0.5,0,0,0.5,0,0],
        ])
        error = abs(numpy.dot(ib, eb.transpose())).max()
        self.assert_(error < 1e-5)

    def test_get_external_basis_new2(self):
        mol = load_molecule_g03fchk("input/ethane/gaussian.fchk")
        eb = mol.get_external_basis_new()
        self.assertEqual(eb.shape[0], 6)
        # orthogonality test: the external basis should be orthogonal
        # to the basis of the internal coordinates. in this case we construct
        # a redundant basis of internal coordinates based on all interatomic
        # distances. this is overkill...
        ib = []
        for i in xrange(mol.size):
            for j in xrange(i):
                b = numpy.zeros((mol.size, 3), float)
                delta = mol.coordinates[j] - mol.coordinates[i]
                delta /= numpy.linalg.norm(delta)
                b[i] = delta
                b[j] = -delta
                ib.append(b.ravel())
        ib = numpy.array(ib)
        error = abs(numpy.dot(ib, eb.transpose())).max()
        self.assert_(error < 1e-5)
