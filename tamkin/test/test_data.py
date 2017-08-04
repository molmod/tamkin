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

import numpy as np
import pkg_resources

from molmod.periodic import periodic
from molmod import angstrom, amu, calorie, avogadro, electronvolt, lightspeed, \
    UnitCell
from molmod.test.common import tmpdir

from tamkin import *


def test_get_submolecule():
    molecule = load_molecule_charmm(
        pkg_resources.resource_filename(__name__, "../data/test/an/ethanol.cor"),
        pkg_resources.resource_filename(__name__, "../data/test/an/ethanol.hess.full"))
    select = [3,4,2]
    molecule2 = molecule.get_submolecule(select)
    for i,at in enumerate(select):
        assert molecule.numbers[at] == molecule2.numbers[i]
        assert abs(molecule.masses[at] - molecule2.masses[i]) < 1e-3
        assert molecule.symbols[at] == molecule2.symbols[i]
        for mu in range(3):
            assert abs(molecule.coordinates[at,mu] - molecule2.coordinates[i,mu]) < 1e-5
            assert abs(molecule.gradient[at,mu] - molecule2.gradient[i,mu]) < 1e-5
            for j,at2 in enumerate(select):
                for nu in range(3):
                    assert abs(molecule.hessian[3*at+mu,3*at2+nu] - molecule2.hessian[3*i+mu,3*j+nu]) < 1e-5
    assert molecule.multiplicity == molecule2.multiplicity
    assert molecule.symmetry_number == molecule2.symmetry_number
    assert molecule.periodic == molecule2.periodic
    assert molecule.energy == molecule2.energy


def test_get_submolecule_cp2k():
    molecule = load_molecule_cp2k(
        pkg_resources.resource_filename(__name__, "../data/test/cp2k/pentane/sp.out"),
        pkg_resources.resource_filename(__name__, "../data/test/cp2k/pentane/freq.out"))
    select = range(5)+[9,11,14]
    molecule2 = molecule.get_submolecule(select, title="this is submol", energy=5., periodic=False, symmetry_number=6)  # just trying out something
    for i,at in enumerate(select):
        assert molecule.numbers[at] == molecule2.numbers[i]
        assert abs(molecule.masses[at] - molecule2.masses[i]) < 1e-3
        for mu in range(3):
            assert abs(molecule.coordinates[at,mu] - molecule2.coordinates[i,mu]) < 1e-5
            assert abs(molecule.gradient[at,mu] - molecule2.gradient[i,mu]) < 1e-5
            for j,at2 in enumerate(select):
                for nu in range(3):
                    assert abs(molecule.hessian[3*at+mu,3*at2+nu] - molecule2.hessian[3*i+mu,3*j+nu]) < 1e-5
    assert molecule.multiplicity == molecule2.multiplicity
    assert molecule2.symmetry_number == 6
    assert molecule2.periodic == False
    assert molecule2.energy == 5.


def test_translate_pbc():
    molecule = load_molecule_cp2k(
        pkg_resources.resource_filename(__name__, "../data/test/cp2k/pentane/sp.out"),
        pkg_resources.resource_filename(__name__, "../data/test/cp2k/pentane/freq.out"))
    assert abs(molecule.unit_cell.matrix[1,1]/angstrom - 30.000) < 1e-3
    assert abs(molecule.unit_cell.matrix[0,1]/angstrom - 0.000) < 1e-3
    assert abs(molecule.coordinates[5,1]/angstrom - 13.928520) < 1e-5
    selected = range(6)+[11,14]
    molecule2 = translate_pbc(molecule, selected, [1,-1,0])
    assert abs(molecule2.coordinates[5,1]/angstrom -  (13.928520-30.0)) < 1e-5
    for i in range(molecule.size):
        assert molecule.numbers[i] == molecule2.numbers[i]
        assert abs(molecule.masses[i] - molecule2.masses[i]) < 1e-3
        for mu in range(3):
            if i not in selected:
                assert abs(molecule.coordinates[i,mu] - molecule2.coordinates[i,mu]) < 1e-3
            else:
                if mu is 0:
                    assert abs(molecule.coordinates[i,mu] + 30.*angstrom - molecule2.coordinates[i,mu]) < 1e-3
                if mu is 1:
                    assert abs(molecule.coordinates[i,mu] - 30.*angstrom - molecule2.coordinates[i,mu]) < 1e-3
                if mu is 2:
                    assert abs(molecule.coordinates[i,mu] - molecule2.coordinates[i,mu]) < 1e-3
            assert abs(molecule.gradient[i,mu] - molecule2.gradient[i,mu]) < 1e-3
            for j in selected:
                for nu in range(3):
                    assert abs(molecule.hessian[3*i+mu,3*j+nu] - molecule2.hessian[3*i+mu,3*j+nu]) < 1e-3
    assert molecule.multiplicity == molecule2.multiplicity
    assert molecule.symmetry_number == molecule2.symmetry_number
    assert abs(molecule2.unit_cell.matrix[0,0]/angstrom - 30.000) < 1e-3
    assert abs(molecule2.unit_cell.matrix[1,2]/angstrom - 0.000) < 1e-3


def test_molecule_checkpoint_basic():
    mol1 = load_molecule_g03fchk(
        pkg_resources.resource_filename(__name__, "../data/test/sterck/aa.fchk"))
    with tmpdir(__name__, 'test_molecule_checkpoint_basic') as dn:
        fn_out = os.path.join(dn, "molecule_checkpoint_basic.chk")
        mol1.write_to_file(fn_out)
        mol2 = Molecule.read_from_file(fn_out)

    assert mol1.numbers.shape == mol2.numbers.shape
    assert abs(mol1.numbers - mol2.numbers).max() == 0
    assert mol1.coordinates.shape == mol2.coordinates.shape
    assert abs(mol1.coordinates - mol2.coordinates).max() < 1e-10
    assert mol1.masses.shape == mol2.masses.shape
    assert abs(mol1.masses - mol2.masses).max() < 1e-10
    assert abs(mol1.energy - mol2.energy) < 1e-5
    assert mol1.gradient.shape == mol2.gradient.shape
    assert abs(mol1.gradient - mol2.gradient).max() < 1e-10
    assert mol1.hessian.shape == mol2.hessian.shape
    assert abs(mol1.hessian - mol2.hessian).max() < 1e-10
    assert mol1.multiplicity == mol2.multiplicity
    assert mol1.symmetry_number == mol2.symmetry_number
    assert mol1.periodic == mol2.periodic


def test_molecule_checkpoint_full():
    mol1 = load_molecule_g03fchk(
        pkg_resources.resource_filename(__name__, "../data/test/sterck/aa.fchk"))
    mol1.set_default_graph()
    mol1.unit_cell = UnitCell(np.identity(3, float)*25)
    mol1.symbols = [periodic[n].symbol for n in mol1.numbers]
    with tmpdir(__name__, 'test_molecule_checkpoint_full') as dn:
        fn_out = os.path.join(dn, "molecule_checkpoint_full.chk")
        mol1.write_to_file(fn_out)
        mol2 = Molecule.read_from_file(fn_out)

    assert mol1.numbers.shape == mol2.numbers.shape
    assert abs(mol1.numbers - mol2.numbers).max() == 0
    assert mol1.coordinates.shape == mol2.coordinates.shape
    assert abs(mol1.coordinates - mol2.coordinates).max() < 1e-10
    assert mol1.masses.shape == mol2.masses.shape
    assert abs(mol1.masses - mol2.masses).max() < 1e-10
    assert abs(mol1.energy - mol2.energy) < 1e-10
    assert mol1.gradient.shape == mol2.gradient.shape
    assert abs(mol1.gradient - mol2.gradient).max() < 1e-10
    assert mol1.hessian.shape == mol2.hessian.shape
    assert abs(mol1.hessian - mol2.hessian).max() < 1e-10
    assert mol1.multiplicity == mol2.multiplicity
    assert mol1.symmetry_number == mol2.symmetry_number
    assert mol1.periodic == mol2.periodic

    assert mol1.graph.edges == mol2.graph.edges
    assert mol1.title == mol2.title
    assert mol1.unit_cell.matrix.shape == mol2.unit_cell.matrix.shape
    assert abs(mol1.unit_cell.matrix - mol2.unit_cell.matrix).max() < 1e-10
    assert mol1.unit_cell.active.shape == mol2.unit_cell.active.shape
    assert (mol1.unit_cell.active == mol2.unit_cell.active).all()
    assert mol1.symbols == mol2.symbols


def test_copy_with():
    mol1 = load_molecule_g03fchk(
        pkg_resources.resource_filename(__name__, "../data/test/sterck/aa.fchk"))
    mol2 = mol1.copy_with(title="foo")
    assert mol2.title == "foo"


def test_get_external_basis_new1():
    mol = load_molecule_g03fchk(
        pkg_resources.resource_filename(__name__, "../data/test/linear/gaussian.fchk"))
    eb = mol.get_external_basis_new()
    assert eb.shape[0] == 5
    # orthogonality test: the external basis should be orthogonal
    # to the basis of the four internal coordinates: stretch1, stertch2,
    # bend1 and bend2
    ib = np.array([
        [0,0,1,0,0,-1,0,0,0],
        [0,0,1,0,0,0,0,0,-1],
        [0,-1,0,0,0.5,0,0,0.5,0],
        [-1,0,0,0.5,0,0,0.5,0,0],
    ])
    error = abs(np.dot(ib, eb.transpose())).max()
    assert error < 1e-5


def test_get_external_basis_new2():
    mol = load_molecule_g03fchk(
        pkg_resources.resource_filename(__name__, "../data/test/ethane/gaussian.fchk"))
    eb = mol.get_external_basis_new()
    assert eb.shape[0] == 6
    # orthogonality test: the external basis should be orthogonal
    # to the basis of the internal coordinates. in this case we construct
    # a redundant basis of internal coordinates based on all interatomic
    # distances. this is overkill...
    ib = []
    for i in xrange(mol.size):
        for j in xrange(i):
            b = np.zeros((mol.size, 3), float)
            delta = mol.coordinates[j] - mol.coordinates[i]
            delta /= np.linalg.norm(delta)
            b[i] = delta
            b[j] = -delta
            ib.append(b.ravel())
    ib = np.array(ib)
    error = abs(np.dot(ib, eb.transpose())).max()
    assert error < 1e-5


def test_rot_scan_ts():
    # The select dihedral angles do not allow an automatic assignment of the top
    # indexes.
    mol = load_molecule_g03fchk(
        pkg_resources.resource_filename(__name__, "../data/test/sterck/paats_1h2o_b_aa.fchk"))
    diheds = [[9, 6, 7, 10], [7, 10, 18, 11], [11, 18, 10, 7]]
    for dihed in diheds:
        try:
            scan = RotScan(np.array(dihed)-1, mol)
            assert False
        except ValueError, e:
            assert e.message == "The rotating top could not be assigned properly. Specify the top_indexes manually."
