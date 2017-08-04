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

from molmod import lightspeed, boltzmann, second, atm, amu, meter, mol, \
    kcalmol, calorie, kelvin, kjmol, planck
from molmod.test.common import tmpdir

from tamkin import *


__all__ = ["PartFunTestCase"]


class PartFunTestCase(unittest.TestCase):
    def test_phva_react_mat(self):
        fixed_atoms = load_fixed_g03com(
            pkg_resources.resource_filename(__name__, "../data/test/mat/Zp_p_react.14mei.com"))
        molecule = load_molecule_g03fchk(
            pkg_resources.resource_filename(__name__, "../data/test/mat/Zp_p_react.28aug.fchk"),
            pkg_resources.resource_filename(__name__, "../data/test/mat/Zp_p_react.14mei.fchk"))
        nma = NMA(molecule, PHVA(fixed_atoms))
        pf = PartFun(nma)

        # values obtained with frektsjek.exe
        self.assertAlmostEqual(pf.log(650)-(-pf.electronic.energy/(boltzmann*650)), np.log(1.69428264E-45),3)
        self.assertAlmostEqual(pf.log(700)-(-pf.electronic.energy/(boltzmann*700)), np.log(3.26776842E-29),3)
        self.assertAlmostEqual(pf.log(750)-(-pf.electronic.energy/(boltzmann*750)), np.log(2.00896366E-14),3)

    def test_phva_rate_const_mat(self):
        fixed_atoms = load_fixed_g03com(
            pkg_resources.resource_filename(__name__, "../data/test/mat/Zp_p_react.14mei.com"))
        mol_react = load_molecule_g03fchk(
            pkg_resources.resource_filename(__name__, "../data/test/mat/Zp_p_react.28aug.fchk"),
            pkg_resources.resource_filename(__name__, "../data/test/mat/Zp_p_react.14mei.fchk"))
        mol_trans = load_molecule_g03fchk(
            pkg_resources.resource_filename(__name__, "../data/test/mat/Zp_p_TS.28aug.fchk"),
            pkg_resources.resource_filename(__name__, "../data/test/mat/5Tp_p_TS.oniom21apr_HF.fchk"))
        pf_react = PartFun(NMA(mol_react, PHVA(fixed_atoms)))
        pf_trans = PartFun(NMA(mol_trans, PHVA(fixed_atoms)))
        km = KineticModel([pf_react], pf_trans)

        # values taken from the fancy excel file...
        temps = np.array([670,680,690,700,710,720,730,740,750,760,770])
        expected_ks = np.array([
            7.9473102E+05, 9.8300444E+05, 1.2085262E+06, 1.4771808E+06,
            1.7955340E+06, 2.1708793E+06, 2.6112829E+06, 3.1256298E+06,
            3.7236678E+06, 4.4160510E+06, 5.2143822E+06
        ])
        for i in xrange(len(temps)):
            k = km.rate_constant(temps[i])
            self.assertAlmostEqual(np.log(k/(1/second)), np.log(expected_ks[i]),5)
            log_k = km.rate_constant(temps[i], do_log=True)
            self.assertAlmostEqual(np.log(k), log_k)

    def test_gas_react_sterck(self):
        # Test both Full and ConstrainExt:
        for treatment, precision_wn, shift in (Full(), 0, 6), (ConstrainExt(), 3, 0):
            ## aa.fchk
            molecule = load_molecule_g03fchk(
                pkg_resources.resource_filename(__name__, "../data/test/sterck/aa.fchk"))
            nma = NMA(molecule, treatment)
            pf = PartFun(nma, [ExtTrans(), ExtRot(1)])
            temp = 298.150

            # check partition function, values from aa.log
            self.assertAlmostEqual(16.973928, pf.translational.log(temp), 6)
            self.assertEqual(pf.rotational.count, 3)
            self.assertAlmostEqual(11.225093, pf.rotational.log(temp), 6)
            vib_contribs = pf.vibrational.log_terms(temp)
            expected_vib_contribs = np.array([
                0.674278, 0.292167, -0.357617, -1.017249, -1.018740, -1.427556,
                -1.428865
            ])
            for i in xrange(len(expected_vib_contribs)):
                self.assertAlmostEqual(vib_contribs[i], expected_vib_contribs[i], 2)
            self.assertAlmostEqual(-53.068692, pf.log(temp)-(-pf.electronic.energy/(boltzmann*temp)), 2)

            ## aa.fchk, rotational symmetry number is computed by molmod
            pf = PartFun(nma, [ExtTrans(), ExtRot()])
            self.assertEqual(pf.rotational.symmetry_number, 1)

            # check partition function, values from aa.log
            self.assertAlmostEqual(16.973928, pf.translational.log(temp), 6)
            self.assertEqual(pf.rotational.count, 3)
            self.assertAlmostEqual(11.225093, pf.rotational.log(temp), 6)
            vib_contribs = pf.vibrational.log_terms(temp)
            expected_vib_contribs = np.array([
                0.674278, 0.292167, -0.357617, -1.017249, -1.018740, -1.427556,
                -1.428865
            ])
            for i in xrange(len(expected_vib_contribs)):
                self.assertAlmostEqual(vib_contribs[i], expected_vib_contribs[i], 2)
            self.assertAlmostEqual(-53.068692, pf.log(temp)-(-pf.electronic.energy/(boltzmann*temp)), 2)

            ## aarad.fchk
            molecule = load_molecule_g03fchk(
                pkg_resources.resource_filename(__name__, "../data/test/sterck/aarad.fchk"))
            nma = NMA(molecule, treatment)
            pf = PartFun(nma, [ExtTrans(), ExtRot(1)])

            # check partition function, values taken from aarad.log
            self.assertAlmostEqual(16.995059, pf.translational.log(temp), 6)
            self.assertEqual(pf.rotational.count, 3)
            self.assertAlmostEqual(11.319073, pf.rotational.log(temp), 6)
            vib_contribs = pf.vibrational.log_terms(temp)
            expected_vib_contribs = np.array([
                0.959168, 0.413328, -0.287477, -0.371573, -0.983851, -1.040910,
                -1.311721, -1.428436,
            ])
            for i in xrange(len(expected_vib_contribs)):
                self.assertAlmostEqual(vib_contribs[i], expected_vib_contribs[i], 2+precision_wn)
            self.assertAlmostEqual(-61.738525, pf.log(temp)-(-pf.electronic.energy/(boltzmann*temp)), 1+precision_wn)

    def test_gas_trans_sterck(self):
        # Test both Full and ConstrainExt:
        for treatment, precision_wn, shift in (Full(), 0, 6), (ConstrainExt(), 3, 0):
            ## paats.fchk
            molecule = load_molecule_g03fchk(
                pkg_resources.resource_filename(__name__, "../data/test/sterck/paats.fchk"))
            nma = NMA(molecule, treatment)
            pf = PartFun(nma, [ExtTrans(), ExtRot(1)])
            temp = 298.150

            # check the natural logarithm partition function (split contributions):
            # (values taken from paats.log)
            self.assertAlmostEqual(18.024251, pf.translational.log(temp), 5)
            self.assertEqual(pf.rotational.count, 3)
            self.assertAlmostEqual(13.615243, pf.rotational.log(temp), 6)
            vib_contribs = pf.vibrational.log_terms(temp)
            expected_vib_contribs = np.array([
                1.749142, 1.335079, 0.789626, 0.704324, 0.365201, 0.061403,
                -0.060710, -0.281250, -0.309673, -0.423131, -0.447578, -0.983086,
                -1.017063, -1.107892, -1.192276, -1.318341, -1.352385, -1.427939,
            ])
            for i in xrange(len(expected_vib_contribs)):
                self.assertAlmostEqual(vib_contribs[i], expected_vib_contribs[i], 2+precision_wn)
            self.assertAlmostEqual(-139.302816, pf.log(temp)-(-pf.electronic.energy/(boltzmann*temp)), 1+precision_wn)

    def test_gas_rate_const_sterck(self):
        mol_react1 = load_molecule_g03fchk(
            pkg_resources.resource_filename(__name__, "../data/test/sterck/aa.fchk"))
        mol_react2 = load_molecule_g03fchk(
            pkg_resources.resource_filename(__name__, "../data/test/sterck/aarad.fchk"))
        mol_trans = load_molecule_g03fchk(
            pkg_resources.resource_filename(__name__, "../data/test/sterck/paats.fchk"))
        pf_react1 = PartFun(NMA(mol_react1, ConstrainExt()), [ExtTrans(cp=False), ExtRot(1)])
        pf_react2 = PartFun(NMA(mol_react2, ConstrainExt()), [ExtTrans(cp=False), ExtRot(1)])
        pf_trans = PartFun(NMA(mol_trans, ConstrainExt()), [ExtTrans(cp=False), ExtRot(1)])
        km = KineticModel([pf_react1, pf_react2], pf_trans)

        # values taken from the fancy excel file...
        temps = np.array([298.15,300,400,500,600,700,800,900,1000,1100])
        expected_ks = np.array([
            6.44881E-03,
            6.87398E-03, 9.74722E-02, 5.36970E-01, 1.82272E+00, 4.65134E+00,
            9.86811E+00, 1.84247E+01, 3.13480E+01, 4.97189E+01,
        ])
        unit = meter**3/mol/second
        for i in xrange(len(temps)):
            k = km.rate_constant(temps[i])
            # Sometimes, the fancy excel files use slightly different constants.
            # Therefore, only expect numbers to be equal up to 2 decimals.
            self.assertAlmostEqual(np.log(k/unit), np.log(expected_ks[i]), 2)

    def test_derivatives(self):
        molecule = load_molecule_g03fchk(
            pkg_resources.resource_filename(__name__, "../data/test/ethane/gaussian.fchk"))
        nma = NMA(molecule)
        rotscan = load_rotscan_g03log(
            pkg_resources.resource_filename(__name__, "../data/test/rotor/gaussian.log"))
        rotor = Rotor(rotscan, molecule, rotsym=3, even=True)
        for classical in True, False:
            for cp in True, False:
                pf = PartFun(nma, [
                    ExtTrans(cp), ExtRot(), rotor,
                    Vibrations(classical, 0.5, 0.7),
                ])
                self.assertEqual(pf.rotational.symmetry_number, 6)

                eps = 0.0001
                temps = np.array([300.0,400.0,500.0,600.0,700.0])
                sfs = [
                    pf.electronic, pf.translational, pf.rotational, pf.vibrational,
                    pf.hindered_rotor_3_4_5, pf
                ]
                for stat_fys in sfs:
                    for temp in temps:
                        # check the first derivative towards temperature with finite differences
                        a = stat_fys.logt(temp)
                        b = (stat_fys.log(temp+eps) -
                             stat_fys.log(temp-eps))/(2*eps)
                        self.assertAlmostEqual(
                            a, b, 6,
                            "error in partial derivative (%s, classical=%s, "
                            "cp=%s): %s!=%s" % (stat_fys.name, classical, cp,
                            a, b)
                        )
                        # check the second derivative towards temperature with finite differences
                        a = stat_fys.logtt(temp)
                        b = (stat_fys.logt(temp+eps) -
                             stat_fys.logt(temp-eps))/(2*eps)
                        self.assertAlmostEqual(
                            a, b, 6,
                            "error in second derivative (%s, classical=%s, "
                            "cp=%s): %s!=%s" % (stat_fys.name, classical, cp,
                            a, b)
                        )
                        # check the helper functions temperature argument
                        self.assertAlmostEqual(stat_fys.helper(temp,1), stat_fys.helper(temp,0)*temp)
                        self.assertAlmostEqual(stat_fys.helpert(temp,1), stat_fys.helpert(temp,0)*temp)
                        self.assertAlmostEqual(stat_fys.helpertt(temp,1), stat_fys.helpertt(temp,0)*temp)

    def test_chemical_potential_cp(self):
        molecule = load_molecule_g03fchk(
            pkg_resources.resource_filename(__name__, "../data/test/ethane/gaussian.fchk"))
        nma = NMA(molecule)
        pf = PartFun(nma, [ExtTrans(cp=True), ExtRot()])
        N0 = 1.0
        p0 = pf.translational.pressure
        N1 = 1.001
        p1 = p0*N1/N0
        self.assertEqual(pf.rotational.symmetry_number, 6)
        temps = np.array([300.0,400.0,500.0,600.0,700.0])
        for temp in temps:
            # compare the derivative of the free energy at N=1 with its
            # finit difference approximation.
            pf.pressure = p0
            F0 = pf.free_energy(temp)*N0
            pf.pressure = p1
            F1 = pf.free_energy(temp)*N1
            pf.pressure = 0.5*(p0+p1)
            mu = pf.chemical_potential(temp)
            error = abs(F1 - F0 - (N1-N0)*mu)
            self.assert_(error < 1e-12)
            # check the helper function and the definition of the chemical
            # potential in an ideal gas
            self.assertAlmostEqual(pf.translational.helpern(temp, 0) - pf.translational.helper(temp, 0), 0.0)
            self.assertAlmostEqual(pf.electronic.helpern(temp, 0) - pf.electronic.helper(temp, 0), 0.0)
            self.assertAlmostEqual(pf.rotational.helpern(temp, 0) - pf.rotational.helper(temp, 0), 0.0)
            self.assertAlmostEqual(pf.vibrational.helpern(temp, 0) - pf.vibrational.helper(temp, 0), 0.0)
            self.assertAlmostEqual(pf.helpern(temp, 0) - pf.helper(temp, 0), 0.0)
            self.assertAlmostEqual(pf.chemical_potential(temp), pf.free_energy(temp))

    def test_chemical_potential_cv(self):
        molecule = load_molecule_g03fchk(
            pkg_resources.resource_filename(__name__, "../data/test/ethane/gaussian.fchk"))
        nma = NMA(molecule)
        pf = PartFun(nma, [ExtTrans(cp=False), ExtRot()])
        N0 = 1.0
        r0 = pf.translational.density
        N1 = 1.001
        r1 = r0*N1/N0
        self.assertEqual(pf.rotational.symmetry_number, 6)
        temps = np.array([300.0,400.0,500.0,600.0,700.0])
        for temp in temps:
            # compare the derivative of the free energy at N=1 with its
            # finit difference approximation.
            pf.translational.density = r0
            F0 = pf.free_energy(temp)*N0
            pf.translational.density = r1
            F1 = pf.free_energy(temp)*N1
            mu = pf.chemical_potential(temp)
            error = abs(F1 - F0 - (N1-N0)*mu)
            self.assert_(error < 1e-8)
            # check the helper function and the definition of the chemical
            # potential in an ideal gas
            self.assertAlmostEqual(pf.translational.helpern(temp, 0) - pf.translational.helper(temp, 0), -1)
            self.assertAlmostEqual(pf.electronic.helpern(temp, 0) - pf.electronic.helper(temp, 0), 0.0)
            self.assertAlmostEqual(pf.rotational.helpern(temp, 0) - pf.rotational.helper(temp, 0), 0.0)
            self.assertAlmostEqual(pf.vibrational.helpern(temp, 0) - pf.vibrational.helper(temp, 0), 0.0)
            self.assertAlmostEqual(pf.helpern(temp, 0) - pf.helper(temp, 0), -1)
            self.assertAlmostEqual(pf.chemical_potential(temp), pf.free_energy(temp) + boltzmann*temp)

    def test_logn(self):
        for cp in False, True:
            molecule = load_molecule_g03fchk(
                pkg_resources.resource_filename(__name__, "../data/test/ethane/gaussian.fchk"))
            nma = NMA(molecule)
            pf = PartFun(nma, [ExtTrans(cp=cp), ExtRot()])
            temps = np.array([300.0,400.0,500.0,600.0,700.0])
            for temp in temps:
                if cp:
                    density = pf.translational.pressure/(boltzmann*temp)
                else:
                    density = pf.translational.density
                self.assertAlmostEqual(pf.translational.logn(temp) - pf.translational.logv(temp), -np.log(density))

    def test_logv(self):
        for cp in False, True:
            molecule = load_molecule_g03fchk(
                pkg_resources.resource_filename(__name__, "../data/test/ethane/gaussian.fchk"))
            nma = NMA(molecule)
            pf = PartFun(nma, [ExtTrans(cp=cp), ExtRot()])
            temps = np.array([300.0,400.0,500.0,600.0,700.0])
            for temp in temps:
                self.assertAlmostEqual(pf.translational.logv(temp), pf.translational._z1(temp))
                self.assertAlmostEqual(pf.electronic.logv(temp) - pf.electronic.log(temp), 0.0)
                self.assertAlmostEqual(pf.rotational.logv(temp) - pf.rotational.log(temp), 0.0)
                self.assertAlmostEqual(pf.vibrational.logv(temp) - pf.vibrational.log(temp), 0.0)

    def test_derived_quantities(self):
        # internal energy, heat capacity and entropy
        molecule = load_molecule_g03fchk(
            pkg_resources.resource_filename(__name__, "../data/test/sterck/aa.fchk"))
        nma = NMA(molecule, ConstrainExt())
        pf = PartFun(nma, [ExtTrans(cp=True), ExtRot(1)])

        # values taken from aa.log
        calmolK = calorie/mol/kelvin
        R = 1.98720649773 # R in calorie/(mol*K)
        temp = 298.15
        RT = R*temp*0.001 # RT in kcal/(mol*K)
        # electronic
        self.assertAlmostEqual(pf.electronic.internal_heat(temp), pf.electronic.energy)
        self.assertAlmostEqual(pf.electronic.heat_capacity(temp)/(calmolK), 0.0)
        self.assertAlmostEqual(pf.electronic.entropy(temp)/(calmolK), 0.0)
        # translational
        # WARNING: internal_heat returns enthalpy, turn it into internal energy
        self.assertAlmostEqual(pf.translational.internal_heat(temp)/(kcalmol)-RT, 0.889, 2)
        # WARNING: heat_capacity returns heat capacity at constant pressure,
        #          turn it into heat capacity at constant volume
        self.assertAlmostEqual(pf.translational.heat_capacity(temp)/(calmolK)-R, 2.981, 2)
        self.assertAlmostEqual(pf.translational.entropy(temp)/(calmolK), 38.699, 2)
        # rotational
        self.assertAlmostEqual(pf.rotational.internal_heat(temp)/(kcalmol), 0.889, 2)
        self.assertAlmostEqual(pf.rotational.heat_capacity(temp)/(calmolK), 2.981, 2)
        self.assertAlmostEqual(pf.rotational.entropy(temp)/(calmolK), 25.287, 2)
        # vibrational
        self.assertAlmostEqual(pf.vibrational.internal_heat(temp)/(kcalmol), 51.343, 2)
        self.assertAlmostEqual(pf.vibrational.heat_capacity(temp)/(calmolK), 13.264, 2)
        self.assertAlmostEqual(pf.vibrational.entropy(temp)/(calmolK), 10.710, 2)
        # total
        # WARNING: internal_heat returns enthalpy, turn it into internal energy
        # WARNING: internal_heat returns enthalpy that includes the electronic energy
        self.assertAlmostEqual((pf.internal_heat(temp)-pf.electronic.energy)/(kcalmol)-RT, 53.121, 2)
        # WARNING: heat_capacity returns heat capacity at constant pressure,
        #          turn it into heat capacity at constant volume
        self.assertAlmostEqual(pf.heat_capacity(temp)/(calmolK)-R, 19.225, 2)
        self.assertAlmostEqual(pf.entropy(temp)/(calmolK), 74.696, 2)
        self.assertAlmostEqual(pf.free_energy(0.0), -247.228749, 2) # Zero-point energy
        self.assertAlmostEqual(pf.internal_heat(temp), -247.222989, 2)
        # switch to constant pressure for the enthalpy and Gibbs free energy.
        self.assertAlmostEqual(pf.internal_heat(temp), -247.222045, 2)
        self.assertAlmostEqual(pf.free_energy(temp), -247.257535, 5)


    def test_classical1(self):
        pf = PartFun(
            NMA(load_molecule_g03fchk(
                pkg_resources.resource_filename(__name__, "../data/test/sterck/aa.fchk"))),
            [Vibrations(classical=True), ExtTrans(), ExtRot(1)],
        )

        for i in xrange(10):
            tmp = pf.vibrational.heat_capacity_terms(np.random.uniform(100,500))
            for value in tmp:
                self.assertAlmostEqual(value, boltzmann, 10)

    def test_classical2(self):
        pfc = PartFun(
            NMA(load_molecule_g03fchk(
                pkg_resources.resource_filename(__name__, "../data/test/sterck/aa.fchk"))),
            [Vibrations(True), ExtTrans(), ExtRot(1)],
        )
        pfq = PartFun(
            NMA(load_molecule_g03fchk(
                pkg_resources.resource_filename(__name__, "../data/test/sterck/aa.fchk"))),
            [Vibrations(False), ExtTrans(), ExtRot(1)],
        )

        lnzcs = pfc.vibrational.log_terms(50000)
        lnzqs = pfq.vibrational.log_terms(50000)
        for lnzc, lnzq in zip(lnzcs, lnzqs):
            self.assertAlmostEqual(lnzc, lnzq, 3)

        lnzcs = pfc.vibrational.logt_terms(50000)
        lnzqs = pfq.vibrational.logt_terms(50000)
        for lnzc, lnzq in zip(lnzcs, lnzqs):
            self.assertAlmostEqual(lnzc, lnzq, 7)

        lnzcs = pfc.vibrational.logtt_terms(50000)
        lnzqs = pfq.vibrational.logtt_terms(50000)
        for lnzc, lnzq in zip(lnzcs, lnzqs):
            self.assertAlmostEqual(lnzc, lnzq, 7)

    def test_external_separate(self):
        calmolK = calorie/mol/kelvin
        R = 1.98720649773 # R in calorie/(mol*K)
        temp = 298.15
        RT = R*temp*0.001 # RT in kcal/(mol*K)

        # values taken from aa.log
        molecule = load_molecule_g03fchk(
            pkg_resources.resource_filename(__name__, "../data/test/sterck/aa.fchk"))
        pf = PartFun(NMA(molecule, ConstrainExt()), [ExtTrans(cp=True)])
        # translational
        # WARNING: internal_heat returns enthalpy, turn it into internal energy
        self.assertAlmostEqual(pf.translational.internal_heat(temp)/(kcalmol)-RT, 0.889, 2)
        # WARNING: heat_capacity returns heat capacity at constant pressure,
        #          turn it into heat capacity at constant volume
        self.assertAlmostEqual(pf.translational.heat_capacity(temp)/(calmolK)-R, 2.981, 2)
        self.assertAlmostEqual(pf.translational.entropy(temp)/(calmolK), 38.699, 2)

        molecule = load_molecule_g03fchk(
            pkg_resources.resource_filename(__name__, "../data/test/sterck/aa.fchk"))
        pf = PartFun(NMA(molecule, ConstrainExt()), [ExtRot(1)])
        # rotational
        self.assertAlmostEqual(pf.rotational.internal_heat(temp)/(kcalmol), 0.889, 2)
        self.assertAlmostEqual(pf.rotational.heat_capacity(temp)/(calmolK), 2.981, 2)
        self.assertAlmostEqual(pf.rotational.entropy(temp)/(calmolK), 25.287, 2)

    def test_linear(self):
        molecule = load_molecule_g03fchk(
            pkg_resources.resource_filename(__name__, "../data/test/linear/gaussian.fchk"))
        pf = PartFun(NMA(molecule, ConstrainExt()), [ExtTrans(),ExtRot(2)])

        # check the inertia moments (eigenvalues only)
        iner_tens = pf.rotational.inertia_tensor
        evals, evecs = np.linalg.eigh(iner_tens)
        expected_evals = np.array([0.00000, 158.27937, 158.27937]) # from paats.log
        for i in 0,1,2:
            self.assertAlmostEqual(
                evals[i]/amu, expected_evals[i], 5,
                "Item %i of inertia tensor eigenvalues is wrong: %s!=%s" % (
                    i, evals[i]/amu, expected_evals[i]
                )
            )

        # check the natural logarithm rotational partition function
        self.assertAlmostEqual(5.607352, pf.rotational.log(298.150), 5)
        # check the count of the external rotation contribution
        self.assertEqual(pf.rotational.count, 2)
        # derived quantities
        calmolK = calorie/mol/kelvin
        self.assertAlmostEqual(pf.rotational.internal_heat(298.15)/(kcalmol), 0.592, 2)
        self.assertAlmostEqual(pf.rotational.heat_capacity(298.15)/(calmolK), 1.987, 2)
        self.assertAlmostEqual(pf.rotational.entropy(298.15)/(calmolK), 13.130, 2)

    def test_equilibrium(self):
        mol_react1 = load_molecule_g03fchk(
            pkg_resources.resource_filename(__name__, "../data/test/sterck/aa_1h2o_a.fchk"))
        mol_react2 = load_molecule_g03fchk(
            pkg_resources.resource_filename(__name__, "../data/test/sterck/aarad.fchk"))
        mol_complex = load_molecule_g03fchk(
            pkg_resources.resource_filename(__name__, "../data/test/sterck/paaprc_1h2o_b_aa.fchk"))
        mol_trans = load_molecule_g03fchk(
            pkg_resources.resource_filename(__name__, "../data/test/sterck/paats_1h2o_b_aa.fchk"))
        pf_react1 = PartFun(NMA(mol_react1), [ExtTrans(), ExtRot(1)])
        pf_react2 = PartFun(NMA(mol_react2), [ExtTrans(), ExtRot(1)])
        pf_complex = PartFun(NMA(mol_complex), [ExtTrans(), ExtRot(1)])
        pf_trans = PartFun(NMA(mol_trans), [ExtTrans(), ExtRot(1)])
        km1 = KineticModel([pf_react1, pf_react2], pf_trans)
        km2 = KineticModel([pf_complex], pf_trans)
        tm = ThermodynamicModel([pf_react1, pf_react2], [pf_complex])
        self.assertAlmostEqual(tm.unit*km2.unit, km1.unit)
        self.assertEqual(tm.unit_name, "m^3 mol^-1")
        self.assertEqual(km1.unit_name, "m^3 mol^-1 s^-1")
        self.assertEqual(km2.unit_name, "s^-1")

        for temp in np.arange(100,1000,10.0):
            K = tm.equilibrium_constant(temp)
            log_K = tm.equilibrium_constant(temp, do_log=True)
            self.assertAlmostEqual(np.log(K), log_K)
            k1 = km1.rate_constant(temp)
            log_k1 = km1.rate_constant(temp, do_log=True)
            self.assertAlmostEqual(np.log(k1), log_k1)
            k2 = km2.rate_constant(temp)
            log_k2 = km2.rate_constant(temp, do_log=True)
            self.assertAlmostEqual(np.log(k2), log_k2)
            self.assertAlmostEqual(log_K + log_k2, log_k1)

    def test_proton_cp(self):
        proton = Proton()
        pf = PartFun(NMA(proton), [ExtTrans(cp=True)])
        temp = 298.15
        self.assertAlmostEqual(pf.heat_capacity(temp), 2.5*boltzmann)
        self.assertAlmostEqual(pf.internal_heat(temp), 2.5*boltzmann*temp)
        self.assertAlmostEqual(pf.translational.helpert(temp, 1), 2.5)
        self.assertAlmostEqual(
            pf.translational.log(temp),
            1.5*np.log(proton.masses[0]*boltzmann*temp/(2*np.pi)) +
            np.log(boltzmann*temp/(1*atm))
        )

    def test_proton_cv(self):
        proton = Proton()
        pf = PartFun(NMA(proton), [ExtTrans(cp=False)])
        temp = 298.15
        self.assertAlmostEqual(pf.heat_capacity(temp), 1.5*boltzmann)
        self.assertAlmostEqual(pf.internal_heat(temp), 1.5*boltzmann*temp)
        self.assertAlmostEqual(pf.translational.helpert(temp, 1), 1.5)
        self.assertAlmostEqual(
            pf.translational.log(temp),
            1.0 + 1.5*np.log(proton.masses[0]*boltzmann*temp/(2*np.pi)) +
            np.log(meter**3/mol)
        )

    def test_pcm_correction(self):
        mol = load_molecule_g03fchk(
            pkg_resources.resource_filename(__name__, "../data/test/sterck/aa_1h2o_a.fchk"))
        pf0 = PartFun(NMA(mol), [])
        pf1 = PartFun(NMA(mol), [PCMCorrection((-5*kjmol,300))])
        pf2 = PartFun(NMA(mol), [PCMCorrection((-5*kjmol,300), (-10*kjmol,600))])
        # real tests
        self.assertAlmostEqual(pf0.free_energy(300)-5*kjmol, pf1.free_energy(300))
        self.assertAlmostEqual(pf0.free_energy(300)-5*kjmol, pf2.free_energy(300))
        self.assertAlmostEqual(pf0.free_energy(600)-5*kjmol, pf1.free_energy(600))
        self.assertAlmostEqual(pf0.free_energy(600)-10*kjmol, pf2.free_energy(600))
        # blind tests
        pf2.heat_capacity(300)
        with tmpdir(__name__, 'test_pcm_correction') as dn:
            pf1.write_to_file(os.path.join(dn, "pcm_partf1.txt"))
            pf2.write_to_file(os.path.join(dn, "pcm_partf2.txt"))

    def test_zero_point_energy(self):
        molecule = load_molecule_g03fchk(
            pkg_resources.resource_filename(__name__, "../data/test/ethane/gaussian.fchk"))
        nma = NMA(molecule)
        rotscan = load_rotscan_g03log(
            pkg_resources.resource_filename(__name__, "../data/test/rotor/gaussian.log"))
        rotor = Rotor(rotscan, molecule, rotsym=3, even=True)
        pf = PartFun(nma, [ExtTrans(), ExtRot(), rotor, Vibrations()])
        # check the vibrational part
        zpe_sum = 0.0
        zpes_vib = pf.vibrational.zero_point_energy_terms()
        freqs = pf.vibrational.positive_freqs
        for i in xrange(len(freqs)):
            self.assertAlmostEqual(0.5*freqs[i]*planck, zpes_vib[i])
            zpe_sum += zpes_vib[i]
        # check the electronic part
        self.assertAlmostEqual(pf.electronic.energy, pf.electronic.zero_point_energy())
        zpe_sum += pf.electronic.energy
        # check the other parts:
        self.assertAlmostEqual(0.0, pf.rotational.zero_point_energy())
        # just take the contribution from the rotor
        zpe_sum += rotor.zero_point_energy()
        # check the sum
        self.assertAlmostEqual(zpe_sum, pf.zero_point_energy())

    def test_freq_threshold(self):
        molecule = load_molecule_g03fchk(
            pkg_resources.resource_filename(__name__, "../data/test/ethane/gaussian.fchk"))
        nma = NMA(molecule)
        pf = PartFun(nma, [ExtTrans(), ExtRot(), Vibrations()])
        assert len(pf.vibrational.zero_freqs) == 6
        pf = PartFun(nma, [ExtTrans(), ExtRot(), Vibrations(freq_threshold=1e-3)])
        assert len(pf.vibrational.zero_freqs) == 12
        assert (pf.vibrational.positive_freqs > 1e-3).all()
