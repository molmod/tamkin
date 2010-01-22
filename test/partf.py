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

from molmod.constants import lightspeed, boltzmann
from molmod.units import cm, s, atm, amu, meter, mol, kcalmol, cal, K, kjmol

import unittest, numpy


__all__ = ["PartFunTestCase"]


class PartFunTestCase(unittest.TestCase):
    def test_phva_react_mat(self):
        fixed_atoms = load_fixed_g03com("input/mat/Zp_p_react.14mei.com")
        molecule = load_molecule_g03fchk("input/mat/Zp_p_react.28aug.fchk", "input/mat/Zp_p_react.14mei.fchk")
        nma = NMA(molecule, PHVA(fixed_atoms))
        pf = PartFun(nma)

        # values obtained with frektsjek.exe
        self.assertAlmostEqual(pf.log_eval(650), numpy.log(1.69428264E-45),3)
        self.assertAlmostEqual(pf.log_eval(700), numpy.log(3.26776842E-29),3)
        self.assertAlmostEqual(pf.log_eval(750), numpy.log(2.00896366E-14),3)

    def test_phva_rate_coeff_mat(self):
        fixed_atoms = load_fixed_g03com("input/mat/Zp_p_react.14mei.com")
        mol_react = load_molecule_g03fchk("input/mat/Zp_p_react.28aug.fchk", "input/mat/Zp_p_react.14mei.fchk")
        mol_trans = load_molecule_g03fchk("input/mat/Zp_p_TS.28aug.fchk", "input/mat/5Tp_p_TS.oniom21apr_HF.fchk")
        pf_react = PartFun(NMA(mol_react, PHVA(fixed_atoms)))
        pf_trans = PartFun(NMA(mol_trans, PHVA(fixed_atoms)))

        # values taken from the fancy excel file...
        temps = numpy.array([670,680,690,700,710,720,730,740,750,760,770])
        expected_ks = numpy.array([
            7.9473102E+05, 9.8300444E+05, 1.2085262E+06, 1.4771808E+06,
            1.7955340E+06, 2.1708793E+06, 2.6112829E+06, 3.1256298E+06,
            3.7236678E+06, 4.4160510E+06, 5.2143822E+06
        ])
        for i in xrange(len(temps)):
            k = compute_rate_coeff([pf_react], pf_trans, temps[i])
            self.assertAlmostEqual(numpy.log(k/(1/s)), numpy.log(expected_ks[i]),5)

    def test_gas_react_sterck(self):
        # Test both Full and ConstrainExt:
        for treatment, precision_wn, shift in (Full(), 0, 6), (ConstrainExt(), 3, 0):
            ## aa.fchk
            molecule = load_molecule_g03fchk("input/sterck/aa.fchk")
            nma = NMA(molecule, treatment)
            pf = PartFun(nma, [ExtTrans(), ExtRot(1)])

            # check partition function, values from aa.log
            self.assertAlmostEqual(16.973928, pf.translational.log_eval(298.150), 6)
            self.assertEqual(pf.rotational.count, 3)
            self.assertAlmostEqual(11.225093, pf.rotational.log_eval(298.150), 6)
            vib_contribs = pf.vibrational.log_eval_terms(298.150)
            expected_vib_contribs = numpy.array([
                0.674278, 0.292167, -0.357617, -1.017249, -1.018740, -1.427556,
                -1.428865
            ])
            for i in xrange(len(expected_vib_contribs)):
                self.assertAlmostEqual(vib_contribs[i], expected_vib_contribs[i], 2)
            self.assertAlmostEqual(-53.068692, pf.log_eval(298.150), 2)

            ## aarad.fchk
            molecule = load_molecule_g03fchk("input/sterck/aarad.fchk")
            nma = NMA(molecule, treatment)
            pf = PartFun(nma, [ExtTrans(), ExtRot(1)])

            # check partition function, values taken from aarad.log
            self.assertAlmostEqual(16.995059, pf.translational.log_eval(298.150), 6)
            self.assertEqual(pf.rotational.count, 3)
            self.assertAlmostEqual(11.319073, pf.rotational.log_eval(298.150), 6)
            vib_contribs = pf.vibrational.log_eval_terms(298.150)
            expected_vib_contribs = numpy.array([
                0.959168, 0.413328, -0.287477, -0.371573, -0.983851, -1.040910,
                -1.311721, -1.428436,
            ])
            for i in xrange(len(expected_vib_contribs)):
                self.assertAlmostEqual(vib_contribs[i], expected_vib_contribs[i], 2+precision_wn)
            self.assertAlmostEqual(-61.738525, pf.log_eval(298.150), 1+precision_wn)

    def test_gas_trans_sterck(self):
        # Test both Full and ConstrainExt:
        for treatment, precision_wn, shift in (Full(), 0, 6), (ConstrainExt(), 3, 0):
            ## paats.fchk
            molecule = load_molecule_g03fchk("input/sterck/paats.fchk")
            nma = NMA(molecule, treatment)
            pf = PartFun(nma, [ExtTrans(), ExtRot(1)])

            # check the natural logarithm partition function (split contributions):
            # (values taken from paats.log)
            self.assertAlmostEqual(18.024251, pf.translational.log_eval(298.150), 5)
            self.assertEqual(pf.rotational.count, 3)
            self.assertAlmostEqual(13.615243, pf.rotational.log_eval(298.150), 6)
            vib_contribs = pf.vibrational.log_eval_terms(298.150)
            expected_vib_contribs = numpy.array([
                1.749142, 1.335079, 0.789626, 0.704324, 0.365201, 0.061403,
                -0.060710, -0.281250, -0.309673, -0.423131, -0.447578, -0.983086,
                -1.017063, -1.107892, -1.192276, -1.318341, -1.352385, -1.427939,
            ])
            for i in xrange(len(expected_vib_contribs)):
                self.assertAlmostEqual(vib_contribs[i], expected_vib_contribs[i], 2+precision_wn)
            self.assertAlmostEqual(-139.302816, pf.log_eval(298.150), 1+precision_wn)

    def test_gas_rate_coeff_sterck(self):
        mol_react1 = load_molecule_g03fchk("input/sterck/aa.fchk")
        mol_react2 = load_molecule_g03fchk("input/sterck/aarad.fchk")
        mol_trans = load_molecule_g03fchk("input/sterck/paats.fchk")
        pf_react1 = PartFun(NMA(mol_react1, ConstrainExt()), [ExtTrans(), ExtRot(1)])
        pf_react2 = PartFun(NMA(mol_react2, ConstrainExt()), [ExtTrans(), ExtRot(1)])
        pf_trans = PartFun(NMA(mol_trans, ConstrainExt()), [ExtTrans(), ExtRot(1)])

        # values taken from the fancy excel file...
        temps = numpy.array([298.15,300,400,500,600,700,800,900,1000,1100])
        expected_ks = numpy.array([
            6.44881E-03,
            6.87398E-03, 9.74722E-02, 5.36970E-01, 1.82272E+00, 4.65134E+00,
            9.86811E+00, 1.84247E+01, 3.13480E+01, 4.97189E+01,
        ])
        unit = meter**3/mol/s
        for i in xrange(len(temps)):
            k = compute_rate_coeff([pf_react1, pf_react2], pf_trans, temps[i])
            # Sometimes, the fancy excel files use slightly different constants.
            # Therefore, only equal up to 2 decimals.
            self.assertAlmostEqual(numpy.log(k/unit), numpy.log(expected_ks[i]),2)

    def test_derivatives(self):
        molecule = load_molecule_g03fchk("input/ethane/gaussian.fchk")
        nma = NMA(molecule)
        rotscan = load_rotscan_g03log("input/rotor/gaussian.log")
        rotor = Rotor(rotscan, molecule, rotsym=3, even=True)
        pf = PartFun(nma, [
            ExtTrans(), ExtRot(6),
            rotor,
        ])

        eps = 0.0001
        temps = numpy.array([300.0,400.0,500.0,600.0,700.0])
        sfs = [
            pf.electronic, pf.translational, pf.rotational, pf.vibrational,
            pf.hindered_rotor_3_4_5, pf
        ]
        for stat_fys in sfs:
            for temp in temps:
                # check the first derivative with finite differences
                a = stat_fys.log_deriv(temp)
                b = (stat_fys.log_eval(temp+eps) - stat_fys.log_eval(temp-eps))/(2*eps)
                self.assertAlmostEqual(
                    a, b, 8,
                    "error in partial derivative (%s): %s!=%s" % (stat_fys.name, a, b)
                )
                # check the second derivative with finite differences
                a = stat_fys.log_deriv2(temp)
                b = (stat_fys.log_deriv(temp+eps) - stat_fys.log_deriv(temp-eps))/(2*eps)
                self.assertAlmostEqual(
                    a, b, 8,
                    "error in second partial derivative (%s): %s!=%s" % (stat_fys.name, a, b)
                )
                # check the helper functions temperature argument
                self.assertAlmostEqual(stat_fys.helper0(temp,1), stat_fys.helper0(temp,0)*temp)
                self.assertAlmostEqual(stat_fys.helper1(temp,1), stat_fys.helper1(temp,0)*temp)
                self.assertAlmostEqual(stat_fys.helper2(temp,1), stat_fys.helper2(temp,0)*temp)


    def test_derived_quantities(self):
        # internal energy, heat capacity and entropy
        molecule = load_molecule_g03fchk("input/sterck/aa.fchk")
        nma = NMA(molecule, ConstrainExt())
        pf = PartFun(nma, [ExtTrans(FixedVolume()), ExtRot(1)])

        # values taken from aa.log
        calmolK = cal/mol/K
        # gaussian seems to add the constant R to the translational part of the
        # entropy for some reason. I hope there is a mistake somewhere.
        R = 1.98720649773 # R in cal/(mol*K)
        # electronic
        self.assertAlmostEqual(pf.electronic.internal_energy(298.15)/(kcalmol), 0.000)
        self.assertAlmostEqual(pf.electronic.heat_capacity(298.15)/(calmolK), 0.000)
        self.assertAlmostEqual(pf.electronic.entropy(298.15)/(calmolK), 0.000)
        # translational
        self.assertAlmostEqual(pf.translational.internal_energy(298.15)/(kcalmol), 0.889, 2)
        self.assertAlmostEqual(pf.translational.heat_capacity(298.15)/(calmolK), 2.981, 2)
        self.assertAlmostEqual(pf.translational.entropy(298.15)/(calmolK), 38.699-R, 2) # corrected
        # rotational
        self.assertAlmostEqual(pf.rotational.internal_energy(298.15)/(kcalmol), 0.889, 2)
        self.assertAlmostEqual(pf.rotational.heat_capacity(298.15)/(calmolK), 2.981, 2)
        self.assertAlmostEqual(pf.rotational.entropy(298.15)/(calmolK), 25.287, 2)
        # vibrational
        self.assertAlmostEqual(pf.vibrational.internal_energy(298.15)/(kcalmol), 51.343, 2)
        self.assertAlmostEqual(pf.vibrational.heat_capacity(298.15)/(calmolK), 13.264, 2)
        self.assertAlmostEqual(pf.vibrational.entropy(298.15)/(calmolK), 10.710, 2)
        # total
        self.assertAlmostEqual((pf.internal_energy(298.15)-pf.energy)/(kcalmol), 53.121, 2)
        self.assertAlmostEqual(pf.heat_capacity(298.15)/(calmolK), 19.225, 2)
        self.assertAlmostEqual(pf.entropy(298.15)/(calmolK), 74.696, 2)
        # free energy of the molecule:
        self.assertAlmostEqual(pf.free_energy(298.15), -247.257535, 5)

    def test_classical(self):
        pf = PartFun(
            NMA(load_molecule_g03fchk("input/sterck/aa.fchk")),
            [Vibrations(classical=True), ExtTrans(), ExtRot(1)],
        )

        for i in xrange(10):
            tmp = pf.vibrational.heat_capacity_terms(numpy.random.uniform(100,500))
            for value in tmp:
                self.assertAlmostEqual(value, boltzmann, 10)

    def test_external_separate(self):
        calmolK = cal/mol/K
        R = 1.98720649773 # R in cal/(mol*K)

        # values taken from aa.log
        molecule = load_molecule_g03fchk("input/sterck/aa.fchk")
        pf = PartFun(NMA(molecule, ConstrainExt()), [ExtTrans(FixedVolume())])
        # translational
        self.assertAlmostEqual(pf.translational.internal_energy(298.15)/(kcalmol), 0.889, 2)
        self.assertAlmostEqual(pf.translational.heat_capacity(298.15)/(calmolK), 2.981, 2)
        self.assertAlmostEqual(pf.translational.entropy(298.15)/(calmolK), 38.699-R, 2) # corrected

        molecule = load_molecule_g03fchk("input/sterck/aa.fchk")
        pf = PartFun(NMA(molecule, ConstrainExt()), [ExtRot(1)])
        # rotational
        self.assertAlmostEqual(pf.rotational.internal_energy(298.15)/(kcalmol), 0.889, 2)
        self.assertAlmostEqual(pf.rotational.heat_capacity(298.15)/(calmolK), 2.981, 2)
        self.assertAlmostEqual(pf.rotational.entropy(298.15)/(calmolK), 25.287, 2)

    def test_linear(self):
        molecule = load_molecule_g03fchk("input/linear/gaussian.fchk")
        pf = PartFun(NMA(molecule, ConstrainExt()), [ExtTrans(),ExtRot(2)])

        # check the inertia moments (eigenvalues only)
        iner_tens = pf.rotational.inertia_tensor
        evals, evecs = numpy.linalg.eigh(iner_tens)
        expected_evals = numpy.array([0.00000, 158.27937, 158.27937]) # from paats.log
        for i in 0,1,2:
            self.assertAlmostEqual(
                evals[i]/amu, expected_evals[i], 5,
                "Item %i of inertia tensor eigenvalues is wrong: %s!=%s" % (
                    i, evals[i]/amu, expected_evals[i]
                )
            )

        # check the natural logarithm rotational partition function
        self.assertAlmostEqual(5.607352, pf.rotational.log_eval(298.150), 5)
        # check the count of the external rotation contribution
        self.assertEqual(pf.rotational.count, 2)
        # derived quantities
        calmolK = cal/mol/K
        self.assertAlmostEqual(pf.rotational.internal_energy(298.15)/(kcalmol), 0.592, 2)
        self.assertAlmostEqual(pf.rotational.heat_capacity(298.15)/(calmolK), 1.987, 2)
        self.assertAlmostEqual(pf.rotational.entropy(298.15)/(calmolK), 13.130, 2)

    def test_equilibrium(self):
        mol_react1 = load_molecule_g03fchk("input/sterck/aa_1h2o_a.fchk")
        mol_react2 = load_molecule_g03fchk("input/sterck/aarad.fchk")
        mol_complex = load_molecule_g03fchk("input/sterck/paaprc_1h2o_b_aa.fchk")
        mol_trans = load_molecule_g03fchk("input/sterck/paats_1h2o_b_aa.fchk")
        pf_react1 = PartFun(NMA(mol_react1), [ExtTrans(), ExtRot(1)])
        pf_react2 = PartFun(NMA(mol_react2), [ExtTrans(), ExtRot(1)])
        pf_comlex = PartFun(NMA(mol_complex), [ExtTrans(), ExtRot(1)])
        pf_trans = PartFun(NMA(mol_trans), [ExtTrans(), ExtRot(1)])

        for temp in numpy.arange(100,1000,10.0):
            log_K = compute_equilibrium_constant([pf_react1, pf_react2], [pf_comlex], temp)
            k = compute_rate_coeff([pf_react1, pf_react2], pf_trans, temp)
            log_k = numpy.log(k)
            k_prime = compute_rate_coeff([pf_comlex], pf_trans, temp)
            log_k_prime = numpy.log(k_prime)
            self.assertAlmostEqual(log_K + log_k_prime, log_k)

    def test_proton(self):
        mol = Proton()
        pf = PartFun(NMA(mol), [ExtTrans(FixedVolume())])
        temp = 298.15
        self.assertAlmostEqual(pf.heat_capacity(temp), 1.5*boltzmann)
        self.assertAlmostEqual(pf.internal_energy(temp), 1.5*boltzmann*temp)
        qt = (mol.masses[0]*boltzmann*temp/(2*numpy.pi))**1.5 * boltzmann*temp/(1*atm)
        self.assertAlmostEqual(pf.translational.helper1(temp, 1), 1.5)
        self.assertAlmostEqual(
            pf.translational.mol_volume.helper0(temp, 0),
            numpy.log(boltzmann*temp/(1*atm))
        )
        self.assertAlmostEqual(
            pf.translational.helper0(temp, 0),
            1.5*numpy.log(mol.masses[0]*boltzmann*temp/(2*numpy.pi)) +
            numpy.log(boltzmann*temp/(1*atm))
        )
        self.assertAlmostEqual(
            pf.entropy(temp)/boltzmann,
            (numpy.log(qt) + 1.5 + 1)
        )

    def test_pcm_correction(self):
        mol = load_molecule_g03fchk("input/sterck/aa_1h2o_a.fchk")
        pf0 = PartFun(NMA(mol), [
            ExtTrans(), ExtRot(1),
        ])
        pf1 = PartFun(NMA(mol), [
            ExtTrans(), ExtRot(1),
            PCMCorrection((-5*kjmol,300)),
        ])
        pf2 = PartFun(NMA(mol), [
            ExtTrans(), ExtRot(1),
            PCMCorrection((-5*kjmol,300), (-10*kjmol,600)),
        ])
        # real tests
        self.assertAlmostEqual(pf0.free_energy(300)-5*kjmol, pf1.free_energy(300))
        self.assertAlmostEqual(pf0.free_energy(300)-5*kjmol, pf2.free_energy(300))
        self.assertAlmostEqual(pf0.free_energy(600)-5*kjmol, pf1.free_energy(600))
        self.assertAlmostEqual(pf0.free_energy(600)-10*kjmol, pf2.free_energy(600))
        # blind tests
        pf1.heat_capacity(300)
        pf2.heat_capacity(300)
        pf1.write_to_file("output/pcm_partf1.txt")
        pf2.write_to_file("output/pcm_partf2.txt")


