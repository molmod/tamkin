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
import unittest

from molmod.constants import lightspeed, boltzmann
from molmod.units import centimeter, atm, amu, meter, mol, kcalmol
from molmod.test.common import tmpdir

from tamkin import *


__all__ = ["NMATestCase"]


class NMATestCase(unittest.TestCase):
    def check_ortho(self, modes):
        if len(modes) == 0:
            return
        unit_matrix = np.dot(modes.transpose(), modes)
        error_max = abs(unit_matrix - np.identity(len(unit_matrix))).max()
        self.assert_(error_max < 1e-5)

    def check_freqs(self, expected_freqs, nma, precision=3, check_zeros=False):
        """Check the frequencies in the partition function against expected values

        The expected values are given in 1/cm while the partition function works
        with atomic units.
        """
        if check_zeros:
            self.assertEqual(len(expected_freqs), len(nma.freqs))
        else:
            self.assertEqual(len(expected_freqs), len(nma.freqs)-len(nma.zeros))
        counter = 0
        for i in xrange(len(nma.freqs)):
            if check_zeros or (i not in nma.zeros):
                freq_in_cm = (nma.freqs[i]/lightspeed)/(1/centimeter)
                expected_freq = expected_freqs[counter]
                self.assertAlmostEqual(
                    expected_freq, freq_in_cm, precision,
                    "Frequency %i does not match, expected - computed = %.1f - %.1f = %f" % (i, expected_freq, freq_in_cm, expected_freq-freq_in_cm)
                )
                counter += 1

    def check_mode(self, expected_eig_mode, molecule, nma, index, precision=1):
        eig_mode = nma.modes[:,index].copy()
        eig_mode /= np.sqrt(molecule.masses3)
        eig_mode /= np.linalg.norm(eig_mode)
        self.assertEqual(len(expected_eig_mode), len(eig_mode))
        if np.dot(eig_mode, expected_eig_mode) < 0:
            expected_eig_mode *= -1
        for i in xrange(len(eig_mode)):
            self.assertAlmostEqual(
                expected_eig_mode[i], eig_mode[i], precision,
                "Component %i does not match, expected - computed = %.3f - %.3f = %e" % (i, expected_eig_mode[i], eig_mode[i], expected_eig_mode[i]-eig_mode[i])
            )

    def test_phva_react_mat(self):
        fixed_atoms = load_fixed_g03com(
            pkg_resources.resource_filename(__name__, "../data/test/mat/Zp_p_react.14mei.com"))
        molecule = load_molecule_g03fchk(
            pkg_resources.resource_filename(__name__, "../data/test/mat/Zp_p_react.28aug.fchk"),
            pkg_resources.resource_filename(__name__, "../data/test/mat/Zp_p_react.14mei.fchk"))
        nma = NMA(molecule, PHVA(fixed_atoms))
        self.check_ortho(nma.modes)

        # from Zp_p_react.28aug.log
        expected_freqs = np.array([ # values in 1/cm !!!
            27.8288, 37.3524, 43.9469, 63.0762, 66.1418, 69.2043, 75.3587,
            77.9613, 80.9633, 81.0899, 85.4568, 89.1801, 91.9183, 93.1096,
            94.3347, 97.5867, 98.0803, 101.7701, 102.8680, 105.0344, 108.0227,
            112.9728, 113.7285, 117.3181, 118.4344, 121.3627, 122.0002,
            123.9126, 124.4411, 125.6629, 129.4535, 131.5315, 132.3212,
            133.6929, 136.5911, 137.0378, 137.9674, 141.4904, 145.3214,
            147.8669, 148.1816, 149.9283, 150.9635, 153.2859, 153.8381,
            160.3900, 161.5367, 164.4291, 169.0828, 171.2338, 175.0859,
            176.5708, 183.0398, 184.0936, 186.1756, 188.3502, 190.9562,
            192.2252, 194.9651, 198.0290, 199.3246, 202.8058, 205.7420,
            207.6902, 212.4975, 219.7480, 219.9852, 221.2309, 224.6316,
            224.8149, 233.4607, 234.3927, 235.3028, 236.9118, 237.9797,
            241.7946, 243.4485, 246.1978, 247.9927, 249.0041, 250.1762,
            250.9420, 251.5648, 254.2162, 254.8433, 256.6327, 257.4815,
            258.4400, 258.7302, 259.6071, 262.8808, 263.8354, 265.8626,
            266.0468, 270.3915, 270.7264, 272.8667, 274.8286, 279.7044,
            280.3368, 281.1415, 281.6211, 282.8733, 284.0038, 285.3430,
            287.5569, 289.8842, 290.6086, 291.9417, 294.9712, 295.4538,
            296.4495, 297.8204, 300.0798, 301.0852, 302.6940, 303.4360,
            304.2641, 304.9022, 307.2555, 307.7802, 309.6201, 311.8870,
            313.6320, 315.0156, 318.3005, 321.1154, 323.0235, 324.8000,
            327.7217, 329.5247, 330.5125, 331.8464, 336.1074, 339.3014,
            340.4516, 342.9348, 346.8101, 348.9697, 352.3876, 354.7276,
            355.9487, 356.8001, 358.9632, 366.4523, 367.0573, 369.5507,
            370.6650, 375.7956, 380.6060, 382.2574, 383.8524, 386.8728,
            388.6802, 391.1230, 395.7112, 396.2303, 400.3917, 402.2810,
            406.1397, 406.5577, 409.8623, 412.0959, 415.6868, 416.5133,
            417.3426, 421.6224, 423.9824, 424.4674, 427.6703, 430.2451,
            434.5627, 437.5288, 438.6750, 439.7655, 445.7135, 447.8034,
            448.4360, 453.2775, 453.8330, 459.5623, 462.9828, 468.6418,
            470.2151, 471.8227, 473.7062, 474.0171, 479.7171, 480.7117,
            482.0891, 483.9965, 490.2403, 493.4729, 494.8948, 498.6155,
            499.3765, 501.4588, 502.3205, 505.5651, 507.5738, 511.6108,
            512.8027, 514.4808, 520.2326, 521.3976, 525.9624, 530.2795,
            533.3707, 536.4974, 569.3762, 571.7824, 573.5930, 579.3130,
            580.5857, 587.4812, 593.6915, 595.8546, 596.7363, 597.6287,
            606.7116, 608.0643, 609.3139, 613.2829, 615.4096, 620.1445,
            627.2731, 629.1332, 636.9593, 640.8470, 645.5420, 662.3197,
            701.7831, 712.0490, 720.4640, 722.3991, 724.2978, 725.8450,
            728.2839, 732.6112, 740.3948, 741.9380, 743.4002, 746.4851,
            747.1192, 751.2785, 752.5632, 752.8419, 756.3768, 759.1793,
            760.2378, 763.6799, 766.5302, 768.6952, 773.5854, 777.1931,
            783.2729, 787.8810, 788.3535, 793.9564, 796.4370, 798.7310,
            801.9063, 804.9096, 808.6233, 808.8822, 810.0362, 811.9579,
            816.3942, 819.6440, 821.6829, 822.7110, 824.5597, 826.3847,
            830.0432, 832.0247, 836.1598, 837.7947, 838.4538, 839.8984,
            844.3143, 847.7515, 849.2375, 851.1902, 852.8806, 854.5031,
            855.0160, 858.5909, 862.2544, 867.9239, 869.2610, 871.3386,
            877.1552, 887.9678, 896.7647, 929.4061, 934.2894, 941.1583,
            957.1980, 1033.3876, 1062.8739, 1072.0005, 1084.6886, 1086.4886,
            1115.5630, 1120.0973, 1141.1535, 1154.2137, 1189.6590, 1196.6184,
            1200.7893, 1203.9797, 1206.1618, 1208.3705, 1211.4513, 1211.8981,
            1216.5096, 1220.5589, 1221.3868, 1222.6332, 1224.5551, 1225.9007,
            1228.8901, 1232.0895, 1233.4362, 1234.1402, 1236.3802, 1238.1117,
            1240.4172, 1264.1351, 1264.6328, 1269.8886, 1271.4419, 1272.5889,
            1273.0191, 1275.6729, 1276.7243, 1279.7357, 1281.2038, 1282.1332,
            1283.2891, 1285.3206, 1285.9726, 1287.2464, 1288.9866, 1289.7828,
            1292.3464, 1292.9224, 1293.6931, 1294.1804, 1295.7870, 1297.3050,
            1297.5375, 1300.0006, 1301.9643, 1303.9580, 1305.2959, 1306.2437,
            1307.8970, 1308.4689, 1309.4465, 1311.6425, 1312.6093, 1313.5695,
            1316.2741, 1317.5007, 1321.8268, 1326.6867, 1330.0227, 1333.1364,
            1334.3544, 1334.8633, 1335.8174, 1338.8094, 1340.0368, 1401.4307,
            1430.9737, 1436.8405, 1443.6757, 1468.4197, 1499.4699, 1506.8307,
            1517.0168, 1523.8745, 1531.1322, 1533.3524, 1715.3981, 3038.9474,
            3046.7532, 3071.6163, 3085.2573, 3126.1789, 3129.6282, 3140.8569,
            3146.8763, 3150.5051, 3153.8649, 3187.9663, 3215.3556, 3231.9653
        ])
        self.check_freqs(expected_freqs, nma)

        # eigenmode 360 (counting from zero) from Zp_p_react.28aug.log
        expected_eig_mode = np.array([
            0, -0.01, -0.01, 0, -0.01, 0, 0, 0, 0, 0, 0.03, 0.01, 0, 0, -0.01,
            0, 0.01, 0, -0.01, 0, 0, 0, -0.01, -0.01, 0, 0.03, -0.01, 0, -0.01,
            0.01, 0, 0, -0.03, 0.02, -0.03, 0.01, -0.02, 0, 0.01, 0.02, -0.01,
            0, 0, 0, 0.01, -0.02, 0.01, 0, 0.02, -0.03, 0, 0, 0.01, 0.03, 0, 0,
            0, 0, -0.02, 0.01, 0, 0, 0.01, -0.01, 0, -0.01, 0.04, -0.02, 0.01,
            0.01, 0.03, -0.04, -0.06, -0.09, 0.12, 0.01, 0.05, -0.02, -0.02,
            -0.04, -0.05, 0, 0.03, 0.01, 0, 0, 0, 0, -0.03, 0.01, -0.02, 0.03,
            0, 0.03, -0.05, -0.01, 0.09, -0.01, -0.03, 0, 0.01, 0.19, -0.09, 0,
            -0.03, 0.37, -0.12, -0.09, -0.07, 0.07, 0.06, -0.15, -0.06, -0.05,
            0.06, 0.04, -0.02, 0, 0, 0.06, -0.08, -0.04, -0.03, 0.18, 0.07,
            -0.05, 0.05, -0.07, 0.06, -0.34, 0.1, -0.1, 0.04, 0.01, 0.01, 0, 0,
            -0.03, 0, 0, 0.01, 0, -0.02, 0, 0, 0, -0.09, 0, -0.01, 0.02, -0.03,
            0.04, 0, -0.01, -0.01, 0, -0.01, -0.11, -0.02, -0.03, 0.02, 0.01,
            0, 0, 0.05, 0.04, -0.01, 0.01, 0.01, 0.07, -0.02, 0.1, 0.01, -0.11,
            -0.03, 0.02, 0.04, 0.03, 0, 0.03, -0.01, 0.01, -0.01, 0, -0.03,
            0.01, 0.02, -0.03, 0, 0, 0, 0, 0.01, 0.01, -0.02, -0.02, -0.03,
            0.01, 0.01, 0.01, 0, -0.01, 0, 0, 0.01, 0.01, 0, 0, -0.05, -0.01,
            0.02, 0.01, 0.01, -0.02, -0.01, 0, 0.01, 0.02, 0.02, -0.04, 0.08,
            -0.03, 0, 0.02, 0, -0.05, -0.07, 0.09, 0.02, -0.03, -0.02, 0.01,
            0.14, -0.03, -0.01, -0.04, 0.01, 0, 0.01, 0, 0.01, 0, 0, 0, -0.01,
            -0.01, 0, 0, 0.02, 0, 0, -0.06, 0, -0.02, 0.03, -0.02, 0.03, 0,
            -0.01, -0.01, 0.01, 0, 0, 0.02, 0, 0, -0.01, 0.01, -0.01, 0, -0.04,
            0.01, -0.01, -0.11, -0.02, -0.12, 0.04, -0.02, 0.02, 0.04, 0.06,
            0.02, -0.02, -0.02, -0.01, 0.02, 0, 0, -0.02, -0.02, 0.02, 0.01,
            0.06, -0.01, 0, -0.02, 0, -0.01, 0, 0, 0, -0.01, 0.02, 0.01, -0.01,
            0, 0.02, 0.09, 0.02, -0.03, -0.3, -0.06, -0.02, 0.07, -0.03, 0.05,
            0.09, 0.11, -0.08, -0.01, 0.08, 0.02, 0.01, -0.03, 0, 0.01, 0, 0,
            -0.01, 0, 0, 0, 0, 0.02, 0.04, -0.05, 0, 0, 0, 0.03, 0.01, -0.04,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.01, 0.02, -0.01, -0.01,
            -0.03, 0.01, 0.02, 0.01, 0.01, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.09,
            0.07, -0.02, -0.2, -0.08, 0.01, -0.03, 0.25, 0.03, -0.03, -0.11,
            0.01, -0.06, 0.01, 0.01, -0.04, 0.05, 0, 0.02, -0.07, 0, 0, 0, 0,
            -0.01, -0.01, -0.01, 0.01, 0.01, 0.02, 0, 0, 0, 0, 0, 0, 0, 0, 0
        ])
        self.check_mode(expected_eig_mode, molecule, nma, 360)

    def test_phva_trans_mat(self):
        fixed_atoms = load_fixed_g03com(
            pkg_resources.resource_filename(__name__, "../data/test/mat/Zp_p_TS.28aug.com"))
        molecule = load_molecule_g03fchk(
            pkg_resources.resource_filename(__name__, "../data/test/mat/Zp_p_TS.28aug.fchk"),
            pkg_resources.resource_filename(__name__, "../data/test/mat/5Tp_p_TS.oniom21apr_HF.fchk"))
        nma = NMA(molecule, PHVA(fixed_atoms))
        self.check_ortho(nma.modes)

        # from Zp_p_TS.28aug.log
        expected_freqs = np.array([
            -277.5011, 23.9092, 41.7120, 61.9660, 63.0962, 69.8307, 73.7926,
            75.0049, 76.9176, 80.6721, 84.2794, 88.8336, 92.1835, 94.9068,
            98.0852, 99.9322, 103.1190, 103.7393, 107.8247, 110.2122, 115.7156,
            117.5809, 118.2052, 121.0087, 122.4923, 123.5816, 124.0811,
            126.9427, 127.7703, 130.6091, 131.3456, 132.3588, 135.2500,
            136.7738, 138.6681, 140.8074, 144.4676, 146.4335, 148.3035,
            149.9876, 150.4780, 150.6423, 152.0692, 154.5734, 155.4264,
            159.0561, 162.2418, 165.1299, 170.0755, 171.2581, 177.1458,
            182.7523, 184.9913, 185.8296, 189.4864, 190.3428, 192.0965,
            195.7371, 196.1645, 197.6649, 200.2669, 202.2563, 205.8011,
            211.3010, 213.0861, 214.9161, 219.4970, 221.4366, 224.1637,
            224.9258, 225.7103, 233.4629, 234.0701, 234.8032, 236.7302,
            237.4706, 240.5476, 242.7368, 247.6557, 248.5794, 250.9117,
            251.0468, 252.7569, 254.2318, 255.7725, 256.8110, 257.6393,
            257.7811, 258.6360, 259.5856, 262.9776, 263.4494, 266.4263,
            269.4301, 270.4638, 272.2295, 276.9440, 279.4317, 280.0615,
            280.7915, 281.8451, 282.6262, 283.1169, 285.0517, 287.4292,
            289.0002, 289.8826, 291.5338, 292.5431, 293.8636, 296.6972,
            297.1196, 299.4641, 300.4013, 301.3917, 302.9005, 303.7075,
            304.7195, 305.9174, 307.5636, 307.9825, 309.9637, 310.3563,
            312.6902, 315.4988, 317.3288, 320.9580, 321.8814, 323.9511,
            327.4935, 328.4590, 331.1124, 332.7342, 338.4337, 339.8494,
            341.3357, 344.9356, 347.1832, 349.9354, 350.9264, 352.6791,
            356.0405, 356.9244, 359.6489, 364.8860, 367.7483, 369.1119,
            371.2764, 375.2813, 379.2187, 381.1047, 384.3394, 385.9994,
            387.6509, 389.5949, 392.6860, 396.4177, 398.7416, 402.9895,
            405.6497, 406.3478, 408.3214, 409.8340, 412.0110, 415.2781,
            416.0211, 417.2053, 418.8546, 422.1958, 424.1089, 424.6220,
            430.8081, 433.0297, 435.3051, 436.9777, 438.4862, 442.1633,
            447.5428, 449.5805, 451.7231, 452.4786, 457.1414, 460.0143,
            462.8036, 470.4479, 471.2113, 473.3745, 474.8016, 475.4379,
            480.1118, 482.3782, 483.0683, 484.4649, 491.8352, 494.2355,
            498.2456, 499.7252, 501.8871, 502.7320, 503.6131, 507.0049,
            511.9813, 513.5206, 513.6919, 520.2267, 523.7444, 526.4146,
            529.8727, 534.4395, 535.7987, 570.7068, 574.4993, 579.9835,
            580.5928, 582.1656, 592.5232, 594.1891, 594.8747, 599.2579,
            607.6835, 608.8170, 610.0326, 614.3049, 616.1731, 620.0777,
            628.3081, 629.3017, 638.0315, 642.8964, 661.8123, 662.3169,
            687.4683, 701.0624, 716.8334, 721.1044, 723.0386, 725.6946,
            727.9533, 731.3815, 732.7977, 739.3924, 742.0206, 743.4530,
            745.3309, 747.0663, 748.1048, 751.6470, 752.1440, 753.1612,
            755.9626, 759.5729, 762.5541, 763.7598, 769.2407, 772.3410,
            777.1195, 782.3273, 782.6473, 787.8514, 794.1608, 795.1825,
            800.2166, 801.1988, 802.8549, 807.0825, 808.5573, 811.3337,
            812.8680, 814.9577, 819.1871, 820.4074, 822.0661, 823.9108,
            826.0580, 828.1174, 833.1488, 835.0466, 836.7346, 838.5173,
            843.5544, 845.7829, 846.7962, 849.6257, 851.7854, 852.8993,
            855.3334, 856.9111, 859.4487, 865.6051, 868.4866, 870.5800,
            873.4649, 878.7202, 924.7622, 930.2570, 945.0211, 960.4800,
            988.6512, 1030.5057, 1044.3807, 1065.9849, 1074.8552, 1076.6902,
            1086.9965, 1103.2362, 1105.6620, 1118.2368, 1152.6625, 1175.1914,
            1194.4755, 1197.2553, 1199.9105, 1207.7522, 1207.8424, 1209.3052,
            1210.0464, 1211.4468, 1215.8036, 1217.3347, 1218.8982, 1223.2471,
            1224.0862, 1225.3685, 1228.0312, 1229.0827, 1230.7054, 1231.2875,
            1234.4712, 1235.8251, 1238.3234, 1259.4533, 1261.6019, 1267.6963,
            1268.0374, 1270.8954, 1271.4609, 1273.5123, 1275.1139, 1275.8049,
            1278.2507, 1280.0891, 1282.0863, 1284.6757, 1285.7967, 1286.6318,
            1286.9970, 1288.2161, 1289.6203, 1290.6982, 1291.8129, 1295.5600,
            1296.4923, 1298.0801, 1298.8114, 1299.7572, 1300.3659, 1302.0762,
            1303.3364, 1305.1519, 1306.1110, 1307.7560, 1308.6144, 1309.8298,
            1310.3611, 1312.5112, 1313.7689, 1318.1737, 1326.5405, 1328.3189,
            1330.5304, 1331.6633, 1332.0576, 1336.1644, 1337.0547, 1338.3461,
            1410.7321, 1434.5536, 1442.6381, 1466.0699, 1471.7556, 1496.3187,
            1503.4529, 1521.2544, 1525.1020, 1535.3449, 1669.8518, 3033.2348,
            3054.8378, 3066.5299, 3077.7394, 3110.8884, 3133.0316, 3144.0042,
            3152.0511, 3162.7781, 3176.4155, 3213.4404, 3251.0431, 3327.1758
        ])
        self.check_freqs(expected_freqs, nma)

    def test_gas_react_sterck(self):
        # Test both Full and ConstrainExt:
        for treatment, precision_wn, shift in (Full(), 0, 6), (ConstrainExt(), 3, 0):
            ## aa.fchk
            molecule = load_molecule_g03fchk(
                pkg_resources.resource_filename(__name__, "../data/test/sterck/aa.fchk"))
            nma = NMA(molecule, treatment)
            self.check_ortho(nma.modes)

            # expected frequencies from aa.log
            expected_freqs = np.array([
                104.4761, 151.3376, 275.5587, 467.4732, 467.9741, 613.6715,
                614.1609, 814.9387, 819.7633, 997.5119, 1014.5531, 1043.8629,
                1118.1225, 1297.5031, 1363.9344, 1455.7262, 1644.8856, 1697.0197,
                1764.3581, 3159.7250, 3175.2247, 3261.5878, 3589.5604, 3717.1382,
            ])
            self.check_freqs(expected_freqs, nma, precision_wn)

            # expected eigenmode (4) from aa.log, counting from zero
            expected_eig_mode = np.array([
                0.03, -0.01, -0.05, 0.04, 0.00, 0.02, 0.05, -0.01, 0.38, 0.03,
                0.02, 0.51, 0.07, 0.00, -0.36, 0.00, 0.02, -0.13, -0.03, 0.02,
                0.05, -0.04, -0.02, 0.02, -0.12, -0.02, -0.24, 0.01, -0.11, 0.59,
            ])
            self.check_mode(expected_eig_mode, molecule, nma, 4+shift)

            ## aarad.fchk
            molecule = load_molecule_g03fchk(
                pkg_resources.resource_filename(__name__, "../data/test/sterck/aarad.fchk"))
            nma = NMA(molecule, treatment)

            # expected frequencies from aa.log
            expected_freqs = np.array([
                78.9330, 134.6847, 259.0440, 278.9383, 456.3221, 475.4515, 570.8836,
                614.0007, 719.9878, 822.1945, 1000.9025, 1039.4987, 1091.4071,
                1158.5616, 1283.5503, 1416.2231, 1451.5736, 1496.9929, 1510.7350,
                1637.4175, 1690.9547, 3008.8989, 3041.7365, 3144.6983, 3169.1683,
                3578.1294, 3701.0701,
            ])
            self.check_freqs(expected_freqs, nma, precision_wn)

            # expected eigenmode (10) from aarad.log, counting from zero
            expected_eig_mode = np.array([
                0.00, 0.00, -0.11, 0.00, 0.00, 0.14, 0.00, 0.00, 0.11, 0.59, 0.22,
                -0.19, 0.00, -0.01, -0.31, 0.00, 0.00, 0.03, 0.00, 0.00, -0.01,
                0.00, 0.00, 0.00, 0.02, 0.01, -0.01, 0.00, 0.00, 0.00, -0.58,
                -0.23, -0.18,
            ])
            self.check_mode(expected_eig_mode, molecule, nma, 10+shift)

    def test_teller_redlich_gas_react_sterck(self):
        for treatment, precision_wn, shift in (Full(), 0, 6), (ConstrainExt(), 3, 0):
            mol1 = load_molecule_g03fchk(
                pkg_resources.resource_filename(__name__, "../data/test/sterck/aa.fchk"))
            nma1 = NMA(mol1, treatment)
            self.check_ortho(nma1.modes)

            mol2 = mol1.copy_with(masses=mol1.masses*np.random.uniform(0.9,1.1,mol1.size))
            nma2 = NMA(mol2, treatment)
            self.check_ortho(nma2.modes)

            ratio12 = 1.0
            for i in xrange(len(nma1.freqs)):
                if i not in nma1.zeros:
                    ratio12 *= nma1.freqs[i]
                if i not in nma2.zeros:
                    ratio12 /= nma2.freqs[i]

            check = (nma1.masses.sum()/nma2.masses.sum())**(1.5)
            imoms1 = np.linalg.eigvalsh(nma1.inertia_tensor)
            imoms2 = np.linalg.eigvalsh(nma2.inertia_tensor)
            check *= (imoms1.prod()/imoms2.prod())**(0.5)
            check /= (nma1.masses/nma2.masses).prod()**(1.5)
            self.assertAlmostEqual(ratio12, check, 2)


    def test_gas_trans_sterck(self):
        # Test both Full and ConstrainExt:
        for treatment, precision_wn, shift in (Full(), 0, 6), (ConstrainExt(), 3, 0):
            ## paats.fchk
            molecule = load_molecule_g03fchk(
                pkg_resources.resource_filename(__name__, "../data/test/sterck/paats.fchk"))
            nma = NMA(molecule, treatment)
            self.check_ortho(nma.modes)

            # expected frequencies from paats.log
            expected_freqs = np.array([
                -413.6054, 35.9959, 54.3725, 93.2932, 101.4448, 141.0851, 188.3344,
                210.9659, 257.6158, 264.1854, 291.6926, 297.8879, 456.0680,
                467.4108, 498.3737, 527.9183, 573.3009, 585.7877, 613.8148,
                678.1017, 719.9223, 790.7100, 802.0216, 833.6055, 837.8822,
                946.2926, 966.0422, 1039.5715, 1063.6764, 1068.3896, 1115.9426,
                1125.0291, 1155.6926, 1285.8051, 1295.6331, 1325.5446, 1421.1412,
                1441.6615, 1458.0357, 1495.3439, 1510.1462, 1568.3896, 1640.5548,
                1658.4111, 1709.9513, 1723.4623, 3011.8044, 3082.1420, 3135.2699,
                3170.8738, 3177.6305, 3180.7047, 3269.1535, 3528.7377, 3579.7690,
                3666.7518, 3704.6548,
            ])
            self.check_freqs(expected_freqs, nma, precision_wn)

            # expected eigenmode (20 from paats.log, counting from zero
            expected_eig_mode = np.array([
                0.01, 0.00, 0.01, 0.00, 0.00, 0.02, 0.00, 0.01, 0.00, 0.00, -0.01,
                0.00, -0.01, 0.00, 0.00, -0.03, -0.02, -0.05, -0.01, 0.00, -0.01,
                0.28, 0.17, 0.24, -0.08, -0.04, -0.04, -0.08, -0.05, -0.10, 0.00,
                0.02, -0.01, -0.02, -0.03, 0.01, -0.04, -0.08, 0.00, -0.01, 0.00,
                0.01, -0.09, -0.11, -0.01, -0.42, -0.39, -0.23, -0.28, -0.18, -0.38,
                0.25, 0.17, 0.18, 0.04, 0.00, 0.02, -0.05, 0.08, 0.00, 0.03, 0.00,
                0.00,
            ])
            self.check_mode(expected_eig_mode, molecule, nma, 20+shift)

            # expected eigenmode (0) from paats.log, counting from zero
            expected_eig_mode = np.array([
                0.47, 0.33, 0.21, 0.07, 0.06, 0.06, 0.03, 0.00, -0.01, 0.01, 0.00,
                0.01, 0.01, -0.01, -0.01, -0.44, -0.29, -0.26, -0.09, -0.01, -0.02,
                -0.09, -0.04, -0.06, 0.01, -0.01, 0.00, 0.01, -0.02, 0.04, 0.00,
                -0.06, 0.01, 0.00, -0.03, 0.00, 0.04, -0.05, 0.04, 0.03, -0.15,
                -0.01, 0.07, 0.01, 0.06, -0.03, 0.01, -0.07, -0.02, 0.01, -0.01,
                0.10, 0.07, 0.11, 0.09, -0.02, 0.01, -0.28, 0.25, 0.10, 0.04, -0.05,
                0.01,
            ])
            self.check_mode(expected_eig_mode, molecule, nma, 0)

            # check the inertia tensor (eigenvalues and eigenvectors)
            # values taken from paats.log
            iner_tens = nma.inertia_tensor
            evals, evecs = np.linalg.eigh(iner_tens)
            expected_evals = np.array([875.32578,2288.50418,2609.96955])
            for i in 0,1,2:
                self.assertAlmostEqual(
                    evals[i]/amu, expected_evals[i], 5,
                    "Item %i of inertia tensor eigenvalues is wrong: %s!=%s" % (
                        i, evals[i]/amu, expected_evals[i]
                    )
                )
            expected_evecs = np.array([
                [ 0.99902, -0.03999, 0.01883],
                [ 0.03964,  0.99904, 0.01875],
                [-0.01956, -0.01799, 0.99965],
            ])
            self.assert_(abs(evecs-expected_evecs).max() < 1e-3)

    def test_gas_pentane(self):
        molecule = load_molecule_cp2k(
            pkg_resources.resource_filename(__name__, "../data/test/cp2k/pentane/sp.out"),
            pkg_resources.resource_filename(__name__, "../data/test/cp2k/pentane/freq.out"),
            is_periodic=False)
        nma = NMA(molecule, ConstrainExt(), do_modes=False)

        expected_freqs = np.array([ # taken from cp2k output file
            125.923216, 132.737491, 251.881675, 258.656045, 260.028863,
            489.690362, 522.430525, 836.589205, 924.800235, 1069.391906,
            1152.109253, 1182.947666, 1226.649618, 1279.047068, 1494.291889,
            1505.963333, 1526.247839, 1526.770449, 1527.220869, 1554.895642,
            1574.963368, 1664.772539, 1726.540705, 1730.035399, 1730.370082,
            1740.966064, 1740.971286, 1758.656380, 1767.704604, 1790.726516,
            1825.670669, 1980.257356, 2107.134204, 4515.156994, 4526.218064,
            4540.443901, 4578.619346, 4578.624300, 4581.926079, 4582.125771,
            4646.530637, 4657.432936, 4671.093295, 4751.240132, 4751.291450
        ])
        self.check_freqs(expected_freqs, nma, 1)

    def test_gas_water_cpmd(self):
        molecule = load_molecule_cpmd(
            pkg_resources.resource_filename(__name__, "../data/test/cpmd/damp.out"),
            pkg_resources.resource_filename(__name__, "../data/test/cpmd/GEOMETRY.xyz"),
            pkg_resources.resource_filename(__name__, "../data/test/cpmd/MOLVIB"),
            is_periodic=True)
        nma = NMA(molecule, Full(), do_modes=False)
        expected_freqs = np.array([
            -333.5480, -86.0913, -39.3998, 49.2809, 63.4021, 190.1188,
            1595.5610, 3724.2110, 3825.7550
        ])
        self.check_freqs(expected_freqs, nma, 4, check_zeros=True)

        nma = NMA(molecule, ConstrainExt(), do_modes=False)
        expected_freqs = np.array([
            -119.052347919, -80.8279424581, 51.861869587, 1595.55084877, 3724.21041319, 3825.75442327 # tamkin values
            #-664.0947, -83.1893, 128.0597, 1620.6557, 3730.2012, 3831.2725 # cpmd values
        ])
        self.check_freqs(expected_freqs, nma, 4)

    def test_vsa(self):
        molecule = load_molecule_charmm(
            pkg_resources.resource_filename(__name__, "../data/test/an/ethanol.cor"),
            pkg_resources.resource_filename(__name__, "../data/test/an/ethanol.hess.full"))

        #  --- single atom as subsystem: results should be three translations
        subs = load_indices(
            pkg_resources.resource_filename(__name__, "../data/test/an/fixed.01.txt"))
        nma = NMA(molecule, VSA(subs))
        self.check_ortho(nma.modes)
        self.assert_(len(nma.zeros)==3)
        expected_freqs = np.array([-0.23313896,  0.01475299,  0.06910995])
        self.check_freqs(expected_freqs, nma, 0, check_zeros=True)

        #  ---  atoms of subsystem are collinear: not yet external_basis implemented...
        subs = load_indices(
            pkg_resources.resource_filename(__name__, "../data/test/an/fixed.02.txt"))
        nma = NMA(molecule, VSA(subs))
        self.check_ortho(nma.modes)
        self.assert_(len(nma.zeros)==5)
        expected_freqs = np.array([ -3.30245504e-01,-2.51869284e-01, -1.16805787e-01,
               1.37273058e-01,1.82394281e-01,1.05116208e+03])
        self.check_freqs(expected_freqs, nma, 0, check_zeros=True)

        #  --- atoms of subsystem are not collinear
        subs = load_indices(
            pkg_resources.resource_filename(__name__, "../data/test/an/fixed.03.txt"))  # atom 1 to atom 7
        nma = NMA(molecule, VSA(subs))
        self.check_ortho(nma.modes)
        self.assert_(len(nma.zeros)==6)
        expected_freqs = np.array([ # taken from previous working python version
               -0.00962070029417,-0.00653077581838,-0.000469690606812,-0.000148860724889,
               0.000379327735493, 0.0142496837765,   # These are the zeros.
               302.048797009,1001.64994625,1045.10286172,1232.6538326,2814.79944849,
               3653.7642078 ])
        self.check_freqs(expected_freqs, nma, 4, check_zeros=True)

    def test_vsa_no_mass(self):
        # Modes are a priori known to be non-orthogonal, so no 'self.check_ortho(nma.modes)'
        molecule = load_molecule_charmm(
            pkg_resources.resource_filename(__name__, "../data/test/an/ethanol.cor"),
            pkg_resources.resource_filename(__name__, "../data/test/an/ethanol.hess.full"))

        #  --- single atom as subsystem: results should be three translations
        subs = load_indices(
            pkg_resources.resource_filename(__name__, "../data/test/an/fixed.01.txt"))
        nma = NMA(molecule, VSANoMass(subs))
        self.assert_(len(nma.zeros)==3)
        expected_freqs = np.array([-0.4205594, 0.03940166, 0.13774798])
        self.check_freqs(expected_freqs, nma, 0, check_zeros=True)

        #  ---  atoms of subsystem are collinear
        subs = load_indices(
            pkg_resources.resource_filename(__name__, "../data/test/an/fixed.02.txt"))
        nma = NMA(molecule, VSANoMass(subs))
        self.assert_(len(nma.zeros)==5)
        expected_freqs = np.array([-6.07975753e-01,-3.47371852e-01,6.34080688e-02,
                1.05589989e-01,1.68657446e-01,1.21704836e+03])
        self.check_freqs(expected_freqs, nma, -1, check_zeros=True)

        #  --- atoms of subsystem are not collinear
        subs = load_indices(
            pkg_resources.resource_filename(__name__, "../data/test/an/fixed.03.txt"))  # atom 1 to atom 7
        nma = NMA(molecule, VSANoMass(subs))
        self.assert_(len(nma.zeros)==6)
        expected_freqs = np.array([ # taken from previous working python version
                        -0.0299898080453, -0.0124485604422, -0.000583010131998, -0.000190084696013,
                        0.000469660947606, 0.0233987667917,    # These are the zeros.
                        388.990852247, 1055.12229682, 1251.18647899, 1352.32135672,
                        2864.80440032, 3680.49886454 ] )
        self.check_freqs(expected_freqs, nma, 4, check_zeros=True)

    def test_mbh(self):
        molecule = load_molecule_charmm(
            pkg_resources.resource_filename(__name__, "../data/test/an/ethanol.cor"),
            pkg_resources.resource_filename(__name__, "../data/test/an/ethanol.hess.full"))
        blocks = load_indices(
            pkg_resources.resource_filename(__name__, "../data/test/an/fixed.07.txt"),
            groups=True)
        nma = NMA(molecule, MBH(blocks))
        self.check_ortho(nma.modes)
        with tmpdir(__name__, 'test_mbh') as dn:
            dump_modes_molden(os.path.join(dn, "ethanol.mbh.molden.log"), nma)
        self.check_ortho(nma.modes)   # write_molden should not have changed this

    def test_mbh_ethane(self):
        molecule = load_molecule_g03fchk(
            pkg_resources.resource_filename(__name__, "../data/test/ethane/gaussian.fchk"))
        blocks = [[1, 0, 2, 6, 7], [1, 0, 3, 4, 5 ]]
        nma = NMA(molecule, MBH(blocks))
        self.assertEqual(len(nma.freqs), 7)
        non_zero = [i for i in xrange(7) if i not in nma.zeros][0]
        self.assertAlmostEqual(nma.freqs[non_zero]/lightspeed*centimeter, 314, 0)

    def test_mbhconstrainext(self):
        # load the plain Hessian
        molecule = load_molecule_g03fchk(
            pkg_resources.resource_filename(__name__, "../data/test/sterck/aa.fchk"))
        blocks = [[3,2,6],[6,7,8]]
        nma1 = NMA(molecule)
        nma2 = NMA(molecule, ConstrainExt(gradient_threshold=1e-2))
        nma3 = NMA(molecule, MBH(blocks))
        nma4 = NMA(molecule, MBHConstrainExt(blocks))
        # construct a molecule with the projected Hessian and projected gradient
        molecule = molecule.constrain_ext()
        nma5 = NMA(molecule)
        nma6 = NMA(molecule, ConstrainExt(gradient_threshold=1e-2))
        nma7 = NMA(molecule, MBH(blocks))
        nma8 = NMA(molecule, MBHConstrainExt(blocks))
        #for i in range(10):
        #    print "%10.6f %10.6f %10.6f %10.6f %10.6f" % (
        #      nma1.freqs[i]/lightspeed*centimeter,  #full
        #      nma3.freqs[i]/lightspeed*centimeter, nma4.freqs[i]/lightspeed*centimeter,
        #      nma7.freqs[i]/lightspeed*centimeter, nma8.freqs[i]/lightspeed*centimeter,)
        self.check_freqs(nma3.freqs/lightspeed*centimeter, nma4, 1, check_zeros=True)
        # the lowest freqs are expected to be different
        self.check_freqs(nma3.freqs[6:]/lightspeed*centimeter, nma7, 0, check_zeros=False)
        self.check_freqs(nma7.freqs/lightspeed*centimeter, nma8, 1, check_zeros=True)

    def test_mbh_raise_ext(self):
        # load the plain Hessian
        molecule = load_molecule_charmm(
            pkg_resources.resource_filename(__name__, "../data/test/an/ethanol.cor"),
            pkg_resources.resource_filename(__name__, "../data/test/an/ethanol.hess.full"))
        blocks = [[3,2,6],[6,7,8]]
        nma1 = NMA(molecule)
        nma2 = NMA(molecule, MBH(blocks))
        nma3 = NMA(molecule, MBHConstrainExt(blocks))
        # raise the eigenvalue of the six global translationas/rotations
        molecule = molecule.raise_ext()
        nma5 = NMA(molecule)
        nma7 = NMA(molecule, MBH(blocks))
        #print "First freqs"
        #for i in range(10):
        #    print "%10.6f %10.6f %10.6f %10.6f %10.6f" % (
        #      nma1.freqs[i]/lightspeed*centimeter, nma5.freqs[i]/lightspeed*centimeter, #full
        #      nma2.freqs[i]/lightspeed*centimeter,  # mbh
        #      nma3.freqs[i]/lightspeed*centimeter, nma7.freqs[i]/lightspeed*centimeter,)
        #print "Last freqs"
        #for i in range(10):
        #    print "%10.6f %10.6f %10.6f %10.6f %10.6f" % (
        #      nma1.freqs[-i]/lightspeed*centimeter, nma5.freqs[-i]/lightspeed*centimeter, #full
        #      nma2.freqs[-i]/lightspeed*centimeter, # mbh
        #      nma3.freqs[-i]/lightspeed*centimeter, nma7.freqs[-i]/lightspeed*centimeter,)
        # compare full, not the zeros
        for i in range(len(nma1.freqs)-6):
            self.assertAlmostEqual( nma1.freqs[i+6],nma5.freqs[i],7)
        # compare mbh, not the zeros
        for i in range(len(nma2.freqs)-6):
            self.assertAlmostEqual( nma2.freqs[i+6],nma3.freqs[i+6],7)
            self.assertAlmostEqual( nma2.freqs[i+6],nma7.freqs[i],7)


    def test_phva_mbh(self):
        molecule = load_molecule_charmm(
            pkg_resources.resource_filename(__name__, "../data/test/an/ethanol.cor"),
            pkg_resources.resource_filename(__name__, "../data/test/an/ethanol.hess.full"))
        blocks = load_indices(
            pkg_resources.resource_filename(__name__, "../data/test/an/fixed.07.txt"),
            groups=True)
        fixed = load_indices(
            pkg_resources.resource_filename(__name__, "../data/test/an/fixed.06.txt"))
        nma = NMA(molecule, PHVA_MBH(fixed,blocks))
        expected_freqs = np.array([214.28936269,   596.97481532,   663.5290044,    787.84964637,
                 859.68400023, 1096.26899018,  1160.4094257, 1224.87768013,  1310.51036299,
                 1411.62142706, 1509.02056509,  2218.9706113,  2762.11931539,  2942.89826034,
                 3173.07796648])
        self.check_freqs(expected_freqs, nma, 6, check_zeros=True)
        self.check_ortho(nma.modes)
        self.assertEqual(nma.modes.shape[0],27)
        self.assertEqual(nma.modes.shape[1],15)

    def test_constrain1(self):
        molecule = load_molecule_charmm(
            pkg_resources.resource_filename(__name__, "../data/test/an/ethanol.cor"),
            pkg_resources.resource_filename(__name__, "../data/test/an/ethanol.hess.full"))
        fixed = [[1,2]]
        nma = NMA(molecule, Constrain(fixed))
        #print nma.freqs/lightspeed*centimeter
        self.check_ortho(nma.modes)
        self.assertEqual(nma.modes.shape[0],27)
        self.assertEqual(nma.modes.shape[1],26)
        with tmpdir(__name__, 'test_constrain1') as dn:
            dump_modes_molden(os.path.join(dn, "ethanol.constr.molden.log"), nma)

    def test_constrain2(self):
        molecule = load_molecule_charmm(
            pkg_resources.resource_filename(__name__, "../data/test/an/ethanol.cor"),
            pkg_resources.resource_filename(__name__, "../data/test/an/ethanol.hess.full"))
        fixed = [[1,2], [0,4], [0,5],[2,8],[3,4],[4,5],[5,6],[6,7],[7,8],[7,1],[7,3]]
        nma = NMA(molecule, Constrain(fixed))
        #print nma.freqs/lightspeed*centimeter
        self.check_ortho(nma.modes)
        self.assertEqual(nma.modes.shape[0],27)
        self.assertEqual(nma.modes.shape[1],16)
        with tmpdir(__name__, 'test_constrain2') as dn:
            dump_modes_molden(os.path.join(dn, "ethanol.constr.molden.log"), nma)

    def test_sandra(self):
        cases = [
            (pkg_resources.resource_filename(__name__, "../data/test/sandra/F_freq.fchk"),
             []),
            (pkg_resources.resource_filename(__name__, "../data/test/sandra/HF_freq.fchk"),
             [4472.7168]),
            (pkg_resources.resource_filename(__name__, "../data/test/sandra/ts_FHF_optfreq.fchk"),
             [-4.1374, 2.0214, 13.3582, 106.6359]),
        ]
        for fn_fchk, expected_freqs in cases:
            molecule = load_molecule_g03fchk(fn_fchk)
            nma = NMA(molecule, ConstrainExt(gradient_threshold=1e-3))
            self.check_ortho(nma.modes)
            self.check_freqs(expected_freqs, nma, 2, check_zeros=False)
            nma = NMA(molecule)
