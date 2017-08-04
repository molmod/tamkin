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


from molmod.test.test_examples import check_example


def test_example_001():
    check_example(__name__, "001_ethane", "thermo.py", ['gaussian.fchk'])

def test_example_002():
    check_example(__name__, "002_linear_co2", "thermo.py", ['gaussian.fchk'])

def test_example_003():
    check_example(__name__, "003_pentane", "thermo.py", ['gaussian.fchk'])

def test_example_005():
    check_example(__name__, "005_acrylamide_reaction", "reaction.py",
                  ['aa.fchk', 'aarad.fchk', 'paats.fchk'])

def test_example_006():
    check_example(__name__, "006_5T_ethyl_ethene_addition", "reaction.py",
                  ['react.fchk', 'ts.fchk'])

def test_example_007():
    check_example(__name__, "007_mfi_propene_reaction", "reaction.py",
                  ['Zp_p_react.28aug.com', 'Zp_p_TS.28aug.com', 'Zp_p_react.28aug.fchk',
                   'Zp_p_react.14mei.fchk', 'Zp_p_TS.28aug.fchk',
                   '5Tp_p_TS.oniom21apr_HF.fchk'])

def test_example_008():
    check_example(__name__, "008_ethane_rotor", "thermo.py",
                  ['freq/gaussian.fchk', 'scan/gaussian.log'])

def test_example_009():
    check_example(__name__, "009_ethyl_ethene", "reaction.py",
                  ['ethyl/freq/gaussian.fchk', 'ethene/freq/gaussian.fchk',
                   'ts_ad1/freq_gauche/gaussian.fchk', 'ts_ad1/freq_trans/gaussian.fchk',
                   'ethyl/scan/gaussian.log', 'ts_ad1/scan1/gaussian.log',
                   'ts_ad1/scan2/gaussian.log'])

def test_example_012():
    check_example(__name__, "012_ethyl_ethene_scaling", "reaction.py",
                  ['ethyl/freq/gaussian.fchk', 'ethene/freq/gaussian.fchk',
                   'ts_ad1/freq_gauche/gaussian.fchk'])

def test_example_013():
    check_example(__name__, "013_butane", "thermo.py",
                  ['freq/gaussian.fchk', 'scan/gaussian.log'])

def test_example_014():
    check_example(__name__, "014_pentane_mbh", "thermo.py", ['gaussian.fchk'])

def test_example_015():
    check_example(__name__, "015_kie", "reaction.py", ['reactant.fchk', 'trans.fchk'])

def test_example_016():
    check_example(__name__, "016_modes", "modes.py", ['PUNCH'])

def test_example_017():
    check_example(__name__, "017_activationkineticmodel", "reaction.py",
                  ['water/gaussian.fchk', 'VO_AA_OH_H2O/gaussian.fchk',
                   'VO_AA_OOtBu/gaussian.fchk', 'TBHP/gaussian.fchk', 'TS/TS.fchk',
                   'cyclohexene/gaussian.fchk'])

def test_example_018():
    check_example(__name__, "018_physisorption", "adsorption.py",
                  ['m062x/argon/gaussian.fchk', 'm062x/benzene/freq/gaussian.fchk',
                   'm062x/complex/freq/gaussian.fchk'])

def test_example_019():
    check_example(__name__, "019_ethyl_ethene_simple", "kinetic.py",
                  ['ethyl.fchk', 'ethene.fchk', 'ts_trans.fchk'])

def test_example_020():
    check_example(__name__, "020_butane_conformers", "equilibrium.py",
                  ['trans.fchk', 'gauche.fchk'])

def test_example_021():
    check_example(__name__, "021_water_formation", "formation.py",
                  ['oxygen.fchk', 'hydrogen.fchk', 'water.fchk'])
