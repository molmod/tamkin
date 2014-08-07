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
'''Functions to load the dispersion correction from file'''


__all__ = ['load_dftd3', 'load_dftd_orca']


def load_dftd3(fn):
    '''Load the dispersion correction (in a.u.) from a DFT-D3 output file

       Arguments:
        | ``fn`` -- The DFT-D3 output file.
    '''
    with open(fn) as f:
        for line in f:
            if line.startswith(' Edisp /kcal,au:'):
                return float(line.split()[-1])


def load_dftd_orca(fn):
    '''Load the dispersion correction (in a.u.) from on Orca output file

       Arguments:
        | ``fn`` -- The DFT-D3 output file.
    '''
    with open(fn) as f:
        for line in f:
            if line.startswith("Van der Waals correction ="):
                return float(line.split()[5])
