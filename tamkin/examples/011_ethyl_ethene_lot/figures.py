#!/usr/bin/env python
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


from lot_basis import lots_list
import html
from kin import *

from molmod.units import kjmol
from molmod.constants import boltzmann

import matplotlib.pyplot as pt, numpy


# ugly hack
import matplotlib.mathtext
matplotlib.mathtext.Parser._function_names.add("theory")


def get_av(fn_template):
    try:
        gauche_k = load_summary(fn_template % "gauche")[0]
        trans_k = load_summary(fn_template % "trans")[0]
    except IOError:
        return None

    gauche_lnk = numpy.log(gauche_k)
    trans_lnk = numpy.log(trans_k)

    lnav = ((gauche_lnk - experimental_lnk).mean() + (trans_lnk - experimental_lnk).mean())/2
    return numpy.exp(lnav)


def get_rmsd(fn_template):
    try:
        gauche_k = load_summary(fn_template % "gauche")[0]
        trans_k = load_summary(fn_template % "trans")[0]
    except IOError:
        return None

    gauche_lnk = numpy.log(gauche_k)
    trans_lnk = numpy.log(trans_k)

    lnav = numpy.sqrt(((gauche_lnk - experimental_lnk)**2 + (trans_lnk - experimental_lnk)**2).sum()/8)
    return numpy.exp(lnav)


#get_data = get_av
get_data = get_rmsd


def do_plot(data, legend, fn_png):
    pt.clf()
    pt.figure(0, (6,8))
    #pt.axes([left, bottom, width, height])
    # in fractional coordinates
    pt.axes([0.20, 0.070, 0.76, 0.84])
    handles = []
    labels = []
    for key, symbol, label in legend:
        x = []
        y = []
        for index, lot in enumerate(lots_list):
            av = data.get((key, lot.key))
            if av is not None:
                x.append(av)
                y.append(-index)
        handle = pt.plot(x, y, symbol)
        handles.append(handle)
        labels.append(label)
    pos = -numpy.arange(len(lots_list))
    ylabels = [lot.key for lot in lots_list]
    pt.yticks(pos, ylabels, size="small")
    pt.xticks(size="small")
    #pt.xlabel(r"Geometric mean of $k_{\theory}/k_{\exp}$", size="small")
    pt.xlabel(r"$\exp(RMSD(\ln(k_{\theory}) - \ln($k_{\exp})))$", size="small")
    pt.semilogx()
    pt.axvline(0.1, color="k", zorder=-5)
    pt.axvline(1.0, color="k", zorder=-5)
    pt.axvline(10.0, color="k", zorder=-5)
    for y in 0, -5, -10, -15, -20, -25, -30, -35, -40:
        pt.axhline(y, color="k", alpha=0.2, lw=10, zorder=-5)
    if len(legend) > 1:
        legend = pt.figlegend(handles, labels, (0.22,0.915), numpoints=1, handlelength=0, prop={"size":"small"})
        frame = legend.get_frame()
        frame.set_lw(0.0)
    #pt.xlim(1e-8, 1e5)
    pt.xlim(1e0, 1e5)
    pt.ylim(pos[-1]-0.5, 0.5)
    pt.savefig(fn_png)


def figure1():
    #  - comparison of basis-set/geometry for each lot
    #  - remains fixed: HO, BSS
    #  - geometric average (over gauche, trans and temps) of ratio between theory and experiment


    dir_templates = [
        ("low", "%s__6-31gd"),
        ("high", "%s__6-311+g3df2p"),
        ("mixed", "GEO__b3lyp__6-31gd__ENERGY__%s__6-311+g3df2p"),
    ]

    data = {}
    for key, dir_template in dir_templates:
        for lot in lots_list:
            dirname = dir_template % lot.key
            data[(key, lot.key)] = get_data(dirname  + "/ho_bss_summary_%s.txt")

    legend = [
        ("low", "ro", "Consistent, 6-31G(d)"),
        ("high", "gD", "Consistent, 6-311G(3df,2p)"),
        ("mixed", "bs", "2-Component"),
    ]
    do_plot(data, legend, "figure1.png")


def figure2():
    #  - comparison of harmonic oscillator versus inernal rotor for each lot
    #  - remains fixed: BSS, B3LYP geo
    #  - geometric average (over gauche, trans and temps) of ratio between theory and experiment

    data = {}
    for key in "ho", "ir":
        for lot in lots_list:
            dirname = "GEO__b3lyp__6-31gd__ENERGY__%s__6-311+g3df2p" % lot.key
            data[(key, lot.key)] = get_data("%s/%s_bss_summary_%%s.txt" % (dirname, key))

    legend = [
        ("ho", "bs", "2-Component, Harmonic oscillator"),
        ("ir", "yH", "2-Component, 1D Internal rotors"),
    ]
    do_plot(data, legend, "figure2.png")


def figure3():
    #  - comparison of bss and cps for each lot
    #  - remains fixed: HO, B3LYP geo
    #  - geometric average (over gauche, trans and temps) of ratio between theory and experiment

    data = {}
    for key in "bss", "cps":
        for lot in lots_list:
            dirname = "GEO__b3lyp__6-31gd__ENERGY__%s__6-311+g3df2p" % lot.key
            data[(key, lot.key)] = get_data("%s/ho_%s_summary_%%s.txt" % (dirname, key))

    legend = [
        ("bss", "bs", "2-Component, Without correction"),
        ("cps", "m^", "2-Component, Counterpoise corrected barrier"),
    ]
    do_plot(data, legend, "figure3.png")


def figure4():
    #  - final ratios
    #  - remains fixed: IR, CPS, B3LYP geo
    #  - geometric average (over gauche, trans and temps) of ratio between theory and experiment

    data = {}
    for lot in lots_list:
        dirname = "GEO__b3lyp__6-31gd__ENERGY__%s__6-311+g3df2p" % lot.key
        data[("foo", lot.key)] = get_data("%s/ho_bss_summary_%%s.txt" % dirname)
        data[("bar", lot.key)] = get_data("%s/ir_cps_summary_%%s.txt" % dirname)

    legend = [
        ("foo", "bs", "2-Component, Harmonic oscillator, without correction"),
        ("bar", "c*", "2-Component, 1D internal rotors, counterpoise corrected barrier"),
    ]
    do_plot(data, legend, "figure4.png")



if __name__ == "__main__":
    figure1()
    figure2()
    figure3()
    figure4()
