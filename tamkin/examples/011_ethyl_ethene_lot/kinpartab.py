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

import numpy as np


def overview(template, title, fn_img, rows):
    if len(rows) == 0:
        rows.append([u"<th>Approach→<br />Functional↓</th>"])
        rows.append(["<th></th>"])
    rows[0].append("<th colspan=2>%s</th>" % title.replace(", ", "<br />"))
    rows[1].append("<td style='font-size:small'>A</td><td style='font-size:small'>E<sub>a</sub></td>")
    counter = 2
    for lot in lots_list:
        if lot.spin == "ROS":
            lot_label = "ro" + lot.label
        else:
            lot_label = lot.label
        if len(rows) <= counter:
            rows.append(["<th>%s</th>" % lot_label])
        try:
            A, Ea = load_summary(template % lot_label)[1]
            color = get_error_color(abs(np.log10(A/experimental_A)))
            rows[counter].append("<td style='background-color:%s'>%.1e</td>" % (color, A))
            color = get_error_color(abs(Ea - experimental_Ea)/4)
            rows[counter].append("<td style='background-color:%s'>%.1f</td>" % (color, Ea))
        except IOError, StopIteration:
            rows[counter].append("<td>&nbsp</td><td>&nbsp</td>")
        counter += 1




f = file("kinpartab.html", "w")
print >> f, html.header % "KIN par Overview"

for do_rotor in False, True:
    ir_str = {True: "ir", False: "ho"}[do_rotor]
    ir_info = {
        True: "internal rotor",
        False: "harmonic oscillator",
    }[do_rotor]
    for do_counterpoise in False, True:
        cp_str = {True: "cps", False: "bss"}[do_counterpoise]
        cp_info = {
            True: "with counterpoise correction",
            False: "without counterpoise correction",
        }[do_counterpoise]

        rows = []

        for ts_conformer in "Gauche", "Trans":
            overview(
                "%%s__6-31gd/%s_%s_summary_%s.txt" % (ir_str, cp_str, ts_conformer.lower()),
                "%s, %s, %s, Consistent, 6-31G(d)" % (ir_str.upper(), cp_str.upper(), ts_conformer),
                "kin_%s_%s_%s_consistent_6-31gd_%%s.pdf" % (ir_str, cp_str, ts_conformer.lower()),
                rows,
            )
        for ts_conformer in "Gauche", "Trans":
            overview(
                "%%s__6-311+g3df2p/%s_%s_summary_%s.txt" % (ir_str, cp_str, ts_conformer.lower()),
                "%s, %s, %s, Consistent, 6-311+G(3df,2p)" % (ir_str.upper(), cp_str.upper(), ts_conformer),
                "kin_%s_%s_%s_consistent_6-311+g3df2p_%%s.pdf" % (ir_str, cp_str, ts_conformer.lower()),
                rows,
            )
        for ts_conformer in "Gauche", "Trans":
            overview(
                "GEO__b3lyp__6-31gd__ENERGY__%%s__6-311+g3df2p/%s_%s_summary_%s.txt" % (ir_str, cp_str, ts_conformer.lower()),
                "%s, %s, %s, GEO=B3LYP/6-31G(d), 6-311+G(3df,2p)" % (ir_str.upper(), cp_str.upper(), ts_conformer),
                "kin_%s_%s_%s_geo_6-311+g3df2p_%%s.pdf" % (ir_str, cp_str, ts_conformer.lower()),
                rows,
            )

        print >> f, "<p>Kinetic parameters (%s, %s)</p>" % (ir_info, cp_info)
        html.print_table(f, rows)

print >> f, html.footer
