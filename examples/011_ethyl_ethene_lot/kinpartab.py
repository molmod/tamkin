#!/usr/bin/env python
# -*- coding:utf8 -*-


from lot_basis import lots_list
import html
from kin import *

from molmod.units import kjmol
from molmod.constants import boltzmann

import numpy


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
            color = get_error_color(abs(numpy.log10(A/experimental_A)))
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
