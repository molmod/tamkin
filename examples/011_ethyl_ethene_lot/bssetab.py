#!/usr/bin/env python
# -*- coding:utf8 -*-


from lot_basis import lots_list
import html

def load_summary(fn):
    f = file(fn)
    result = tuple(float(word) for word in f.next().split())
    f.close()
    return result[:-4], result[-4:-2], result[-2:]


def ts_bsse(fn_template, title, rows):
    if len(rows) == 0:
        rows.append([u"<th>Approach→<br />Functional↓</th>"])
        rows.append(["<th></th>"])
    rows[0].append("<th colspan=3>%s</th>" % title.replace(", ", "<br />"))
    rows[1].append(u"<td>ΔE</td><td>ΔE'</td><td>ΔE'-ΔE</td>")
    counter = 2
    for lot in lots_list:
        if lot.spin == "ROS":
            lot_label = "ro" + lot.label
        else:
            lot_label = lot.label
        if len(rows) <= counter:
            rows.append(["<th>%s</th>" % lot_label])
        try:
            (A, Ea), (dE0, dE) = load_summary(fn_template % (lot_label, "bss"))[1:]
            (Ap, Eap), (dE0p, dEp) = load_summary(fn_template % (lot_label, "cps"))[1:]
            rows[counter].append("<td>%.0f</td>" % dE)
            rows[counter].append("<td>%.0f</td>" % dEp)
            rows[counter].append("<td style='background-color:#DDD;'>%.0f</td>" % (dEp-dE))
            #rows[counter].append("<td>%.0f</td>" % Ea)
            #rows[counter].append("<td>%.0f</td>" % Eap)
            #rows[counter].append("<td style='background-color:#DDD;'>%.0f</td>" % (Eap-Ea))
        except IOError, StopIteration:
            #rows[counter].append("<td>&nbsp</td><td>&nbsp</td><td>&nbsp</td>")
            rows[counter].append("<td>&nbsp</td><td>&nbsp</td><td>&nbsp</td>")
        counter += 1

f = file("bssetab.html", "w")
print >> f, html.header % "BSSE Overview"

rows = []
ts_bsse("%s__6-31gd/ho_%s_summary_gauche.txt", "Gauche, HO, Consistent, 6-31G(d)", rows)
ts_bsse("%s__6-31gd/ho_%s_summary_trans.txt", "Trans, HO, Consistent, 6-31G(d)", rows)
ts_bsse("%s__6-311+g3df2p/ho_%s_summary_gauche.txt", "Gauche, HO, Consistent, 6-311+G(3df,2p)", rows)
ts_bsse("%s__6-311+g3df2p/ho_%s_summary_trans.txt", "Trans, HO, Consistent, 6-311+G(3df,2p)", rows)
ts_bsse("GEO__b3lyp__6-31gd__ENERGY__%s__6-311+g3df2p/ho_%s_summary_gauche.txt", "Gauche, HO, GEO=B3LYP/6-31G(d), 6-311+G(3df,2p)", rows)
ts_bsse("GEO__b3lyp__6-31gd__ENERGY__%s__6-311+g3df2p/ho_%s_summary_trans.txt", "Trans, HO, GEO=B3LYP/6-31G(d), 6-311+G(3df,2p)", rows)

print >> f, "<p>BSSE corrections on the transition state.</p>"
html.print_table(f, rows)

print >> f, html.footer
