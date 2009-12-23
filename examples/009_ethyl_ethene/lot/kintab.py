#!/usr/bin/env python
# -*- coding:utf8 -*-


from lot_basis import lots_list

from molmod.units import kjmol
from molmod.constants import boltzmann

import pylab, numpy

from matplotlib.backend_bases import GraphicsContextBase
GraphicsContextBase.dashd["dashed"] = (0, (6.0, 3.0))
GraphicsContextBase.dashd["dashdot"] = (0, (4.0, 2.0, 1.0, 2.0))
GraphicsContextBase.dashd["dotted"] = (0, (1.5, 1.5))

def load_summary(fn):
    f = file(fn)
    result = tuple(float(word) for word in f.next().split())
    f.close()
    return result

temps = numpy.array([300, 400, 500, 600], float)
invtemps = 1/temps
experimental_lnk = numpy.array([-0.6103, 2.3416, 4.2361, 5.5809], float)
experimental_k = numpy.exp(experimental_lnk)

# fit experimental A and Ea and also sensitivity to errors

design_matrix = numpy.array([numpy.ones(4), -invtemps/boltzmann*kjmol]).transpose()
expected_values = experimental_lnk

A = numpy.dot(design_matrix.transpose(), design_matrix)
B = numpy.dot(design_matrix.transpose(), expected_values)

experimental_parameters = numpy.linalg.solve(A, B)
experimental_A = numpy.exp(experimental_parameters[0])
experimental_Ea = experimental_parameters[1]
covariance_lnk = numpy.ones((4,4),float)*numpy.log(10)**2 + numpy.identity(4,float)*numpy.log(2)**2
#covariance_lnk = numpy.identity(4,float)*numpy.log(10)**2
sensitivity = numpy.linalg.solve(A, design_matrix.transpose())
covariance_parameters = numpy.dot(numpy.dot(sensitivity, covariance_lnk), sensitivity.transpose())
print covariance_parameters
#sys.exit()
del design_matrix
del expected_values
del A
del B

# done fitting

def get_error_color(ratio):
    x = abs(ratio)
    if x < 1:
        rgb = (x,1.0,0.0)
    elif x < 2:
        rgb = (1.0,2.0-x,0.0)
    elif x < 3:
        rgb = (3.0-x,0.0,0.0)
    else:
        rgb = (0.0,0.0,0.0)
    result = "#%s" % "".join(hex(int(c*255))[2:].rjust(2,"0") for c in rgb)
    return result

def overview(template, title, fn_img, rows):
    if len(rows) == 0:
        rows.append([u"<th>Approach→<br />Functional↓</th>"])
        rows.append(["<th></th>"])
    rows[0].append("<th colspan=4>%s</th>" % title.replace(", ", "<br />"))
    for temp in temps:
        rows[1].append("<td style='font-size:small'>%.0fK</td>" % temp)
    lines = []
    labels = []
    pylab.clf()
    line = pylab.plot(invtemps, experimental_k, color="k", linestyle="-", lw=4)
    lines.append(line)
    labels.append("experiment")
    counter = 2
    for lot in lots_list:
        if lot.spin == "ROS":
            label = "ro" + lot.label
        else:
            label = lot.label
        if len(rows) <= counter:
            rows.append(["<th>%s</th>" % label])
        try:
            values = load_summary(template % label)
            line = pylab.plot(invtemps, values[:-2], color=lot.color, linestyle=lot.linestyle, lw=2)
            lines.append(line)
            labels.append(label)
            for j in xrange(4):
                ln10ratio = numpy.log10(values[j]/experimental_k[j])
                color = get_error_color(ln10ratio)
                rows[counter].append("<td style='background-color:%s'>%.0f</td>" % (color, ln10ratio*10))
            counter += 1
            continue
        except IOError, StopIteration:
            pass
        rows[counter].append("<td>&nbsp</td><td>&nbsp</td><td>&nbsp</td><td>&nbsp</td>")
        counter += 1
    pylab.semilogy()
    pylab.fill(
        numpy.concatenate((invtemps, invtemps[::-1])),
        numpy.concatenate((experimental_k*10, experimental_k[::-1]/10)),
        "k", alpha=0.2, lw=0,
    )
    pylab.xticks(
        1/numpy.array([300, 350, 400, 450, 500, 550, 600], float),
        ["300", "350", "400", "450", "500", "550", "600"],
    )
    pylab.title(title)
    pylab.xlabel("T [K]")
    pylab.ylabel("k [(m**3/mol)/s]")
    pylab.ylim(1e-8,1e7)
    legend = pylab.figlegend(
        lines, labels, (0.07,0.06), ncol=3, prop={"size":11},
        handlelength=3, labelspacing=0.1, columnspacing=1
    )
    #legend.get_frame().set_linewidth(0)
    legend.get_frame().set_fill(True)
    legend.get_frame().set_alpha(0.5)
    pylab.savefig(fn_img % "rates")

    pylab.clf()
    lines = []
    labels = []
    line = pylab.plot([experimental_Ea], [experimental_A], color="k", marker="o", ms=11, mew=2, lw=0, ls=" ")
    lines.append(line)
    labels.append("experiment")
    for lot in lots_list:
        if lot.spin == "ROS":
            label = "ro" + lot.label
        else:
            label = lot.label
        try:
            A, Ea = load_summary(template % label)[-2:]
            marker = {"-": "o", "--": "s", ":": "v", "-.": "h"}[lot.linestyle]
            line = pylab.plot([Ea], [A], color=lot.color, marker=marker, ms=11, mew=2, lw=0, ls=" ")
            lines.append(line)
            labels.append(label)
            continue
        except IOError, StopIteration:
            pass

    pylab.title(title)
    pylab.xlabel("Activation energy [kJ/mol]")
    pylab.ylabel("Pre-exponential factor [(m**3/mol)/s]")
    pylab.semilogy()
    # error margin around experimental data point
    x = []
    y = []
    evals, evecs = numpy.linalg.eigh(covariance_parameters)
    angles = numpy.arange(0.0,360.5,1.0)/180*numpy.pi
    data = numpy.outer(evecs[:,0],numpy.cos(angles))*numpy.sqrt(evals[0]) + \
           numpy.outer(evecs[:,1],numpy.sin(angles))*numpy.sqrt(evals[1])
    pylab.fill(
        experimental_parameters[1] + data[1],
        numpy.exp(experimental_parameters[0] + data[0]),
        color="k", alpha=0.2, lw=0
    )
    # end error margin
    legend = pylab.legend(
        lines, labels, loc=4, ncol=4, prop={"size":11},
        handlelength=1, labelspacing=0.2, columnspacing=2,
        numpoints=1,
    )
    pylab.xlim(0,90)
    pylab.ylim(1e3,1e7)
    #legend.get_frame().set_linewidth(0)
    legend.get_frame().set_fill(True)
    legend.get_frame().set_alpha(0.5)
    pylab.savefig(fn_img % "params")




f = file("kintab.html", "w")
print >> f, "<?xml version='1.0' encoding='UTF-8'?>"
print >> f, "<html><head><title>LOT Overview</title>"
print >> f, "<style type='text/css'>"
print >> f, "body { font-family: sans; }"
print >> f, "td { text-align: right; width: 50px; padding: 0px 10px; }"
print >> f, "table, td, th, tr { border: solid 1px black; }"
print >> f, "table { table-layout: fixed; border-collapse:collapse; }"
print >> f, "</style>"
print >> f, "</head><body>"

for do_rotor in False, True:
    ir_str = {True: "ir", False: "ho"}[do_rotor]
    for do_counterpoise in False, True:
        cp_str = {True: "cps", False: "bss"}[do_counterpoise]

        rows = []

        for ts_conformer in "Gauche", "Trans":
            overview(
                "%%s__6-31gd/%s_%s_summary_%s.txt" % (ir_str, cp_str, ts_conformer.lower()),
                "%s, %s, %s, Consistent, 6-31G(d)" % (ir_str.upper(), cp_str.upper(), ts_conformer),
                "kin_%s_%s_%s_consistent_6-31gd_%%s.png" % (ir_str, cp_str, ts_conformer.lower()),
                rows,
            )
        for ts_conformer in "Gauche", "Trans":
            overview(
                "GEO__b3lyp__6-31gd__ENERGY__%%s__6-311+g3df2p/%s_%s_summary_%s.txt" % (ir_str, cp_str, ts_conformer.lower()),
                "%s, %s, %s, GEO=B3LYP/6-31G(d), 6-311+G(3df,2p)" % (ir_str.upper(), cp_str.upper(), ts_conformer),
                "kin_%s_%s_%s_geo_6-311+g3df2p_%%s.png" % (ir_str, cp_str, ts_conformer.lower()),
                rows,
            )
        for ts_conformer in "Gauche", "Trans":
            overview(
                "%%s__6-311+g3df2p/%s_%s_summary_%s.txt" % (ir_str, cp_str, ts_conformer.lower()),
                "%s, %s, %s, Consistent, 6-311+G(3df,2p)" % (ir_str.upper(), cp_str.upper(), ts_conformer),
                "kin_%s_%s_%s_consistent_6-311+g3df2p_%%s.png" % (ir_str, cp_str, ts_conformer.lower()),
                rows,
            )

        print >> f, "<p>10*Log10 of ratio between theoretical and experimental rate constant.</p>"
        print >> f, "<table style='border-color:black'>"
        for row in rows:
            print >> f, "<tr>%s</tr>" % ("".join(row).encode('UTF-8'))
        print >> f, "</table>"

print >> f, "</body></html>"

