#!/usr/bin/env python
# -*- coding:utf8 -*-


from lot_basis import lots_list

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
experimental_k = numpy.exp(numpy.array([-0.6103, 2.3416, 4.2361, 5.5809], float))

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
    pylab.clf()
    lines = []
    labels = []
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
    legend.get_frame().set_linewidth(0)
    legend.get_frame().set_fill(False)
    pylab.savefig(fn_img)

    #print >> f, "<p><img src='%s'></p>" % png


f = file("summary.html", "w")
print >> f, "<?xml version='1.0' encoding='UTF-8'?>"
print >> f, "<html><head><title>LOT Overview</title>"
print >> f, "<style type='text/css'>"
print >> f, "body { font-family: sans; }"
print >> f, "td { text-align: right; width: 50px; padding: 0px 10px; }"
print >> f, "table, td, th, tr { border: solid 1px black; }"
print >> f, "table { table-layout: fixed; border-collapse:collapse; }"
print >> f, "</style>"
print >> f, "</head><body>"


rows = []

overview(
    "%s__6-31gd/ho_summary_gauche.txt",
    "Gauche, HO, Consistent, 6-31G(d)",
    "summary_gauche_ho_consistent_6-31gd.png",
    rows,
)
overview(
    "%s__6-31gd/ho_summary_trans.txt",
    "Trans, HO, Consistent, 6-31G(d)",
    "summary_trans_ho_consistent_6-31gd.png",
    rows,
)
overview(
    "GEO__b3lyp__6-31gd__ENERGY__%s__6-311+g3df2p/ho_summary_gauche.txt",
    "Gauche, HO, GEO=B3LYP/6-31G(d), 6-311+G(3df,2p)",
    "summary_gauche_ho_geo_6-311+g3df2p.png",
    rows,
)
overview(
    "GEO__b3lyp__6-31gd__ENERGY__%s__6-311+g3df2p/ho_summary_trans.txt",
    "Trans, HO, GEO=B3LYP/6-31G(d), 6-311+G(3df,2p)",
    "summary_trans_ho_geo_6-311+g3df2p.png",
    rows,
)
overview(
    "%s__6-311+g3df2p/ho_summary_gauche.txt",
    "Gauche, HO, Consistent, 6-311+G(3df,2p)",
    "summary_gauche_ho_consistent_6-311+g3df2p.png",
    rows,
)
overview(
    "%s__6-311+g3df2p/ho_summary_trans.txt",
    "Trans, HO, Consistent, 6-311+G(3df,2p)",
    "summary_trans_ho_consistent_6-311+g3df2p.png",
    rows,
)

print >> f, "10*Log10 of ratio between theoretical and experimental rate constant."
print >> f, "<table style='border-color:black'>"
for row in rows:
    print >> f, "<tr>%s</tr>" % ("".join(row).encode('UTF-8'))
print >> f, "</table>"

rows = []

overview(
    "%s__6-31gd/ir_summary_gauche.txt",
    "Gauche, IR, Consistent, 6-31G(d)",
    "summary_gauche_ir_consistent_6-31gd.png",
    rows,
)
overview(
    "%s__6-31gd/ir_summary_trans.txt",
    "Trans, IR, Consistent, 6-31G(d)",
    "summary_trans_ir_consistent_6-31gd.png",
    rows,
)
overview(
    "GEO__b3lyp__6-31gd__ENERGY__%s__6-311+g3df2p/ir_summary_gauche.txt",
    "Gauche, IR, GEO=B3LYP/6-31G(d), 6-311+G(3df,2p)",
    "summary_gauche_ir_geo_6-311+g3df2p.png",
    rows,
)
overview(
    "GEO__b3lyp__6-31gd__ENERGY__%s__6-311+g3df2p/ir_summary_trans.txt",
    "Trans, IR, GEO=B3LYP/6-31G(d), 6-311+G(3df,2p)",
    "summary_trans_ir_geo_6-311+g3df2p.png",
    rows,
)
overview(
    "%s__6-311+g3df2p/ir_summary_gauche.txt",
    "Gauche, IR, Consistent, 6-311+G(3df,2p)",
    "summary_gauche_ir_consistent_6-311+g3df2p.png",
    rows,
)
overview(
    "%s__6-311+g3df2p/ir_summary_trans.txt",
    "Trans, IR, Consistent, 6-311+G(3df,2p)",
    "summary_trans_ir_consistent_6-311+g3df2p.png",
    rows,
)

print >> f, "10*Log10 of ratio between theoretical and experimental rate constant."
print >> f, "<table style='border-color:black'>"
for row in rows:
    print >> f, "<tr>%s</tr>" % ("".join(row).encode('UTF-8'))
print >> f, "</table>"

print >> f, "</body></html>"

