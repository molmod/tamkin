#!/usr/bin/env python
# -*- coding:utf8 -*-

from lot_basis import lots_list

from molmod.io.gaussian03.fchk import FCHKFile
from molmod.units import angstrom, deg
import molmod.ic as ic

def bond_length(mol, i0, i1):
    c = mol.coordinates
    return ic.bond_length(c[i0], c[i1])[0]

def bend_angle(mol, i0, i1, i2):
    c = mol.coordinates
    return ic.bend_angle(c[i0], c[i1], c[i2])[0]


def generic_geom(ic_descs, fn_fchk, name, basis_label, rows):
    if len(rows) == 0:
        rows.append([u"<th>Approach→<br />Functional↓</th>"])
        rows.append(["<th></th>"])
    rows[0].append("<th colspan=%i>%s<br />%s</th>" % (len(ic_descs), name, basis_label))
    for ic_desc in ic_descs:
        rows[1].append(ic_desc[0])
    counter = 2
    for lot in lots_list:
        if lot.spin == "ROS":
            lot_label = "ro" + lot.label
        else:
            lot_label = lot.label
        if len(rows) <= counter:
            rows.append(["<th>%s</th>" % lot_label])
        try:
            mol = FCHKFile(fn_fchk % (lot_label, basis_label), field_labels=[]).molecule
        except IOError:
            mol = None
        if mol is None:
            for ic_desc in ic_descs:
                rows[counter].append("<td>&nbsp</td>")
        else:
            for ic_desc in ic_descs:
                cell_format, ic_fn, indexes, unit = ic_desc[1:]
                rows[counter].append(cell_format % (ic_fn(mol, *indexes)/unit))
        counter += 1

def ts_geom(ts_name, basis_label, rows):
    fn_fchk = "%%s__%%s/ts_ad1_%s__opt/gaussian.fchk" % ts_name.lower()
    ic_descs = [
        (u"<td>d<sub>12</sub>/Å</td>", "<td>%.2f</td>", bond_length, (1, 2), angstrom),
        (u"<td>a<sub>012</sub>/°</td>", "<td>%.1f</td>", bend_angle, (0, 1, 2), deg),
        (u"<td>a<sub>123</sub>/°</td>", "<td>%.1f</td>", bend_angle, (1, 2, 3), deg),
    ]
    generic_geom(ic_descs, fn_fchk, ts_name, basis_label, rows)

def ra_geom(mol_name, basis_label, rows):
    fn_fchk = "%%s__%%s/%s__opt/gaussian.fchk" % mol_name.lower()
    ic_descs = [
        (u"<td>d<sub>01</sub>/Å</td>", "<td>%.2f</td>", bond_length, (0, 1), angstrom),
    ]
    generic_geom(ic_descs, fn_fchk, mol_name, basis_label, rows)


f = file("geotab.html", "w")
print >> f, "<?xml version='1.0' encoding='UTF-8'?>"
print >> f, "<html><head><title>Geometry Overview</title>"
print >> f, "<style type='text/css'>"
print >> f, "body { font-family: sans; }"
print >> f, "td { text-align: right; width: 50px; padding: 0px 10px; }"
print >> f, "table, td, th, tr { border: solid 1px black; }"
print >> f, "table { table-layout: fixed; border-collapse:collapse; }"
print >> f, "</style>"
print >> f, "</head><body>"

rows = []
ts_geom("Gauche", "6-31gd", rows)
ts_geom("Trans", "6-31gd", rows)
ts_geom("Gauche", "6-311+g3df2p", rows)
ts_geom("Trans", "6-311+g3df2p", rows)

print >> f, "<p>Geometrical parameters related to the transition state.</p>"
print >> f, "<table style='border-color:black'>"
for row in rows:
    print >> f, "<tr>%s</tr>" % ("".join(row).encode('UTF-8'))
print >> f, "</table>"

rows = []
ra_geom("Ethene", "6-31gd", rows)
ra_geom("Ethyl", "6-31gd", rows)
ra_geom("Ethene", "6-311+g3df2p", rows)
ra_geom("Ethyl", "6-311+g3df2p", rows)

print >> f, "<p>Geometrical parameters related to the reactants.</p>"
print >> f, "<table style='border-color:black'>"
for row in rows:
    print >> f, "<tr>%s</tr>" % ("".join(row).encode('UTF-8'))
print >> f, "</table>"

print >> f, "</body></html>"

