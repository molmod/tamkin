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
#!/usr/bin/env python

# This script takes a list of Gaussian fchk files as arguments and will strip
# all sections that are irrelevant for TAMkin. Also the trailing white space on
# each line is stripped. This is used to reduce the size of the example fchk
# files in TAMkin.


keep_fields = set([
    "Multiplicity", "Total Energy", "Atomic numbers",
    "Current cartesian coordinates", "Real atomic weights",
    "Cartesian Gradient", "Cartesian Force Constants",
])

def strip(fn):
    # a list with all lines that we'd like to keep
    lines = []
    # load the good lines
    f = file(fn, "r")
    busy = False
    keep = False
    for line in f:
        line = line.rstrip()
        if len(lines) < 2:
            # keep the first two lines
            lines.append(line)
        else:
            if busy:
                if not line.startswith(" "):
                    keep = False
                    busy = False
                if keep:
                    lines.append(line)
            if not busy:
                busy = line[47:49] == "N="
                field = line[:43].strip()
                if field in keep_fields:
                    lines.append(line)
                    keep = True
                else:
                    keep = False
    f.close()
    # print stuff back into the same file
    f = file(fn, "w")
    for line in lines:
        print >> f, line
    f.close()


if __name__ == "__main__":
    import sys
    fns_fchk = sys.argv[1:]
    for fn in fns_fchk:
        if fn.endswith(".fchk"):
            print "Stripping", fn
            strip(fn)
        else:
            print "Skipping", fn, "(wrong extension)"
