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


import glob, os
from optparse import OptionParser

from states import *

usage = """
%prog create lot_label basis_label
%prog clone geo_lot_label geo_basis_label sp_lot_label sp_basis_label
"""

def main():
    parser = OptionParser(usage)
    parser.add_option(
        "-r", "--random", default=False, action='store_true',
        help="Randomly distort the initial geometries with uniform "
             "displacements of [-0.1*angstrom,0.1*angstrom] in x, y, and z "
             "directions."
    )
    parser.add_option(
        "-s", "--suffix", default="",
        help="Append __SUFFIX to the directory name"
    )
    options, args = parser.parse_args()

    if len(args) == 0:
        parser.error("Expecting at least one argument")

    if args[0] == "create":
        if len(args) != 3:
            parser.error("create needs two extra arguments")
        lot_label, basis_label = args[1:]
        root = get_root(lot_label, basis_label, options.suffix)
        for state in states:
            for job in state.jobs:
                if job.name == "sp":
                    continue
                job.write_input(state, root, lot_label, basis_label, options.suffix, options.random)
    elif args[0] == "clone":
        if len(args) != 5:
            parser.error("create needs four extra arguments")
        lot0_label, basis0_label, lot1_label, basis1_label = args[1:]
        root0 = get_root(lot0_label, basis0_label, options.suffix)
        root1 = "GEO__%s__ENERGY__%s" % (
            get_root(lot0_label, basis0_label, ""),
            get_root(lot1_label, basis1_label, options.suffix)
        )
        if not os.path.isdir(root1):
            os.mkdir(root1)
        for dirname in glob.glob("%s/*" % root0):
            if os.path.isdir(dirname) and not (dirname.endswith("sp") or dirname.endswith("bsse")):
                jobdir = os.path.basename(dirname)
                source = "../%s/%s" % (root0, jobdir)
                destination = "%s/%s" % (root1, jobdir)
                if not os.path.exists(destination):
                    os.symlink(source, destination)
        for state in states:
            for job in state.jobs:
                if job.name == "bsse" or job.name == "sp" or job.name.startswith("cps"):
                    job.write_input(state, root1, lot1_label, basis1_label, options.suffix, options.random)
    else:
        parser.error("Unknown command: %s" % args[0])



if __name__ == "__main__":
    main()
