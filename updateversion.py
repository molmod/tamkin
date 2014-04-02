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

import re, sys

rules = [
    ('setup.py', '^    version=\'(...+)\',$'),
    ('tamkin/__init__.py', '^__version__ = \'(...+)\'$'),
    ('doc/conf.py', '^version = \'(...+)\'$'),
    ('doc/conf.py', '^release = \'(...+)\'$'),
    ('doc/tutorial/install.rst', '^    http://github.com/molmod/tamkin/releases/download/v(...+)/TAMkin-(...+).tar.gz'),
    ('doc/tutorial/install.rst', '^    wget http://github.com/molmod/tamkin/releases/download/v(...+)/TAMkin-(...+).tar.gz'),
    ('doc/tutorial/install.rst', '^    tar -xvzf TAMkin-(...+).tar.gz$'),
    ('doc/tutorial/install.rst', '^    cd TAMkin-(...+)$'),
]

if __name__ == '__main__':
    newversion = sys.argv[1]

    for fn, regex in rules:
        r = re.compile(regex)
        with open(fn) as f:
            lines = f.readlines()
        for iline in xrange(len(lines)):
            line = lines[iline]
            m = r.match(line)
            if m is not None:
                for igroup in xrange(m.lastindex, 0, -1):
                    line = line[:m.start(igroup)] + newversion + line[m.end(igroup):]
                lines[iline] = line
        with open(fn, 'w') as f:
            f.writelines(lines)
