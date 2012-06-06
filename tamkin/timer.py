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
"""Timer to keep track of wall and cpu time"""


import time, sys


__all__ = ["Timer"]


class Timer(object):
    """Timer object to keep track of the time spent in several parts of TAMKin.

       This timer works like a stopwatch. Each time the :meth:`sample` method is
       called the current cpu and wall times are recorded toghether with a label

       The methods :meth:`dump` and :meth:`write_to_file` can be used to
       generate a report.
    """

    def __init__(self):
        self.cpu_times = []
        self.wall_times = []
        self.labels = []

    def sample(self, label):
        """Record the current timings and associate them with the given label

           Argument:
            | label  --  A string describing this point in the code.
        """
        self.cpu_times.append(time.clock())
        self.wall_times.append(time.time())
        self.labels.append(label)

    def dump(self, f=sys.stdout):
        """Dump the logfile with timing information, to screen or to a file stream.

           Optional argument:
            | f  --  the stream to write to. [default=sys.stdout]
        """
        print >> f, "-------------------"
        print >> f, "Printing LOG jobtimer"
        print >> f, '%12s %12s %21s %16s %30s' %("cpu times [s]", "diff [s]", "wall times [s]", "diff [s]", "labels" )
        for i,label in enumerate(self.labels[:-1]):
            print >> f, '%12.3f %12.3f %21.3f %16.3f %30s' %(self.cpu_times[i],
                                         self.cpu_times[i+1]-self.cpu_times[i],
                                         self.wall_times[i],
                                         self.wall_times[i+1]-self.wall_times[i],
                                         label)
        print >> f, '%12.3f %12s %21.3f %16s %30s' %(self.cpu_times[-1], "",
                                         self.wall_times[-1], "",
                                         self.labels[-1])
        print >> f, "-------------------"

    def write_to_file(self, filename):
        """Write the logfile with timing information to filename.

           Argument:
            | filename  --  the file to write to.
        """
        f = file(filename, 'w')
        self.dump(f)
        f.close()
