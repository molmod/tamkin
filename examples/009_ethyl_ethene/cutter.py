#!/usr/bin/env python

# this tiny program cuts away all unused lines from a gaussian log file of
# a torsional scan computation. All torsional scan files in TAMkin are reduced
# in size with this program.

import sys

active = 3
for line in sys.stdin:
    line = line[:-1]
    if line == "                          Input orientation:                          ":
        active = 3
    if active > 0 or line.startswith(" SCF Done:") or line.startswith("    -- Stationary point found."):
        print line
    if line == " ---------------------------------------------------------------------":
        active -= 1
