# TAMkin is a post-processing toolkit for thermochemistry and kinetics analysis.
# Copyright (C) 2008 Toon Verstraelen <Toon.Verstraelen@UGent.be>, 
# Matthias Vandichel <Matthias.Vandichel@UGent.be> and
# An Ghysels <An.Ghysels@UGent.be>
#
# This file is part of TAMkin.
#
# TAMkin is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# TAMkin is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
# --


def load_fixed_g03com(filename):
    f = file(filename)
    for line in f:
        # iterate until we reach the title line
        line = line.strip()
        if not (len(line)==0 or line[0] == "#" or line[0] == "%"):
            break
    for line in f:
        # skip lines until charge/multiplicty is reached
        words = line.split()
        if len(words) == 2:
            try:
                int(words[0])
                int(words[1])
                break
            except ValueError:
                pass
    counter = 0
    fixed_atoms = []
    for line in f:
        # go trough the list of atoms and store the fixed ones.
        words = line.split()
        if len(words) == 0:
            break
        if words[1] == "-1":
            fixed_atoms.append(counter)
        counter += 1
    return fixed_atoms

