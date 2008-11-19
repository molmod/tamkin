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


from tamkin.data import Molecule

from molmod.io.gaussian03.fchk import FCHKFile
from molmod.units import amu

import numpy


__all__ = ["load_fixed_g03com", "load_molecule_g03fchk"]


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


def load_molecule_g03fchk(filename_freq,filename_ener=None): # if one file contains every information, give the name twice
    fchk_freq = FCHKFile(filename_freq, ignore_errors=True, field_labels=[
        "Cartesian Force Constants", "Real atomic weights", "Total Energy", "Multiplicity"
    ])
    if filename_ener is None:
        fchk_ener = fchk_freq
    else:
        fchk_ener = FCHKFile(filename_ener, ignore_errors=True, field_labels=[
            "Total Energy"
        ])
    return Molecule(
        fchk_freq.molecule.numbers,
        fchk_freq.molecule.coordinates,
        fchk_freq.fields["Real atomic weights"]*amu,
        fchk_ener.fields["Total Energy"],
        fchk_freq.get_hessian(),
        fchk_freq.fields["Multiplicity"],
    )

