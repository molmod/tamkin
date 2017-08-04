#!/bin/bash
# TAMkin is a post-processing toolkit for thermochemistry and kinetics analysis.
# Copyright (C) 2008-2010 Toon Verstraelen <Toon.Verstraelen@UGent.be>,
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


if [ -z "$1" ]; then
  echo "One arument is required, e.g. alkane_2"
  exit
fi

cd $1/md

# convert the cp2k output to a binary tracks database
tr-from-cp2k-ener md-1.ener

# make plots of the energies and the temperature
tr-plot --xunit=ps --xlabel="Time" --yunit=kjmol --ylabel="Energy" \
  :line tracks/time tracks/potential_energy -l E_pot \
  :line tracks/time tracks/conserved_quantity -l E_cons \
  ../energies.png
tr-plot --xunit=ps --xlabel="Time" --yunit=K --ylabel="Temperature" \
  :line tracks/time tracks/temperature \
  ../temperatures.png

# compute the variance on the kinetic energy
tr-fluct tracks/kinetic_energy tracks/kinetic_energy.fluct
echo "The parameters below are required to compute the heat capacity."
echo "Variance on the kinetic energy [a.u.], statistical error [a.u.], sampling efficiency, correlation time [ps]"
tr-blav tracks/kinetic_energy.fluct tracks/time -t ps

# compute the average temperature
echo "Average temperature [K], statistical error [K], sampling efficiency, correlation time [ps]"
tr-blav tracks/temperature tracks/time -t ps -u K



