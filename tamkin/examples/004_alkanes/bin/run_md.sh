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

if [ ! -f $1/opt.xyz ]; then
  echo "Run the static simulation first with run_static.sh"
  exit
fi

echo "Starting MD simulation. This will take a few minutes. Sit back and relax."

# We assume that the optimized geometry is in $1/opt.xyz and that the 
# file init.psf has been generated before.
rm -rf $1/md # clean up the previous run
mkdir $1/md # create a directory for the new run
cp templates/md.inp $1/md/ # copy the template file
(
   cd $1/md/; 
   cp2k.sopt -i md.inp -o md.out;
   rm *.restart
)



