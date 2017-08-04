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

python bin/xyz2psf.py $1/init.xyz $1/init.psf

rm -rf $1/opt # clean up the previous run
mkdir $1/opt # create a directory for the new run
cp templates/opt.inp $1/opt/ # copy the template file
(
   cd $1/opt/; 
   cp2k.sopt -i opt.inp -o opt.out;
   tail -n $(( $(head -n1 opt-pos-1.xyz) +2 )) opt-pos-1.xyz > ../opt.xyz
)

rm -rf $1/sp
mkdir $1/sp
cp templates/sp.inp $1/sp/
(
   cd $1/sp/; 
   cp2k.sopt -i sp.inp -o sp.out;
)

rm -rf $1/freq
mkdir $1/freq
cp templates/freq.inp $1/freq/
(
   cd $1/freq/; 
   cp2k.sopt -i freq.inp -o freq.out;
)


