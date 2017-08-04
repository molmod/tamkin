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


( cd 001_ethane/
  rm -v partfun.txt thermo.csv 
)

( cd 003_pentane/
  rm -v partfun1.txt thermo1.csv 
  rm -v partfun2.txt thermo2.csv 
)

( cd 002_linear_co2/
  rm -v partfun.txt thermo.csv 
)

( cd 004_alkanes/
  rm -v alkane_*/*.png
  rm -v alkane_*/*.csv
  rm -v alkane_*/*.txt
  rm -v alkane_*/init.psf
  rm -v alkane_*/opt.xyz
  rm -rv alkane_*/opt
  rm -rv alkane_*/sp
  rm -rv alkane_*/freq
  rm -rv alkane_*/md
)

( cd 005_acrylamide_reaction/
  rm -v arrhenius.png parameters.png reaction.txt
)

( cd 006_5T_ethyl_ethene_addition/
  rm -v arrhenius.png reaction.txt
)

( cd 007_mfi_propene_reaction/
  rm -v arrhenius.png parameters.png reaction.txt
)

( cd 008_ethane_rotor/
  rm -v energy_levels.png partfun.txt thermo.csv
)

( cd 009_ethyl_ethene/
  rm -v *.png *.csv *.txt
)

( cd 012_ethyl_ethene_scaling
  rm -v *.png *.csv *.txt
)

( cd 013_butane
  rm -v *.png *.csv *.txt
)

( cd 014_pentane_mbh
  rm -v *.png *.csv *.txt *.vib
)

( cd 017_activationkineticmodel
  rm -v *.png *.txt
)

( cd 018_physisorption
  rm -v *.png *.txt
)

( cd 019_ethyl_ethene_simple
  rm -v *.png *.txt *.csv
)

( cd 020_butane_conformers
  rm -v *.txt *.csv
)

( cd 021_water_formation
  rm -v *.txt *.csv
)
