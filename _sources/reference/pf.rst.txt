..
    : TAMkin is a post-processing toolkit for normal mode analysis, thermochemistry
    : and reaction kinetics.
    : Copyright (C) 2008-2012 Toon Verstraelen <Toon.Verstraelen@UGent.be>, An Ghysels
    : <An.Ghysels@UGent.be> and Matthias Vandichel <Matthias.Vandichel@UGent.be>
    : Center for Molecular Modeling (CMM), Ghent University, Ghent, Belgium; all
    : rights reserved unless otherwise stated.
    :
    : This file is part of TAMkin.
    :
    : TAMkin is free software; you can redistribute it and/or
    : modify it under the terms of the GNU General Public License
    : as published by the Free Software Foundation; either version 3
    : of the License, or (at your option) any later version.
    :
    : In addition to the regulations of the GNU General Public License,
    : publications and communications based in parts on this program or on
    : parts of this program are required to cite the following article:
    :
    : "TAMkin: A Versatile Package for Vibrational Analysis and Chemical Kinetics",
    : An Ghysels, Toon Verstraelen, Karen Hemelsoet, Michel Waroquier and Veronique
    : Van Speybroeck, Journal of Chemical Information and Modeling, 2010, 50,
    : 1736-1750W
    : http://dx.doi.org/10.1021/ci100099g
    :
    : TAMkin is distributed in the hope that it will be useful,
    : but WITHOUT ANY WARRANTY; without even the implied warranty of
    : MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    : GNU General Public License for more details.
    :
    : You should have received a copy of the GNU General Public License
    : along with this program; if not, see <http://www.gnu.org/licenses/>
    :
    : --

Partition functions
===================

The molecular partion function
------------------------------

**Inheritance diagram**

.. inheritance-diagram:: tamkin.partf
   :parts: 1

.. automodule:: tamkin.partf
   :members:
   :show-inheritance:

Chemical models based on partition functions
--------------------------------------------

The term `chemical models` is used as a group name for thermodynamic models of
chemical equilibria and kinetic models of chemical reactions.

.. inheritance-diagram:: tamkin.chemmod
   :parts: 1

.. automodule:: tamkin.chemmod
   :members:

Tools to analyze partition functions and kinetic models
-------------------------------------------------------

.. automodule:: tamkin.pftools
   :members:

The 1-D rotor
-------------

The description of one-dimensional hindered rotors is also implemented in
TAMkin, in the module ``rotor.py``. It calculates the partition function of
a 1-D rotor.

.. automodule:: tamkin.rotor
   :members:

Tunneling effects in chemical reactions
---------------------------------------

This module calculates the correction to the partition function due
to (quantum) tunneling through the reaction barrier.

.. automodule:: tamkin.tunneling
   :members:
