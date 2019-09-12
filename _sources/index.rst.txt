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

The TAMkin Manual
=================

TAMkin is a post-processing toolkit for normal mode analysis, thermochemistry
and reaction kinetics. It uses a Hessian computation from a standard
computational chemistry program as its input. CHARMM, CP2K, CPMD, GAMESS,
GAUSSIAN, QCHEM and VASP are supported. Multiple methods are implemented to
perform a normal mode analysis (NMA). The frequencies from the NMA can be used
to construct a molecular partition function to derive thermodynamic and kinetic
parameters.

If you use TAMkin for your research or for the preparation of publications,
please cite the following paper: Ghysels, A.; Verstraelen, T.; Hemelsoet, K.;
Waroquier, M.; Van Speybroeck, V. *J. Chem. Inf. Model.* **2010**, 50,
1736-1750.
(`http://dx.doi.org/10.1021/ci100099g <http://dx.doi.org/10.1021/ci100099g>`_).

A release history can be found here: :ref:`releases`.


Tutorial
--------

.. toctree::
   :maxdepth: 1
   :numbered:

   tutorial/install.rst
   tutorial/support.rst
   tutorial/cite.rst
   tutorial/getting_started.rst
   tutorial/data.rst
   tutorial/chemphys_theory.rst
   tutorial/chemphys_tamkin.rst
   tutorial/chemphys_examples.rst
   tutorial/chemphys_advanced.rst
   tutorial/development.rst

Library Reference
-----------------

.. toctree::
   :maxdepth: 1
   :numbered:

   reference/tamkin-driver.rst
   reference/io.rst
   reference/nma.rst
   reference/pf.rst
   reference/aux.rst
   references.rst
   releases.rst


Indices, tables and refereces
-----------------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
