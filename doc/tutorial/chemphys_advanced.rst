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

Chemical Physics -- Advanced TAMkin examples
============================================


**TODO**: add introductory text.

Intrinsic rate constants of reactions in zeolites
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Intrinsic rate constants are typically modeled in a two-step approach where one
first describes the reaction in a small (e.g. 5T) cluster model of the zeolite.
Given the reactant and transition state geometry with the small cluster, one can
construct an initial estimate of the corresponding geometries in a larger
cluster using `Zeobuilder <http://molmod.ugent.be/code/wiki/Zeobuilder>`_. These
structures can be refined with ONIOM or another type of QM/MM computations.

In the example below we study the addition of ethene to an ethyl alkoxide on an
acid site in ZSM-5. The product of this reaction is methylcyclopropane. For the
forward reaction rate no computations on the product state are required.


Small cluster computation (5T)
------------------------------

The reactant and the transition state are depicted in the figures below

.. image:: 5t.png

The two checkpoint files in the example are derived from frequency computations
with Gaussian03 using the B3LYP/6-31G(d) level of theory at the reactant and
transition state.

.. literalinclude:: ../../tamkin/examples/006_5T_ethyl_ethene_addition/reaction.py
   :lines: 37-
   :linenos:
   :caption: tamkin/examples/006_5T_ethyl_ethene_addition/reaction.py

Take note of the following (subtle) choices in the input file:

Line 8:
    The transition state geometry is not perfectly optimized. This would be a
    time-consuming computation and not really required for setting up a draft
    model. The gradient_threshold parameter is increased from ``1e-4`` to
    ``2e-4``.

Line 11-12:
    The translational and rotational contribution can be omitted as these
    reactions take place on a fixed active site in a zeolite catalyst. Free
    rotation and translation of the 5T fragment with the guest molecules is
    not possible. In principle these missing degrees of freedom in the partition
    function should be replaced by a vibrational coupling of the cluster with
    the surrounding framework. This is coupling is essentially neglected in
    a small cluster computation.




Large cluster computation
-------------------------




DIY Kinetic parameter estimation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**TODO**: show how TAMkin can be used to make a powerful ``freqchk``-like
program.


Many similar reactions with a single script
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**TODO**: explain the for-loop.


Basic kinetics level of theory study
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**TODO**


Smart TAMkin script
~~~~~~~~~~~~~~~~~~~

**TODO**: explain how to write a smart script with TAMkin that detects, based on
a few filename convetions, what type of computation must be carried out.
(Specific for Gaussian.)


Heat of formation
~~~~~~~~~~~~~~~~~

**TODO**: general script for the computation of the heat of formation.


Pattern matching
~~~~~~~~~~~~~~~~

**TODO**: how to use the ``molmod.graphs`` module to locate a set of frozen
atoms in a series of similar systems with unknown atom order.
