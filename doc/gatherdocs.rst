Input/output
============

.. automodule:: tamkin.io
   :members:

Tools to read data from general packages
----------------------------------------

.. automodule:: tamkin.io.charmm
   :members:

.. automodule:: tamkin.io.cp2k
   :members:

.. automodule:: tamkin.io.cpmd
   :members:

.. automodule:: tamkin.io.gaussian
   :members:

.. automodule:: tamkin.io.qchem
   :members:

.. automodule:: tamkin.io.vasp
   :members:


Tools for internal usage in Tamkin
----------------------------------

.. automodule:: tamkin.io.internal
   :members:

Tools to generate trajectories
------------------------------

.. automodule:: tamkin.io.trajectory
   :members:


Normal mode analysis
====================

The different NMA models
------------------------

.. automodule:: tamkin.nma
   :members:

Tools to analyze frequencies and modes
--------------------------------------

Some additional tools are in the ``nmatools.py`` code.

.. automodule:: tamkin.nmatools
   :members:


Partition functions
===================

The molecular partion function
------------------------------

.. automodule:: tamkin.partf
   :members:

Tools to analyze partition functions
------------------------------------

For instance, a thermodynamical analysis, the
calculation of the equilibrium constant of a
reaction, the reaction rate constant, the reaction
parameters of the Arrhenius plot, ...

.. automodule:: tamkin.pftools
   :members:

The 1-D rotor
-------------

The description of one-dimensional hindered rotors is also implemented in
TAMkin, in the module ``rotor.py``. It calculates the partition function of
a 1-D rotor.

.. automodule:: tamkin.rotor
   :members:

Tunneling effects
-----------------

This module calculates the correction to the partition function due
to (quantum) tunneling through the reaction barrier.

.. automodule:: tamkin.tunneling
   :members:


Auxiliary modules
=================

Some more functions/classes.

.. automodule:: tamkin.__init__
   :members:

.. automodule:: tamkin.timer
   :members:

.. automodule:: tamkin.data
   :members:

.. automodule:: tamkin.geom
   :members:

