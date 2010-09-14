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
