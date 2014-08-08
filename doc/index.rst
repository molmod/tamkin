
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


Indices, tables and refereces
-----------------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
* :ref:`references`
