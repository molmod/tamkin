.. image:: https://travis-ci.org/molmod/tamkin.svg?branch=master
    :target: https://travis-ci.org/molmod/tamkin

TAMkin is a post-processing toolkit for normal mode analysis, thermochemistry
and reaction kinetics. It uses a Hessian computation from a standard
computational chemistry program as its input. CHARMM, CP2K, CPMD, GAMESS,
GAUSSIAN, QCHEM and VASP are supported. Multiple methods are implemented to
perform a normal mode analysis (NMA). The frequencies from the NMA can be used
to construct a molecular partition function to derive thermodynamic and kinetic
parameters.

More information about TAMkin can be found on the CMM Software website:
http://molmod.ugent.be/software/

TAMkin is distributed as open source software under the conditions of the GPL
license version 3.  Read the file COPYING for more details, or visit
http://www.gnu.org/licenses/


Installation
============

TAMkin can be installed with pip (system wide or in a virtual environment):

.. code:: bash

    pip install tamkin

Alternatively, you can install TAMkin in your home directory:

.. code:: bash

    pip install tamkin --user

Lastly, you can also install TAMkin with conda. (See
https://www.continuum.io/downloads)

.. code:: bash

    conda install -c tovrstra tamkin


Testing
=======

The tests can be executed as follows:

.. code:: bash

    nosetests tamkin
