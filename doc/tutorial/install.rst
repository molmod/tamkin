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

Installation of TAMkin software
###############################


Disclaimer
==========

TAMkin is mainly developed and tested on Linux systems. If you run any other
operating system, some of the instructions below may not work.


External depedencies
====================

Some other software packages should be installed before TAMkin can be installed
or used. It is recommended to use the software package management of your Linux
distribution to install these dependencies.

The following software must be installed for TAMkin:

* Python >=2.7,<3.0: http://www.python.org/doc/
* Numpy >=1.0: http://numpy.scipy.org/
* Scipy >=0.17.1: http://www.scipy.org/
* MatPlotLib >=1.1: http://matplotlib.org/
* MolMod >=1.3.1: https://github.com/molmod/molmod
* Nosetests >=0.11: http://somethingaboutorange.com/mrl/projects/nose/0.11.2/

Most Linux distributions can install most required software with a single command. Except
for Python, missing dependencies will get automatically installed by pip. (See below.)

* **Ubuntu**

    sudo apt-get install python python-numpy python-scipy python-matplotlib

* **Debian 5** or newer. You first have to become root because the sudo program is not
  configured by default::

    su -
    apt-get install python python-numpy python-scipy python-matplotlib
    exit

* **Fedora 17 to 22**

    sudo yum install python numpy scipy python-matplotlib


Installation
============

TAMkin can be installed with pip (system wide or in a virtual environment):

.. code:: bash

    pip install numpy
    pip install tamkin

Alternatively, you can install TAMkin in your home directory:

.. code:: bash

    pip install numpy --user
    pip install tamkin --user

Lastly, you can also install MolMod with conda. (See
https://www.continuum.io/downloads)

.. code:: bash

    conda install -c molmod tamkin


Testing your installation
=========================

Just run, after installation, the following command to test your copy of TAMkin:

.. code:: bash

    nosetests -v tamkin

If some tests fail, post an issue on https://github.com/molmod/tamkin/issues
