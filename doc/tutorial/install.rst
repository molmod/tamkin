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

Most Linux distributions can install most required software with a single command.

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

Execute the following command to install TAMkin in your home directory.

.. code:: bash

    pip install tamkin


Testing your installation
=========================

Just run, after installation, the following command to test your copy of TAMkin:

.. code:: bash

    nosetests -v tamkin

If some tests fail, post an issue on https://github.com/molmod/tamkin/issues
