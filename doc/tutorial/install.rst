Installation of TAMkin software
###############################


Disclaimer
==========

TAMkin is mainly developed and tested on Linux systems. If you run any other
operating system, some of the instructions below may not work.


MolMod dependency
=================

`MolMod <http://molmod.github.com/molmod/>`_ is a Python library used by most
Python programs developed at the CMM. It must be installed before TAMkin can
be used or installed. Installation and download instructions can be found in the
`molmod documentation <http://molmod.github.com/molmod/tutorial/install.html>`_.
The instructions below only work of the MolMod package is installed.


External depedencies
====================

Some other software packages should be installed before TAMkin can be installed
or used. It is recommended to use the software package management of your Linux
distribution to install these dependencies.

The following software must be installed for TAMkin:

* Python 2.6 or 2.7: http://www.python.org/doc/
* Numpy 1.0 or newer: http://numpy.scipy.org/
* Scipy 0.6 or newer: http://www.scipy.org/
* MatPlotLib 1.1 or newer: http://matplotlib.org/

Most Linux distributions can install this software with just a single command.
They are listed below for several popular Linux distributions:

* **Ubuntu 12.4** or newer::

    sudo apt-get install python python-numpy python-scipy python-matplotlib

* **Debian 5** or newer. You first have to become root because the sudo program is not
  configured by default::

    su -
    apt-get install python python-numpy python-scipy python-matplotlib
    exit

* **Fedora 17** or newer::

    sudo yum install python numpy scipy python-matplotlib


Download the code
=================

Stable release (recommended)
----------------------------

The latest stable source code release can be downloaded here:

    http://github.com/molmod/tamkin/releases/download/v1.0.9/TAMkin-1.0.9.tar.gz

Choose a suitable directory, e.g. ``~/build``, download and unpack the archive::

    mkdir -p ~/build
    cd ~/build
    wget http://github.com/molmod/tamkin/releases/download/v1.0.9/TAMkin-1.0.9.tar.gz
    tar -xvzf TAMkin-1.0.9.tar.gz
    cd TAMkin-1.0.9


Latest development code (experts only)
--------------------------------------

In order to get the latest development version of the source code, you need to
work with git. Git is a version control system
that makes life easy when a group of people are working on a common source code.
All information about git (including downloads and tutorials) can be found here:
http://git-scm.com/. The official git URL of TAMkin is:
git://github.com/molmod/tamkin.git. In order to `clone` the public TAMkin
repository, run this command::

    git clone git://github.com/molmod/tamkin.git
    cd tamkin

The version history can be updated with the latest patches with the following
command::

    git pull

There is also a web interface to TAMkin's git repository:
https://github.com/molmod/tamkin


Installing the latest version of TAMkin
=======================================

Execute the following command in the TAMkin source directory to install TAMkin
in your home directory. ::

    ./setup.py install --user

You are now ready to start using TAMkin!


A few quick checks
==================

The Python modules should be accessible from any Python session. This can be
checked by starting Python interactively and loading the modules manually. There
should be no errors when importing the modules::

    $ python
    Python 2.6.5 (r265:79063, Apr 16 2010, 13:57:41)
    [GCC 4.4.3] on linux2
    Type "help", "copyright", "credits" or "license" for more information.
    >>> import tamkin
    >>> import molmod
    >>> quit()


Testing your installation
=========================

For the development and testing, one needs to install one additional package:

* Nosetests >= 0.11: http://somethingaboutorange.com/mrl/projects/nose/0.11.2/

Most Linux distributions can install this software with just a single terminal command:

* **Ubuntu 12.4** or newer::

    sudo apt-get install python-nose

* **Debian 5** or newer::

    su -
    apt-get install python-nose
    exit

* **Fedora 17** or newer::

    sudo yum install python-nose


Once these dependecies are installed, execute the following command in the
TAMkin source tree to run the tests::

    nosetests -v test

If some tests fail, post the output of the tests on the `TAMkin
mailing list <https://groups.google.com/forum/#!forum/tamkin>`_.
