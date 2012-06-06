Installation of TAMkin software
###############################


Disclaimer
==========

TAMkin is developed and tested in modern Linux environments. The
installation and usage will therefore be relatively easy on Linux. If you want
to use TAMkin on other operating systems such as Windows or OSX, you should
have a minimal computer geek status to get it working, and this document may
provide at best some clues. We are always interested in hearing from your
installation adventures.

Also note that TAMkin undergoes regular updates and improvements. You are
supposed to simply check out the latests development version and try it out.
There may be obvious and obscure bugs, although we try hard to get the code as
reliable as possible. TAMkin is bundled with a large set of unit tests to
validate the code, but that does not mean bugs are impossible. Note that this is
Open Source software and that no warranty of any kind is implied or expressed.


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

* Python 2.5, 2.6 or 2.7 (including the header files): http://www.python.org/doc/
* Numpy 1.0 or later: http://numpy.scipy.org/
* Scipy 0.6 or later: http://www.scipy.org/
* A Fortran and a C compiler supported by the F2PY module in Numpy, e.g.
  gfortran and gcc: http://gcc.gnu.org/
* Git: http://git-scm.com/

Most Linux distributions can install this software with just a single command.
They are listed below for several popular Linux distributions:

* Ubuntu 12.4::

    sudo apt-get install python python-dev python-numpy python-scipy gfortran gcc git-core

* Debian 5. You first have to become root because the sudo program is not
  configured by default. ::

    su -
    apt-get install python python-dev python-numpy python-scipy gfortran gcc git-core
    exit

* Fedora 17.::

    sudo yum install python-devel numpy numpy-f2py scipy gcc-gfortran gcc git

* Suse 11.2::

    sudo zypper ar http://download.opensuse.org/repositories/devel:/languages:/python/openSUSE_11.2/devel:languages:python.repo
    sudo zypper install python-devel python-numpy python-scipy gcc gcc-fortran git


Installing the latest version of TAMkin
=======================================

The following series of commands will download the latest version of TAMkin, and
will then install it into your home directory.  ::

    cd ~/build/
    git clone git://github.com/molmod/tamkin.git
    (cd tamkin; ./setup.py install --home=~)

You are now ready to start using TAMkin!


A few quick checks
==================

It may be interesting to double check your installation before proceeding,
unless you `feel lucky`. The TAMkin and MolMod files are installed in the
following directories:

* Scripts: ``~/bin``
* Modules: ``~/lib/python`` or ``~/lib64/python``
* Data: ``~/share``

There should be at least some files present in these directories.

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


Upgrading to the latest version of MolMod and TAMkin
====================================================

In case you want to upgrade TAMkin to the latests development version after a
previous install, then execute the following commands (in the same directory
that was originall used to install TAMkin)::

    cd ~/build/
    (cd molmod; git pull; rm -r ~/lib*/python/molmod*; ./setup.py install --home=~)
    (cd tamkin; git pull; rm -r ~/lib*/python/tamkin*; ./setup.py install --home=~)


Testing your installation
=========================

For the development and testing one needs to install one additional package:

* Nosetests >= 0.11: http://somethingaboutorange.com/mrl/projects/nose/0.11.2/

Most Linux distributions can install this software with just a single command on
the command line by the administrator.

* Ubuntu 10.4::

    sudo apt-get install python-nose

* Debian 5::

    su -
    apt-get install python-nose
    exit

* Fedora 17::

    sudo yum install python-nose

* Suse 11.2::

    sudo zypper install python-nose

Once these dependecies are installed, execute the following commands to run the
tests::

    cd ~/build/
    cd tamkin
    nosetests -v

If some tests fail, post the output of the tests on the `TAMkin
mailing list <https://groups.google.com/forum/#!forum/tamkin>`_.

