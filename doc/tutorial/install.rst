Installation of TAMkin software
===============================


Preparing your mind
-------------------

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


Preparing your Linux system
---------------------------

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

Most Linux distributions can install this software with just a single command
on the command line by the administrator. They are listed below for several
popular Linux distributions:

* Ubuntu 10.4::

    sudo apt-get install python python-dev python-numpy python-scipy gfortran gcc git-core

* Debian 5. You first have to become root because the sudo program is not
  configured by default. ::

    su -
    apt-get install python python-dev python-numpy python-scipy gfortran gcc git-core
    exit

* Fedora 12 and 13. You first have to become root because the sudo program is
  not configured by default. ::

    su -
    pkcon install python-devel numpy numpy-f2py scipy gcc-gfortran gcc git
    exit

* Suse 11.2::

    sudo zypper ar http://download.opensuse.org/repositories/devel:/languages:/python/openSUSE_11.2/devel:languages:python.repo
    sudo zypper install python-devel python-numpy python-scipy gcc gcc-fortran git

  There seems to be something odd going on with the default Python configuration
  on Suse installations. You have to edit the file
  ``/usr/lib64/python2.4/distutils/distutils.cfg`` or
  ``/usr/lib32/python2.4/distutils/distutils.cfg``, depending on the CPU
  architecture, to comment out the line ``prefix=/usr/local`` with a ``#``
  symbol. Otherwise it is impossible to install Python packages in the home
  directory, as we will do below.

In order to enable the installation and usage of Python packages in the home
directory, as we will do in the next section, one must configure a few
environment variables:

* Bash users: add the following two lines to your ``~/.bashrc`` file::

    export PYTHONPATH=$HOME/lib/python:$HOME/lib64/python:$PYTHONPATH
    export PATH=$HOME/bin:$PATH

* TC Shell users: add the lines to your ``~/.tcshrc`` file::

    setenv PYTHONPATH $HOME/lib/python:$HOME/lib64/python:$PYTHONPATH
    setenv PATH $HOME/bin:$PATH

If you don't know which shell you are using, you are probably using Bash. Note
that some of these lines may already be present. **These settings are only
loaded in new terminal sessions, so close your terminal and open a new one
before proceeding.**


Installing the bleeding edge version of TAMkin
----------------------------------------------

The following series of commands will download the latest versions of the
MolMod package (required) and TAMkin, and will then install them into your
home directory. Make sure you execute these commands in some sort of temporary
directory. ::

    git clone git://molmod.ugent.be/git/molmod.git
    git clone git://molmod.ugent.be/git/tamkin.git
    (cd molmod; ./setup.py install --home=~)
    (cd tamkin; ./setup.py install --home=~)

You are now ready to start using TAMkin!

A few quick checks
------------------

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

    toon@poony ~> python
    Python 2.6.5 (r265:79063, Apr 16 2010, 13:57:41)
    [GCC 4.4.3] on linux2
    Type "help", "copyright", "credits" or "license" for more information.
    >>> import tamkin
    >>> import molmod
    >>> quit()
    toon@poony ~>


Upgrading to the bleeding edge version of TAMkin
------------------------------------------------

In case you want to upgrade TAMkin to the latests development version after a
previous install, then execute the following commands (in the same directory)::

    (cd molmod; git pull; rm -r ~/lib*/python/molmod*; ./setup.py install --home=~)
    (cd tamkin; git pull; rm -r ~/lib*/python/tamkin*; ./setup.py install --home=~)


Testing your installation
-------------------------

The unit tests included in the source tree can be used to check if various
components of TAMkin produces the correct results::

    cd tamkin/test
    ./test.py

If some tests fail, post the output of the tests on the `mailing list
<http://molmod.ugent.be/code/wiki/TAMkin/MailingList>`_.
