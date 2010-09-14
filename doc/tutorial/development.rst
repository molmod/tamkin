TAMkin programmer's guide
=========================

Work in progres ...

Have a look at the structure of the package
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It may be useful to understand the basic structure of the TAMkin package. Assume
TAMkin is installed in a directory ``~/code/``, then one can have a look at the
TAMkin files in the following way::

    $ cd ~/code/tamkin
    $ ls
    cleancode.sh   COPYING   HEADER.py   lib       test
    cleanfiles.sh  examples  install.sh  setup.py  uninstall.sh

* The ``install.sh`` script uses the ``setup.py`` script to install TAMkin
  (see the https://molmod.ugent.be/code/wiki/TAMkin/InstallationGuide).
  Similarly, the ``uninstall.sh`` script can be used to uninstall TAMkin.
* The ``cleancode.sh``, ``cleanfiles.sh`` and ``HEADER.py`` files are only
  relevant for code developers, but not for the general user.
* The file ``COPYING`` contains the license under which TAMKin is distributed.
* The directory ``test/`` contains the testing routines. This is mainly intented
  for developers, but a regular user can also use it to test the validity of its
  installation. ::

    $ ls test/
    input  nma.py       partf.py    rotor.py  timer.py
    io.py  nmatools.py  pftools.py  test.py   tunneling.py

* The directory ``examples/`` contains worked out examples (see further).
* The directory ``lib/`` contains the TAMkin source code itself (the files with
  lines of code, for developers). ::

    $ ls lib/
    data.py  __init__.py  nma.py       partf.py    rotor.py  tunneling.py
    geom.py  io           nmatools.py  pftools.py  timer.py
