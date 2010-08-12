The TAMKin Tutorial
===================

Getting started
---------------

Download and install TAMkin
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Start by downloading TAMkin on the download page. Install the TAMkin package.
Check if TAMkin is installed correctly by typing::

    $ python
    >>> import tamkin


If this results in an error, you should go through the installation process again
(see the https://molmod.ugent.be/code/wiki/TAMkin/InstallationGuide installation manual).

Have a look at the structure of the package
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It may be useful to understand the basic structure of the TAMkin package. Assume
TAMkin is installed in a directory ``~/code/``, then you can have a look at the
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


Try out the examples
~~~~~~~~~~~~~~~~~~~~

A good way to get started, is to have a look at the examples in the ``examples/``
directory of the distribution. Assume TAMkin is installed in a directory
``~/code/``, then you will find the examples on the following location::

    $ cd ~/code/tamkin/examples
    $ ls
    001_ethane               006_5T_ethene_reaction    011_ethyl_ethene_lot
    002_linear_co2           007_mfi_propene_reaction  012_ethyl_ethene_scaling
    003_pentane              008_ethane_rotor          014_pentane_mbh
    004_alkanes              009_ethyl_ethene          clean.sh
    005_acrylamide_reaction  010_adk


How a typical scripts looks like
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here is the script `thermo.py` of the first example ``examples/001_ethane/``.
In Python scripts, comments are preceded by a ``#`` sign.

::

    # Typical script
    from tamkin import *

    molecule = load_molecule_g03fchk("gaussian.fchk")
    nma = NMA(molecule)
    pf = PartFun(nma, [ExtTrans(), ExtRot()])

    # Write some general information about the molecule
    # and the partition function to a file.
    pf.write_to_file("partfun.txt")

    # Write an extensive overview of the thermodynamic properties to a file:
    ta = ThermoAnalysis(pf, [300,400,500,600])
    ta.write_to_file("thermo.csv")


The typical ingredients of your script are the following:

1. Load the TAMkin package::

        from tamkin import *

   This loads all TAMkin classes and functions, such that you can use them in
   the rest of the script.

2. Load the data::

        molecule = load_molecule_g03fchk("gaussian.fchk")

   The masses, coordinates, energy, gradient and Hessian are read from the file
   ``gaussian.fchk``. The data are stored in a ``Molecule`` object, which we give
   here the name ``molecule``.

3. Perform normal mode analysis::

        nma = NMA(molecule)

   Calculate the frequencies and normal modes of the molecule.

4. Construct a partition function::

        pf = PartFun(nma, [ExtTrans(), ExtRot()])

   Construct the partition function: construct a ``PartFun`` object, which has
   all thermodynamical properties implemented as attributes or functions. Here
   we gave it the name ``pf``. The translational and rotational contributions
   are also requested by adding ``[ExtTrans(), ExtRot()]`` as an argument. The
   vibrational contribution is always included implicitely.

5. Generate some output. For instance, ::

        pf.write_to_file("partfun.txt")

   will write the information about the partition function to a file
   ``partfun.txt``.


Generating data for TAMkin
--------------------------

TAMkin can read information from ADF, CHARMM, Gaussian, CP2K, CPMD, Q-CHEM and VASP.
Some of these programs offer freedom in the amount of output that is printed to
the output file. To produce output that TAMkin can read from, we suggest the following
settings.


CHARMM
~~~~~~

See http://www.charmm.org/

The CHARMM-input file should contain a ``VIBRAN`` section::

    CALC k = 3 * ?NATOM
    VIBRan NMODes @k
    OPEN WRITe CARD UNIT 20 NAME filename.hessian
    WRITe SECOnd CARD UNIT 20
    * title - Full Hessian
    *
    END

where the ``WRITe SECOnd`` command will write the Hessian to the file
``filename.hessian``. In addition, the statement

::

    SCALar MASS SHOW

should be called before writing the coordinates to a coordinate file,
such that the atomic masses will be included in the coordinate file.


CP2K
~~~~

See http://cp2k.berlios.de/

a. The optimization::

        &GLOBAL
          RUN_TYPE GEO_OPT
          PRINT_LEVEL LOW
        &END GLOBAL

b. The calculation of the energy (one SCF optimization)::

        &GLOBAL
          RUN_TYPE  ENERGY_FORCE
          PRINT_LEVEL HIGH
        &END GLOBAL

c. The frequency calculation::

        &VIBRATIONAL_ANALYSIS
          FULLY_PERIODIC T
        &END VIBRATIONAL_ANALYSIS
        ...
        &GLOBAL
          RUN_TYPE  vibrational_analysis
          PRINT_LEVEL medium
        &END GLOBAL


Gaussian 03
~~~~~~~~~~~

See http://www.gaussian.com/

In order to have the Hessian being printed, the following specification is
required for *large* molecules::

    iop(7/33=1) iop(6/7=3)



Q-Chem
~~~~~~

See http://www.q-chem.com/

(TODO)

Q-Chem/CHARMM (QM/MM)
~~~~~~~~~~~~~~~~~~~~~

(TODO)
