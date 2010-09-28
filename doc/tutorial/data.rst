Generating data for TAMkin
==========================

TAMkin can read information from CHARMM, Gaussian, GAMESS, CP2K, CPMD,
Q-CHEM and VASP. Some of these programs offer freedom in the amount of output
that is printed to the output file. To produce output that TAMkin can read from,
we suggest the following settings.

Although most packages can do the standard frequency computation, the actual
frequencies are never read by TAMkin. Instead, TAMkin reads all the ingredients
(such as the Hessian and the atomic masses) required for the computation the
frequencies. This makes it possible to compute frequencies with different types
of normal mode analysis methods.


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

CPMD
~~~~

**TODO**

GAMESS
~~~~~~

**TODO**

Gaussian 98/03/09
~~~~~~~~~~~~~~~~~

See http://www.gaussian.com/

Always let Gaussian keep the checkpoint file by adding the following line on top
of the Gaussian input::

    %chk=prefix.chk

Replace ``prefix`` by the prefix of the ``.com`` file. After the computation,
immediately transform the checkpoint file into a formatted checkpoint file::

    toon@poony ~> formchk prefix.chk prefix.fchk

The original checkpoint file is in binary format and is not transferable from
one CPU architecture to another, so the ``formchk`` command must be executed on
the (type of) machine where the computation was done. The formmated checkpoint
file is a text file that is used by TAMkin to extract all the results from a
Gaussian computation.

For optimizations we recommend the ``opt(tight)`` option as it improves the
reproducibility of the frequencies. For frequency calculations, one can use
``freq(noraman)`` to skip the computation of Raman intensities. They are not
used by TAMkin and require about 10% of the computation time in a frequency job.
It is not recommended to work with link jobs because the formatted checkpoint
file only contains information about the last link.


Q-Chem
~~~~~~

See http://www.q-chem.com/

**TODO**

Q-Chem/CHARMM (QM/MM)
~~~~~~~~~~~~~~~~~~~~~

**TODO**

VASP
~~~~

**TODO**


Loading data into TAMkin
========================

Once the proper output is created with a computational chemistry package, it
can be loaded into TAMkin to perform a normal mode analysis and thermochemistry
computations. Detailed instructions per package are listed below.

CHARMM
~~~~~~

For the details, see :mod:`tamkin.io.charmm`.

**TODO**

CP2K
~~~~

For the details, see :mod:`tamkin.io.cp2k`.

**TODO**

CPMD
~~~~

For the details, see :mod:`tamkin.io.cpmd`.

**TODO**

GAMESS
~~~~~~

For the details, see :mod:`tamkin.io.gamess`.

**TODO**

Gaussian 98
~~~~~~~~~~~

For the details, see :func:`tamkin.io.gaussian.load_molecule_g98fchk`.

Given a formatted checkpoint file, it is loaded as follows::

    molecule = load_molecule_g98fchk("freq.fchk")

where ``"freq.fchk"`` is the name of the formatted checkpoint file of a
frequency computation in Gaussian98. One may also provide a second formatted
checkpoint with a refined energy computation::

    molecule = load_molecule_g98fchk("freq.fchk", "ener.fchk")

It is also possible to give a numerical value for the refined energy (in
internal units, i.e. Hartree)::

    molecule = load_molecule_g98fchk("freq.fchk", energy=-135.12597)

Gaussian98 does not write the atomic masses to the formatted checkpoint file.
Therefore the atomic masses used by Gaussian98 are added in a rather artificial
way inside the ``load_molecule_g98fchk`` routine. If you wish to override these
masses with the IUPAC 2005 values, use the following snippet::

    from molmod.periodci import periodic
    import numpy

    molecule = load_molecule_g98fchk("gaussian.fchk")
    new_masses = numpy.array([periodic[n].mass for n in molecule.numbers])
    molecule = molecule.copy_with(masses=new_masses)

Gaussian 03/09
~~~~~~~~~~~~~~

For the details, see :func:`tamkin.io.gaussian.load_molecule_g03fchk`.

`Note`: Formatted checkpoint files of Gaussian03 and Gaussian09 can both be
read with ``load_molecule_g03fchk``

Given a formatted checkpoint file, it is loaded as follows::

    molecule = load_molecule_g03fchk("freq.fchk")

where ``"freq.fchk"`` is the name of the formatted checkpoint file of a
frequency computation in Gaussian03 or Gaussian09. One may also provide a second
formatted checkpoint with a refined energy computation::

    molecule = load_molecule_g03fchk("freq.fchk", "ener.fchk")

It is also possible to give a numerical value for the refined energy (in
internal units, i.e. Hartree)::

    molecule = load_molecule_g03fchk("freq.fchk", energy=-135.12597)

Torsional potentials
--------------------

For the details, see :func:`tamkin.io.gaussian.load_rotscan_g03log`.

For the treatment of hindered internal rotors one must load the torsional
potential data from a relaxed potential energy surface scan. ::

    load_rotscan_g03log("scan.log")

where ``"scan.log"`` refers to the Gaussian log file of the PES scan. (The
formatted checkpoint file is not used here as it does not specify the dihedral
angle that was used for the scan.) The routine ``load_rotscan_g03log`` tries
to figure out which atoms belong to the rotor `top`. It will always take the
smallest of the two possibilities. In case this does not work for some reason,
one may manually specify the atom indexes of the rotor `top`::

    top_indexes = [0, 1, 3, 5] # Counting starts from zero.
    load_rotscan_g03log("scan.log", top_indexes)



Q-Chem
~~~~~~

For the details, see :mod:`tamkin.io.qchem`.

**TODO**

Q-Chem/CHARMM (QM/MM)
~~~~~~~~~~~~~~~~~~~~~

For the details, see :mod:`tamkin.io.qchem`.

**TODO**

VASP
~~~~

For the details, see :mod:`tamkin.io.vasp`.

**TODO**
