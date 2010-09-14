Generating data for TAMkin
==========================

TAMkin can read information from ADF, CHARMM, Gaussian, GAMESS, CP2K, CPMD,
Q-CHEM and VASP. Some of these programs offer freedom in the amount of output
that is printed to the output file. To produce output that TAMkin can read from,
we suggest the following settings.


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


Gaussian 03/09
~~~~~~~~~~~~~~

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


Q-Chem
~~~~~~

See http://www.q-chem.com/

(TODO)

Q-Chem/CHARMM (QM/MM)
~~~~~~~~~~~~~~~~~~~~~

(TODO)
