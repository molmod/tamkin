
The TAMkin driver script
########################

The TAMkin driver script (``tamkin-driver``) runs a standardized TAMkin
computation based on the files and directories present in the current directory.
This script needs no command line arguments. (When command line arguments are
present, this documentation is printed on screen.)


Assumed layout of the current directory
=======================================

This script assumes that the current directory has subdirectories as follows::

    re_*/            # Reactant molecule subdirectories (at least one).
    ts_*/            # Transition state molecule subdirectory (optional, at most one).
    pr_*/            # Reaction product molecule subdirectories (optional).
    kinetics.cfg     # Contains parameters for reaction kinetics computation.
    equilibrium.cfg  # Contains parameters for chemical equilibrium computation.

The layout of each `molecule` subdirectory is documented below.

Reactants are mandatory. One transition state or one or more products may
be present. Depending on the available directories, the following two
computations may take place:

* When a transition state is present, the kinetic parameters are computed. In
  this case, a config file ``kinetics.cfg`` must be present. (Details given below.)

* When reaction products are present, the equilibrium constant is computed. In
  this case, one may add a file ``equilibrium.cfg``. (Details given below.)

Even when no transition state or reaction products are present, useful computation can be
performed, e.g. by providing a file ``thermo.cfg`` telling at which temperatures all
properties of the reactants must be computed. (Details given below.)


Assumed layout of a molecule directory (re_*/, ts_*/ or pr_*/)
==============================================================

At most one of the following sets of files must be present in the directory::

    freq/gaussian.fchk     # the formatted checkpoint file of a Gaussian03/09 computation.

    freq/POSCAR            # input geometry and frequency output of a VASP calculation
    freq/OUTCAR

    sp/cp2k.out            # CP2K output files for frequency and single-point calculations.
    freq/cp2k.out          # The single-point must contain energy and forces.

The following files may be present in the directory::

    molecule.cfg           # specifies additional parameters of the molecule.
                           # (Details given below.)

    dftd3/dftd3.out        # the screen output of the dftd3 program redirected
                           # to a file. this adds an a posteriori Grimme
                           # Dispersion correction to the energy.

    sp/gaussian.fchk       # a Gaussian03/09 formatted checkpoint file to
                           # replace the energy from the frequency job by a more
                           # accurate one.

    sp/OUTCAR              # The output of a VASP single-point calculation with a
                           # more accurate energy compared to the frequency calculation.

    rotor_g_*/gaussian.log # adds a rotor based on a Gaussian 03/09 relaxed
    rotor_g_*/rotor.cfg    # scan. (Details given below. More than one allowed.)

    rotor_f_*/rotor.cfg    # specifies a free rotor. (Details given
                           # below. More than one allowed.)

    rotor_c_*/rotor.dat    # specifies a custom hindered rotor. (Details given
    rotor_c_*/rotor.cfg    # below. More than one allowed.)

All other files are simply ignored.


Config and data files
=====================

All configuration files (``*.cfg``) have the following format. Every non-empty
line consists of a key followed one or more values, all separated by whitespace.
Comments can be added with a #, just as in Python source code. Each file has its
specific keys that are processed. Unknown keys are ignored.

For example, the following is a valid ``molecule.cfg`` file::

    freq_scaling 0.95
    symnum 3
    # This line is ignored.
    unkown_option is ignored as well

The following config files are read by the ``tamkin-driver`` script:

* **kinetics.cfg**: ``temp_low``, ``temp_high``, ``temp_step``, ``tunneling``

    ``temp_low`` and ``temp_high`` (mandatory)
        These specifiy the minimum and maximum temperature for the Arrhenius
        plot.

    ``temp_step`` (optional)
        The temperature interval for the datapoints for the Arrhenius plot.
        [default=10]

    ``tunneling`` (optional)
        When set to True, the Eckart tunneling is applied to the computation
        of reaction rates.

* **equilibrium.cfg**: ``temps``

    ``temps`` (optional)
        A list of temperatures at which the equilibrium constant is computed.

* **thermo.cfg**: ``temps``

    ``temps`` (optional)
        A list of temperatures at which all properties of the partition functions must be
        computed.

* **molecule.cfg**: ``symnum``, ``freq_scaling``, ``zp_scaling``, ``periodic``

    ``freq_scaling`` (optional)
        The frequency scaling factor to correct for systematic deviations of the
        level theory used to compute the Hessian. [default=1.0]

    ``symnum`` (optional)
        This keyword can be used to assign the rotational symmetry number. For
        molecules with less than 10 atoms, this number is estimated
        automatically when not given. For larger molecules, the default value is
        1.

    ``zp_scaling`` (optional)
        The zero-point scaling factor to correct for systematic deviations of the
        level theory used to compute the Hessian. [default=1.0]

    ``periodic`` (optional)
        The default value is True for VASP calculations. It is ignored for Gaussian
        calculations.

* **rotor_g_*/rotor.cfg**: ``dofmax``, ``even``, ``fortran``, ``num_levels``,
  ``rotsym``, ``top``
* **rotor_f_*/rotor.cfg**: ``dihed``, ``dofmax``, ``fortran``, ``num_levels``,
  ``rotsym``, ``top``
* **rotor_c_*/rotor.cfg**: ``even``, ``dihed``, ``dofmax``, ``fortran``,
  ``num_levels``, ``rotsym``, ``top``

    ``even`` (optional)
        A boolean (True or False) to indicate that the torsional potential is
        even. [default=False]

    ``dihed`` (mandatory)
        A list of four atom indexes that define the dihedral angle, separated by
        whitespace.

    ``dofmax`` (optional)
        The maximum number of cosines used to represent the torsional potential.
        if the potential is not even, the same number of sines is also used.
        [default=5]

    ``fortran`` (optional)
        A boolean (True or False) to indicate that the atom indexes are given in
        Fortran convention. (Counting starts from one instead of zero). This is
        option relevant for the keys ``dihed`` and ``top``. [default=False]

    ``num_levels`` (optional)
        The number of energy levels considered in the QM treatment of the rotor.
        [default=50]

    ``rotsym`` (optional)
        The rotational symmetry of the internal rotor. [default=1]

    ``top`` (optional)
        The atoms in the rotating top. When not given, an attempt is made to
        derive this top from the choice of the dihedral angle and the molecular
        topology. (This attempt is often not successful for structures
        containing multiple molecules. In that case, top must be
        provided.

* **rotor_c_*/rotor.dat**

    The file ``rotor_c_*/rotor.dat`` just contains two columns of data, angles
    (radians) and energies (hartree), that specify the custom torsional potential.
    It does not follow the ``*.cfg`` format.

* **fixed.txt**

    A file with atomic indices, which are to be considered fixed in space. When fixed
    atoms are present, the PHVA method is used and external degrees of freedom are no
    longer projected out. The format of files with atomic indices is discussed below.

* **blocks.txt**

    A file with groups of atoms that should be treated as rigid blocks in the normal-mode
    analysis. The format of files with atomic indices is discussed below.


Notes
=====

* The energy levels of a hindered rotor are found by solving the Schroedinger
  equation in a plane wave basis. A truncated Fourier series is used to expand
  the potential energy. The truncation can be controlled with the ``dofmax``
  parameter. When the RMSD between the Fourier series and the data is larger
  than 1 kJ/mol or 1%, the driver will stop with an error message. The simplest
  solution is to increase ``dofmax`` (above the default of 5). However, one also
  has to make sure that the potential from the relaxed scan is sensible. If it
  contains a rotational symmetry, limit the scan to one period and add the
  appropriate ``rotsym`` keyword in the ``rotor.cfg`` file. If the scan is even,
  one can again halve the range of the scan and add ``even true`` to the file
  ``rotor.cfg``. For example, for a standard methyl top, the scan of the
  dihedral angle must be limited to the interval [0 deg, 60 deg] and the following
  lines must be added to the file ``rotor.cfg`` ::

    rotsym 3
    even true

* Format of files with atomic indices:

  - Comments start with #. All text after this character is ignored (except for shift, see
    below).
  - Empty lines are ignored.
  - A non-empty line contains a series of integers or ranges defined by two integers. A
    range is defined as either [N1,N2] or [N1,N2[. In the former case, the N2 is included
    in the range while in the latter case it is not.
  - If a line is present of the form "#shift=N" (without spaces), N is used to shift all
    indices upon reading. When it is not encountered, N is assumed to be -1. When N=0,
    the indices are C-style (counting starts from zero). When N=-1, the indices are
    Fortran-style (counting starts from 1).
  - When defining multiple rigid blocks, the atom indices must be grouped in paragraphs
    separated by an empty line.
