Chemical Physics with TAMkin
============================


Macroscopic properties of molecular systems
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TAMkin describes the partition function of molecular systems starting from the
the ideal gas and the harmonic oscillator approximation. Once the partition
function is specified, all thermodynamic quantities can be derived as a function
of temperature and pressure (or density). Several extensions to the ideal gas
and harmonic oscillator approximation are supported to improve the agreement
between computational predictions and experimental observations.


Definition of the partition function
------------------------------------

In order to compute thermodynamic quantities with TAMkin, it is crucial to
specify the correct parameters that define the partition function for the
molecular system under scrutiny. One must pay special attention to the
compatibility of the normal mode analysis (NMA, to compute the vibrational
spectrum) and the definition of the partition function. TAMkin will not complain
about `unusual` combinations of NMA's and Parition functions.

In the tutorial below, the microscopic computations are all carried out with
Gaussian03, but one can replace the ``load_molecule_...`` line with anything
suitable to load data from other quantum chemistry packages. The same
instructions also work starting from Gaussian09 computations.

The text below shows `snippets` of Python code to perform certain tasks. When
you use these `snippets` in a TAMkin script, do not forget to add the line
``from tamkin import *`` on top of the file.


Ideal gas molecules
^^^^^^^^^^^^^^^^^^^

The following three lines of Python code construct a partition function for
an ideal gas molecule using the standard options::

    mol = load_molecule_g03fchk("gaussian.fchk")
    nma = NMA(mol, ConstrainExt())
    pf = PartFun(nma, [ExtTrans(), ExtRot()])

In this example, the file ``"gaussian.fchk"`` comes from a Gaussian frequency
computation on the optimized geometry of the molecule of interest.

The following options can be used to specify the details, but the standard
values for these options are suitable for most applications.

* **Options for ConstrainExt**. See :class:`tamkin.nma.ConstrainExt` for the
  details. The ``ConstrainExt()`` object specifies that the normal mode analysis
  must be performed in 3N-6 (or 3N-5) internal coordinates. This is the default
  procedure in most programs.

  The method only gives the correct result when the geometry is sufficiently
  optimized, i.e. when the gradient of the energy is approximately zero. This
  condition will be checked by comparing the maximum absolute value of the
  gradient vector with a threshold value. If the trheshold is violated, TAMkin
  raises an error and your script will abort. The threshold can be specified as
  an option to ``ConstrainExt``. The following would relax the default
  threshold::

      nma = NMA(mol, ConstrainExt(gradient_threshold=0.01))

  The method also checks if the molecule is linear, which is based on the
  angular moments of inertia. Any angular momentum below the parameter
  im_threshold is treated as zero. The following setting would allow more
  deviations from linearity in linear molecules::

      nma = NMA(mol, ConstrainExt(im_threshold=10.0))

  *Note:* that the im_threshold value is given here in internal (atomic) units.
  For comparison, the lowest moment of inertia in ethane is about 40435 atomic
  units.

  One can always check the log file (see below) of a partition function to see
  if the molecule was considered to be linear or not. When combining both
  options, they must be separated by a comma::

      nma = NMA(mol, ConstrainExt(gradient_threshold=0.01, im_threshold=10.0))

* **Arguments for the PartFun object**. See :class:`tamkin.partf.PartFun` for
  the details. The partition function object is merely a definition of the
  partition function. It can be used to compute the numerical value of the
  partition function (and many other properties) at a given temperature and
  pressure.

  The contributions to the partition function must be specified when the object
  is created. By default, the vibrational and the electronic part are included.
  The first argument is a normal mode analysis object. It is used to set up the
  vibrational and the electronic contribution. Additional contributions are
  givin in a list (square brackets). For gas phase molecules, one typically adds
  the ``ExtTrans`` and ``ExtRot`` contributions.

* **Options for ExtTrans**. See :class:`tamkin.partf.ExtTrans` for the details.
  The ``ExtTrans`` class defines the translational contribution to the partition
  function. For practical reasons, this contribution also includes the many
  body effects of the partition function of the classical (ideal) gas, and the
  optional contributions for the constant pressure ensemble. The following
  options can be used:

  ``cp`` (default ``cp=True``).
    By default the translational contribution includes the corrections to
    describe an NpT partition function instead of an NVT partition function. If
    you are interested in the constant volume boundary conditions, use the
    option ``cp=False``

  ``gaslaw`` (default ``gaslaw=IdealGasLaw()``).
    The ideal gas law is the only gas law implemented so far, and is therefore
    also used by default. One may specify some options for the IdealGasLaw
    object. They are discussed below.

  ``dim`` (default ``dim=3``).
    The dimension of the gas. For ordinary gases, the dimension is three.

  ``mobile`` (default ``mobile=None``).
    One can optionally specify that only a part of the system is translational
    freedom. This is not relevant for molecules in the gas phase.

* **Options for ExtRot**. See :class:`tamkin.partf.ExtRot` for the details.

  ``symmetry_number`` (default ``symmetry_number=None``).
    When the symmetry number is not given, it is computed from the molecular
    geometry and topology. This may not work properly or very slowly for
    gigantic systems. In that case, specify symmetry_number=1, or whatever the
    number it should be.

  ``im_threshold`` (default ``im_threshold=1.0``).
    The threshold to determine if the molecule is linear or not. If one of the
    moments of inertia drops below this number, the molecule is considered to be
    linear. The value 1.0 is in internal (atomic) units.

* **Options for IdealGasLaw**. See :class:`tamkin.partf.IdealGasLaw` for the
  details. The ideal gas law has two optional parameters.

  ``pressure`` (default ``pressure=None``).
    The default value of ``pressure`` is 1 bar for 3D gases, 4.86e-05 atomic
    units for 2D gases (surface tension of water) and 1.0 atomic units for any
    other dimension. *Note:* several quantities derived from the partition
    function do not explicitly depend on the pressure in the case of ideal
    gases. In case you want to see the pressure dependence, use the method
    ``ExtTrans.set_pressure()`` and compute the thermodynamic quantities
    afterwards.

  ``dim`` (default ``dim=3``).
    The dimension of the gas. This must match the option ``dim`` given to
    ``ExtTrans``. When the ideal gas law is not specified in ExtTrans, the
    default value will have automatically the proper dimension.


Immobile adsorbed molecules
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Make sure you first read and understand the section on partition functions for
ideal gas molecules.

In this section, we show how one defines a partition function for a particle
that is adsorbed on a surface (flat or inside a porous material) and that it can
not rotate or displace over the surface once adsorbed. If it has to adsorb at
another place, or somewhere else, it first has to desorb and adsorb again.

We assume that the adsorption energy is computed with Gaussian using a cluster
approximation for the surface. This means that some the cluster is terminated
and that the atoms at the termination are fixed in space with constraints during
the geometry optimization. We also assume that the adsorbed molecule is free to
rotate as it can do in the gas phase.

The following code can be used to define the partition function for such a
system::

    fixed = [0, 1, 2, ...] # atom indexes of the fixed atoms, counting from zero
    mol_both = load_molecule_g03fchk("gaussian_both.fchk")
    nma_both = NMA(mol_both, PHVA(fixed))
    pf_both = PartFun(nma_both, [])

Compared to the gas phase, external translation and rotation are removed. Yhe
file ``"gaussian_both.fchk"`` comes from a frequency computation of the adsorbed
molecule on the cluster model of the surface.

The partition function of the surface without absorbed species is defined as
follows::

    fixed = [0, 1, 2, ...] # atom indexes of the fixed atoms, counting from zero
    mol_surf = load_molecule_g03fchk("gaussian_surf.fchk")
    nma_surf = NMA(mol_surf, PHVA(fixed))
    pf_surf = PartFun(nma_surf, [])

The surface is treated as a cluster fixed in space, i.e. there are not external
rotation and translation contributions to its partition function. The file
``"gaussian_surf.fchk"`` comes from a frequency computation on the surface
cluster model. The geometry of the cluster must be optimized with constraints on
the atoms that terminate the cluster.

One may load the indexes of the fixed atoms from a Gaussian ``.com`` file as
follows::

    fixed = load_fixed_g03com("gaussian.com")

Be aware that the fixed atom indexes may be different in the two computations,
but we recommend some consistency in this context. The following convention
avoids a lot of confusion: put all your surface atoms in the beginning of the
geometry definition, and within this group of atoms, put all fixed atoms first,
then the free atoms.

Mobile adsorbed molecules
^^^^^^^^^^^^^^^^^^^^^^^^^

Make sure you first read and understand the section on partition functions for
ideal gas molecules.

In this section, we show how one defines a partition function for a particle
that is adsorbed on a surface. We assume that the particle can still hover over
the surface and that this translational motion can be modeled with a 2D ideal
gas partition function with a constant surface area.

Further we assume that the adsorption energy is computed with Gaussian using
a cluster approximation for the surface. This means that some the cluster is
terminated and that the atoms at the termination are fixed in space with
constraints during the geometry optimization. We also assume that the adsorbed
molecule is free to rotate as it can do in the gas phase.

The following code can be used to define the partition function for such a system::

    fixed = [0, 1, 2, ...] # atom indexes of the fixed atoms, counting from zero
    mobile = [5, 6, 7, ...] # atom indexes of the mobile atoms, counting from zero
    mol_both = load_molecule_g03fchk("gaussian_both.fchk")
    nma_both = NMA(mol_both, PHVA(fixed))
    pf_both = PartFun(nma_both, [ExtTrans(cp=False, dim=2, mobile=mobile), ExtRot()])

In this code, the file ``"gaussian_both.fchk"`` comes from a frequency
computation of the adsorbed molecule on the cluster model of the surface. The
partition function of the surface without the adsorbed molecule is constructed
as follows::

    fixed = [0, 1, 2, ...] # atom indexes of the fixed atoms, counting from zero
    mol_surf = load_molecule_g03fchk("gaussian_surf.fchk")
    nma_surf = NMA(mol_surf, PHVA(fixed))
    pf_surf = PartFun(nma_surf, [])

The surface is treated as a cluster fixed in space, i.e. there are not external
rotation and translation contributions to its partition function. The file
``"gaussian_surf.fchk"`` comes from a frequency computation on the surface
cluster model. The geometry of the cluster must be optimized with constraints on
the atoms that terminate the cluster.

One may load the indexes of the fixed atoms from a Gaussian ``.com`` file as
follows::

    fixed = load_fixed_g03com("gaussian.com")

Be aware that the fixed atom indexes may be different in the two computations,
but we recommend some consistency in this context. The following convention
avoids a lot of confusion: put all your surface atoms in the beginning of the
geometry definition, and within this group of atoms, put all fixed atoms first,
then the free atoms.


Free or hindered internal rotors
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Make sure you first read and understand the section on partition functions for
ideal gas molecules.

TODO


The Partition function dump file
--------------------------------

After a partition function is defined in your script, one can write the entire
description to a text file for later reference::

    pf.write_to_file("partfun.txt")

It is recommended to double check the contents of the file.


Computation of thermodynamic quantities
---------------------------------------

Once the partition function of a system is defined, one can start computing
thermodynamic quantities at different temperatures and pressures (or densities).


Overview of standard quantities
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Thermodynamic quantities can be computed for a given ``PartFun`` object by calling
the appropriate methods. All extensive quantities, i.e. all quantities except
the chemical potential, are transformed into intensive quantities by dividing
through the number of particles. The following table relates the methods to the
meaning of the returned numbers for two common ensembles.

========================= ====================== ====================================================== ====================================================
``PartFun`` method        Internal unit          NVT Ensemble (3D gas)                                  NpT Ensemble (3D gas)
========================= ====================== ====================================================== ====================================================
``internal_energy``       Hartree/particle       Internal energy (per particle)                         Enthalpy (per particle)
``heat_capacity``         Hartree/(K*particle)   Heat capacity at constant volume (per particle)        Heat capacity at constant pressure (per particle)
``free_energy``           Hartree/particle       Helmholtz free energy (per particle)                   Gibbs free energy (per particle)
``chemical_potential``    Hartree/particle       Chemical potential                                     (idem)
``entropy``               Hartree/particle       Entropy (per particle)                                 (idem)
``log``                   1/particle             Logarithm of the partition function (per particles)    (idem)
``logt``                  1/(K*particle)         First derivative of ``log`` towards temperature        (idem)
``logtt``                 1/(K^2*particle)       Second derivative of ``log`` towards temperature       (idem)
========================= ====================== ====================================================== ====================================================

One can print out these values in a TAMkin script::

    from molmod import *  # for the unit conversion
    pf = ...
    print "The internal energy at 300K [kJ/mol]", pf.internal_energy(300)/kjmol
    print "The heat capacity at 300K [J/mol/K]", pf.heat_capacity(300)/(joule/(mol*kelvin))

Unit conventions and reference values
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Four thermodynamic functions in the table above have a poorly defined reference
value. The problematic cases are ``free_energy``, ``chemical_potential``,
``entropy``, ``log``. The translational contribution to these quantities
contains a term that is proportional to

.. math:: \ln\left(\frac{V}{N}\right).

In principle one can only define the logarithm of a dimensionless number.
Because the program works in some unit system, it actually computes

.. math:: \ln\left(\frac{V}{NV_0}\right),

where :math:`V_0` is the internal unit of volume. From the numerical
perspective, the contribution due to the unit can be considered as a separate
factor,

.. math:: \ln\left(\frac{V}{N}\right) - \ln(V_0),

which reveals that this unit convention affects the reference value of the four
quantities. This reference `correction` may in practice be temperature
dependent. For example, the free energy per particle is defined as:

.. math:: F_1 = -k_BT\frac{\ln(Z_N)}{N},

where :math:`Z_N` is the many body partition function, and :math:`N` is the
number of particles. This gives a temperature dependent reference correction:

.. math:: k_BT\ln(V_0).

It is clear that genuinely physical, i.e. measurable, quantities should not have
this weakness. This issue will return in the computation of equilibrium and rate
constants.

Poking under the hood
^^^^^^^^^^^^^^^^^^^^^

Besides the standard thermodynamic functions, all internal quantities of the
partition function and its contributions are also accessible. For example, one
computes the translational contribution to the free energy as follows::

    from molmod import *  # for the unit conversion
    pf = ...
    print "The free energy at 300K due to translation [kJ/mol]", pf.translational.internal_energy(300)/kjmol

A complete overview of internals can be found in the reference documentation
of the :mod:`tamkin.partf` module, or by reading the source code.


Generating tables
^^^^^^^^^^^^^^^^^

Tables of thermodynamic quantities can be computed for given temperatures and
sorting out all contributions from the components of the partition function to
each quantity. The example below generates a CSV file that can be loaded into
spreadsheet software. ::

    from tamkin import *
    molecule = load_molecule_g03fchk("gaussian.fchk")
    nma = NMA(molecule, ConstrainExt())
    pf = PartFun(nma, [ExtTrans(), ExtRot()])
    ta = ThermoAnalysis(pf, [300, 400, 500, 600])
    ta.write_to_file("thermo.csv")


The CSV file contains tables with thermodynamic quantities, at the temperatures
in the second argument of the ThermoAnalysis constructor, corresponding to the
PartFun methods as explained the table below.

==================== ============ ==========================
Name in CSV file     Unit         ``PartFun`` method name
==================== ============ ==========================
Energy               kJ/mol       ``internal_energy``
Heat capacity        J/(mol*K)    ``heat_capacity``
Free energy          kJ/mol       ``free_energy``
Chemical potential   kJ/mol       ``chemical_potential``
Entropy              J/(mol*K)    ``entropy``
log(q)               1/mol        ``log``
d log(q) / dT        1/(mol*K)    ``logt``
d^2 log(q) / dT^2    1/(mol*K^2)  ``logtt``
==================== ============ ==========================


Thermodynamic equilibrium
~~~~~~~~~~~~~~~~~~~~~~~~~

Definition of the equilibrium constant
--------------------------------------

Unit conventions
----------------

Computation of the equilibrium constant
---------------------------------------

ThermodynamicModel objects
--------------------------

Reaction kinetics
~~~~~~~~~~~~~~~~~

Definition of the equilibrium constant
--------------------------------------

Unit conventions
----------------

Computation of the equilibrium constant
---------------------------------------

KineticModel objects
--------------------

Tunneling corrections
^^^^^^^^^^^^^^^^^^^^^

ReactionAnalysis objects -- fitting kinetic parameters A and E\ :sub:`a`
------------------------------------------------------------------------

Error analysis
^^^^^^^^^^^^^^
