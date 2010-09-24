Chemical Physics -- Working with TAMkin
=======================================


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
about `unusual` combinations of NMA's and Partition functions.

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
  gradient vector with a threshold value. If the threshold is violated, TAMkin
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

* **Full instead of ConstrainExt**. See :class:`tamkin.nma.Full` for the
  details. Instead of constraining the external degrees of freedom, one may
  also perform the normal mode analysis in 3N coordinates and hope that the
  Hessian is accurate enough to produce 6 zero frequencies. They will be
  elminated automatically for the vibrational contribution to the partition
  function.

  For a proper identification of the zero frequencies, the ``Full`` treatment
  will also check the linearity of the molecule, just like the ``ConstrainExt``
  treatment. One can relax the linearity check as follows::

      nma = NMA(mol, Full(im_threshold=10.0))

  There are no other options for the full treatment.

* **Arguments for the PartFun object**. See :class:`tamkin.partf.PartFun` for
  the details. The partition function object is merely a definition of the
  partition function. It can be used to compute the numerical value of the
  partition function (and many other properties) at a given temperature and
  pressure.

  The contributions to the partition function must be specified when the object
  is created. By default, the vibrational and the electronic part are included.
  The first argument is a normal mode analysis object. It is used to set up the
  vibrational and the electronic contribution. Additional contributions are
  given in a list (square brackets). For gas phase molecules, one typically adds
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
    you are interested in the constant volume boundary conditions (NVY
    ensemble), use the option ``cp=False``.

  ``gaslaw`` (default ``gaslaw=IdealGasLaw()``).
    The ideal gas law is the only gas law implemented so far, and is therefore
    also used by default. One may specify some options for the IdealGasLaw
    object. They are discussed below.

  ``dim`` (default ``dim=3``).
    The dimension of the gas. For ordinary gases, the dimension is three.

  ``mobile`` (default ``mobile=None``).
    One can optionally specify that only a part of the system is translational
    freedom. This is not relevant for molecules in the gas phase.

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

* **Options for Vibrations**. See :class:`tamkin.partf.Vibrations` for the
  details. One does not have to add a ``Vibrations`` object to the list of
  contributions in the partition function, unless one wants to modify the
  options of the vibrational contribution.

  ``classical`` (default ``classical=False``)
    When True, the vibrations are treated classically. The QM treatment is the
    default. This may be useful to compare TAMkin results with Monte Carlo or
    Molecular Dynamics simulations.

  ``freq_scaling`` (default ``freq_scaling=1.0``)
    Scaling factor for the classical part of the vibrational partition function,
    ie. excluding the zero-point term.

  ``zp_scaling``  (default ``zp_scaling=1.0``)
    Scaling factor for the zero-point term in the vibrational contribution [default=1]


Immobile adsorbed molecules
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Make sure you first read and understand the section on partition functions for
ideal gas molecules.

In this section, we show how one defines a partition function for a particle
that is adsorbed on a surface (flat or inside a porous material) and that it can
not rotate or displace over the surface once adsorbed. If it has to adsorb at
another place, or somewhere else, it first has to desorb and adsorb again.

We assume that the adsorption energy is computed with Gaussian using a cluster
approximation for the surface. This means that the cluster is terminated
and that the atoms at the termination are fixed in space with constraints during
the geometry optimization. We also assume that the adsorbed molecule is free to
rotate as it can do in the gas phase.

The following code can be used to define the partition function for such a
system::

    fixed = [0, 1, 2, ...] # atom indexes of the fixed atoms, counting from zero
    mol_both = load_molecule_g03fchk("gaussian_both.fchk")
    nma_both = NMA(mol_both, PHVA(fixed))
    pf_both = PartFun(nma_both, [])

Compared to the gas phase, external translation and rotation are removed. The
file ``"gaussian_both.fchk"`` comes from a frequency computation of the adsorbed
molecule on the cluster model of the surface.

The partition function of the surface without absorbed species is defined as
follows::

    fixed = [0, 1, 2, ...] # atom indexes of the fixed atoms, counting from zero
    mol_surf = load_molecule_g03fchk("gaussian_surf.fchk")
    nma_surf = NMA(mol_surf, PHVA(fixed))
    pf_surf = PartFun(nma_surf, [])

The surface is treated as a cluster fixed in space, i.e. there are no external
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
a cluster approximation for the surface. This means that the cluster is
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
``log``                   1/particle             Logarithm of the partition function (per particle)     (idem)
``logt``                  1/(K*particle)         First derivative of ``log`` towards temperature        (idem)
``logtt``                 1/(K^2*particle)       Second derivative of ``log`` towards temperature       (idem)
========================= ====================== ====================================================== ====================================================

One can print out these values in a TAMkin script::

    from molmod import *  # for the unit conversion
    pf = ...
    print "The internal energy at 300K [kJ/mol]", pf.internal_energy(300)/kjmol
    print "The heat capacity at 300K [J/mol/K]", pf.heat_capacity(300)/(joule/(mol*kelvin))


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

=============================================================================== ============ ==========================
Name in CSV file                                                                Unit         ``PartFun`` method name
=============================================================================== ============ ==========================
Energy                                                                          kJ/mol       ``internal_energy``
Heat capacity                                                                   J/(mol*K)    ``heat_capacity``
Free energy                                                                     kJ/mol       ``free_energy``
Chemical potential                                                              kJ/mol       ``chemical_potential``
Entropy                                                                         J/(mol*K)    ``entropy``
log= :math:`\frac{log(Z_N)}{N}`                                                 1/mol        ``log``
logt= :math:`\frac{\partial}{\partial T}\left(\frac{log(Z_N)}{N}\right)`        1/(mol*K)    ``logt``
logtt= :math:`\frac{\partial^2}{\partial T^2}\left(\frac{log(Z_N)}{N}\right)`   1/(mol*K^2)  ``logtt``
=============================================================================== ============ ==========================


Thermodynamic equilibrium
~~~~~~~~~~~~~~~~~~~~~~~~~


The method ``PartFun.logv``
---------------------------

To guarantee the numerical stability of the results obtained with TAMkin,
logarithms of partition functions are computed in the ``PartFun`` object and
its contributions. These can be used to compute the logarithm of the equilibrium
constant:

.. math:: \ln(K_c(T)) = \nu_C\ln(Z'_C(1, \ldots)) + \nu_D\ln(Z'_D(1, \ldots))
                       -\nu_A\ln(Z'_A(1, \ldots)) - \nu_B\ln(Z'_B(1, \ldots))

The method ``PartFun.logv`` computes the quantity :math:`\ln(Z'_X(1, \ldots))`.
The same method can be found in all the contributions to the partition function.
For all contributions, except the translational one, the method ``logv`` and
``log`` are identical.


``ThermodynamicModel`` objects
------------------------------

Given a list of partition functions of reactants (``pfs_react``) and a list of
product partition functions (``pfs_prod``), one must first construct a
``ThermodynamicModel`` object. ::

    tm = ThermodynamicModel(pfs_react, pfs_prod)

Then the equilibrium constant is computed at a certain temperature, ``temp``, as
follows::

    kc = tm.equilibrium_constant(temp)

This function takes one optional argument: ``do_log``, which is by default
``False``. When set to True, the logarithm of the partition function is
returned. The name of the SI unit of the equilibrium constant and the
corresponding conversion factor are attributes of the ``ThermodymanicModel``
object. This can be used to facilitate the output of equilibrium constants. For
example::

    print "K_c at 350K [%s] = %.5e" % (tm.unit_name, kc/tm.unit)

Currently TAMkin only supports ideal gases for the translational contribution to
the partition function, which means that :math:`K_c` does not depend on the
pressure set in ``ExtTrans.gaslaw.pressure``.

The change in free energy is computed (and printed) as follows::

    drf = tm.free_energy_change(temp)
    print "Change in free energy [kJ/mol] = %.3f" % (drf/kjmol)

The change in free energy does depend on the pressure parameter in the
translational part of the partition functions of the reactants and products.


Reaction kinetics
~~~~~~~~~~~~~~~~~

KineticModel objects
--------------------

Tunneling corrections
^^^^^^^^^^^^^^^^^^^^^

ReactionAnalysis objects -- fitting kinetic parameters A and E\ :sub:`a`
------------------------------------------------------------------------

Error analysis
^^^^^^^^^^^^^^
