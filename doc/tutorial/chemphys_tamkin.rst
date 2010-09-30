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
specify the correct parameters for the partition function of the molecular
system under scrutiny. One must pay special attention to the compatibility of
the normal mode analysis (NMA, to compute the vibrational spectrum) and the
definition of the partition function. TAMkin will not complain about `unusual`
combinations of NMA's and Partition functions.

In the tutorial below, all the microscopic computations are all carried out with
Gaussian03, but one can replace the ``load_molecule_...`` line with anything
suitable to load data from other (quantum) chemistry packages. The same
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

  ``pressure`` (default ``pressure=None``, only allowed when ``cp==True``).
    The default value of ``pressure`` is 1 bar for 3D gases, 75.64 mNewton/m for
    2D gases (surface tension of water) and 1.0 atomic units for any other
    dimension. *Note:* several quantities derived from the partition function do
    not explicitly depend on the pressure in the case of ideal gases. If you
    want to see the (absence of) pressure dependence, use the method
    ``ExtTrans.set_pressure()`` and recompute the thermodynamic quantities
    afterwards.

  ``density`` (default ``density=None``, only allowed when ``cp==False``).
    The default value of ``density`` is 1 mol/m\ :sub:`dim`, where dim is the
    dimension of the gas. *Note:* several quantities derived from the partition
    function do not explicitly depend on the density in the case of ideal gases.
    If you want to see the (absence of) density dependence, use the method
    ``ExtTrans.set_density()`` and recompute the thermodynamic quantities
    afterwards.

  ``dim`` (default ``dim=3``).
    The dimension of the gas. For ordinary gases, the dimension is three.

  ``mobile`` (default ``mobile=None``).
    One may specify that only a part of the system has translational freedom.
    This is not relevant for molecules in the gas phase.

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
that is adsorbed on a surface (flat or inside a porous material) such that it
can not rotate or translate over the surface once adsorbed. It can only adsorb
at another site, if it first desorbs from the surface.

We assume that the adsorption energy is computed with Gaussian using a cluster
approximation for the surface. This means that the cluster is terminated
and that the atoms at the termination are fixed in space with constraints during
the geometry optimization.

The following code can be used to define the partition function for such a
system::

    fixed = [0, 1, 2, ...] # atom indexes of the fixed atoms, counting from zero
    mol_both = load_molecule_g03fchk("gaussian_both.fchk")
    nma_both = NMA(mol_both, PHVA(fixed))
    pf_both = PartFun(nma_both, [])

Compared to the gas phase, external translation and rotation are removed. The
file ``"gaussian_both.fchk"`` comes from a frequency computation of the adsorbed
molecule on the cluster model of the surface. The list ``fixed`` contains all
atom indexes that fixed in space during the optimization.

The partition function of the surface without absorbed species is defined as
follows::

    fixed = [0, 1, 2, ...] # atom indexes of the fixed atoms, counting from zero
    mol_surf = load_molecule_g03fchk("gaussian_surf.fchk")
    nma_surf = NMA(mol_surf, PHVA(fixed))
    pf_surf = PartFun(nma_surf, [])

The surface is treated as a cluster fixed in space, i.e. there are no external
rotation and translation contributions to the partition function. The file
``"gaussian_surf.fchk"`` comes from a frequency computation on the surface
cluster model. The geometry of the cluster must be optimized with constraints on
the atoms that terminate the cluster. The list ``fixed`` contains all
atom indexes that fixed in space during the optimization.

One may load the indexes of the fixed atoms from a Gaussian ``.com`` file as
follows::

    fixed = load_fixed_g03com("gaussian.com")

Be aware that the fixed atom indexes may be different in the two computations,
but we recommend some consistency in this context. The following convention
avoids a lot of confusion: put all your surface atoms in the beginning of the
geometry definition, and within this group of surface atoms, put all fixed atoms
first, then the free atoms.

Mobile adsorbed molecules
^^^^^^^^^^^^^^^^^^^^^^^^^

Make sure you first read and understand the section on partition functions for
ideal gas molecules.

In this section, we show how one defines a partition function for a particle
that is adsorbed on a surface. We assume that the particle can still hover over
the surface and that this translational motion can be modeled with a 2D ideal
gas partition function with a constant surface area. We also assume that the
adsorbed molecule is free to rotate as it can do in the gas phase.

Further we assume that the adsorption energy is computed with Gaussian using
a cluster approximation for the surface. This means that the cluster is
terminated and that the atoms at the termination are fixed in space with
constraints during the geometry optimization.

The following code can be used to define the partition function for such a system::

    fixed = [0, 1, 2, ...] # atom indexes of the fixed atoms, counting from zero
    mobile = [5, 6, 7, ...] # atom indexes of the mobile atoms, counting from zero
    mol_both = load_molecule_g03fchk("gaussian_both.fchk")
    nma_both = NMA(mol_both, PHVA(fixed))
    pf_both = PartFun(nma_both, [ExtTrans(cp=False, dim=2, mobile=mobile), ExtRot()])

In this code, the file ``"gaussian_both.fchk"`` comes from a frequency
computation of the adsorbed molecule on the cluster model of the surface. The
list ``fixed`` contains all atom indexes that fixed in space during the
optimization. The ``mobile`` atoms refer to those of the adsorbed species that
has 2D translational motion on the surface.

The partition function of the surface without the adsorbed molecule is
constructed as follows::

    fixed = [0, 1, 2, ...] # atom indexes of the fixed atoms, counting from zero
    mol_surf = load_molecule_g03fchk("gaussian_surf.fchk")
    nma_surf = NMA(mol_surf, PHVA(fixed))
    pf_surf = PartFun(nma_surf, [])

The surface is treated as a cluster fixed in space, i.e. there are not external
rotation and translation contributions to its partition function. The file
``"gaussian_surf.fchk"`` comes from a frequency computation on the surface
cluster model. The geometry of the cluster must be optimized with constraints on
the atoms that terminate the cluster. The list ``fixed`` contains all atom
indexes that fixed in space during the optimization.

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

Once a partition function is defined in your script, one can write an extensive
description to a text file for later reference::

    pf.write_to_file("partfun.txt")

It is recommended to double check the contents of the file.


Computation of thermodynamic quantities
---------------------------------------

Once a partition function object is created, one can start computing
thermodynamic quantities at different temperatures and pressures (or densities).


Overview of standard quantities
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Thermodynamic quantities can be computed for a given ``PartFun`` object by
calling the appropriate methods. All extensive quantities, i.e. all quantities
except the chemical potential, are transformed into intensive quantities by
dividing through the number of particles. The following table relates the
methods to the meaning of the returned numbers for two common ensembles.

========================= ====================== ======================================================================================== ====================================================
``PartFun`` method        Internal unit          NVT Ensemble (3D gas)                                                                    NpT Ensemble (3D gas)
========================= ====================== ======================================================================================== ====================================================
``internal_heat``         Hartree/particle       Internal energy (per particle)                                                           Enthalpy (per particle)
``heat_capacity``         Hartree/(K*particle)   Heat capacity at constant volume (per particle)                                          Heat capacity at constant pressure (per particle)
``free_energy``           Hartree/particle       Helmholtz free energy (per particle)                                                     Gibbs free energy (per particle)
``chemical_potential``    Hartree/particle       Chemical potential                                                                       (idem)
``entropy``               Hartree/particle       Entropy (per particle)                                                                   (idem)
``log``                   1/particle             Logarithm of the partition function (per particle)                                       (idem)
``logt``                  1/(K*particle)         First derivative of ``log`` towards temperature                                          (idem)
``logtt``                 1/(K^2*particle)       Second derivative of ``log`` towards temperature                                         (idem)
``logn``                  1/particle             Derivative of the logarithm of the partition function towards the number of particles    (idem)
``logv``                  1/particle             ``logn`` - :math:`\ln(V/N)`                                                              (idem)
========================= ====================== ======================================================================================== ====================================================

One can print any of these quantities in a TAMkin script::

    from molmod import *  # for the unit conversion
    pf = ...
    print "The internal heat at 300K [kJ/mol]", pf.internal_heat(300)/kjmol
    print "The heat capacity at 300K [J/mol/K]", pf.heat_capacity(300)/(joule/(mol*kelvin))


Poking under the hood
^^^^^^^^^^^^^^^^^^^^^

Besides the standard thermodynamic functions, all internal quantities of the
partition function and its contributions are also accessible. For example, one
computes the translational contribution to the free energy as follows::

    from molmod import *  # for the unit conversion
    pf = ...
    print "The free energy at 300K due to translation [kJ/mol]", pf.translational.internal_heat(300)/kjmol

A complete overview of internals can be found in the reference documentation
of the :mod:`tamkin.partf` module, or by reading the source code.


Generating tables
^^^^^^^^^^^^^^^^^

See :class:`tamkin.pftools.ThermoAnalysis` for the details.

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
in the second argument of the ``ThermoAnalysis`` constructor, corresponding to
the ``PartFun`` methods as explained the table below.

=============================================================================== ===================== ==========================
Name in CSV file                                                                Unit                  ``PartFun`` method name
=============================================================================== ===================== ==========================
Internal Heat                                                                   kJ/mol                ``internal_heat``
Heat capacity                                                                   J/(mol*K)             ``heat_capacity``
Free energy                                                                     kJ/mol                ``free_energy``
Chemical potential                                                              kJ/mol                ``chemical_potential``
Entropy                                                                         J/(mol*K)             ``entropy``
log= :math:`\frac{log(Z_N)}{N}`                                                 1                     ``log``
logt= :math:`\frac{\partial}{\partial T}\left(\frac{log(Z_N)}{N}\right)`        1/K                   ``logt``
logtt= :math:`\frac{\partial^2}{\partial T^2}\left(\frac{log(Z_N)}{N}\right)`   1/K^2                 ``logtt``
logn= :math:`\frac{\partial log(Z_N)}{\partial N}`                              1                     ``logn``
logv= :math:`\frac{\partial log(Z_N)}{\partial N} - \ln(V/N)`                   1 or ln(bohr^-dim)    ``logv``
=============================================================================== ===================== ==========================


Thermodynamic equilibrium
~~~~~~~~~~~~~~~~~~~~~~~~~

See :class:`tamkin.chemmod.ThermodynamicModel` for the details.

Given a list of partition functions of reactants (``pfs_react``) and a list of
product partition functions (``pfs_prod``), one can construct a
``ThermodynamicModel`` object. ::

    tm = ThermodynamicModel(pfs_react, pfs_prod)

With a ``ThermodynamicModel`` object one can compute the following quantities:

* **The equilibrium constant** is computed at a given temperature, ``temp``, as
  follows::

      kc = tm.equilibrium_constant(temp)

  This function takes one optional argument: ``do_log``, which is by default
  ``False``. When set to True, the logarithm of the equilibrium constant is
  returned. The name of the SI unit of the equilibrium constant and the
  corresponding conversion factor are attributes of the ``ThermodymanicModel``
  object. This can be used to facilitate the output of equilibrium constants.
  For example::

      print "K_c at 350K [%s] = %.5e" % (tm.unit_name, kc/tm.unit)

  Currently TAMkin only supports ideal gases in the translational partition
  function which means that :math:`K_c` does not depend on the pressure set in
  ``ExtTrans.pressure`` or the density set in ``ExtTrans.density``. (When no
  translational freedom is included in the partition function, there is no
  pressure or density to worry about.)

* **The change in free energy** is computed (and printed) as follows::

      delta_fr = tm.free_energy_change(temp)
      print "Change in free energy [kJ/mol] = %.3f" % (delta_fr/kjmol)

  The change in free energy does depend on the pressure (or density) in the
  translational part of the partition functions of the reactants and products.

* **The electronic energy difference** between the reactants (-) and the
  products (+) can be computed as follows::

      delta_e = tm.energy_difference(temp)

* **The zero-point energy difference** between the reactants (-) and the
  products (+) can also be computed::

      delta_zpe = tm.zero_point_energy_difference(temp)

A complete overview of the thermodynamic model, including the specification of
the partition functions, is written to file as follows::

    tm.write_to_file("thermodynamic_model.txt")


Reaction kinetics
~~~~~~~~~~~~~~~~~

See :class:`tamkin.chemmod.KineticModel` for the details.

All kinetic properties in TAMkin are computed with the ``KineticModel`` object.
For a given list of reactant partition functions, ``pfs_react`` and a transition
state partition function, ``pf_trans``, the ``KineticModel`` object is created
as follows::

      km = KineticModel(pfs_react, pf_trans)

The following kinetic properties can be computed with a ``KineticModel`` object:

* **The rate constant** is the main property of interest. The following two
  lines print the rate constant at 303K in SI units for a given kinetic model.
  ::

      rc = km.rate_constant(303)
      print "Rate constant [%s] at 303K = %.5e" % (km.unit_name, rc/km.unit)

  In the case of ideal gases, the rate constant does not depend on the pressure
  (or density) set in the translation partition functions.

* **The change in free energy** when going from reactants to products is
  computed as follows for a given temperature ``temp``::

      delta_fr = km.free_energy_change(temp)

* **The electronic energy difference** between the reactants (-) and the
  transition state (+) is computed as follows::

      delta_e = km.energy_difference(temp)

* **The zero-point energy difference** between the reactants (-) and the
  transition state (+) can also be computed::

      delta_zpe = km.zero_point_energy_difference(temp)

A complete overview of the kinetic model, including the specification of the
partition functions, is written to file as follows::

    km.write_to_file("kinetic_model.txt")


Tunneling corrections
---------------------

Three tunneling correction models are implemented in :mod:`tamkin.tunneling`.
One can include tunneling corrections in the kinetic model by giving a
``TunnelingCorrection`` object as the third argument to the ``KineticModel``
constructor. One must first create the tunneling object::

    # Just take one of the three following
    tunneling = Eckart(pfs_react, pf_trans, pfs_prod)
    tunneling = Wigner(pf_trans)
    tunneling = Miller(pf_trans)
    # Then create a kinetic model
    km = KineticModel(pfs_react, pf_trans, tunneling)

In this snippet ``pfs_prod`` is the list of product partition functions. All
rate constants computed with such a kinetic model include the tunneling
correction. The other properties computed by the kinetic model are not affected.


ReactionAnalysis objects -- fitting kinetic parameters A and E\ :sub:`a`
------------------------------------------------------------------------

See :class:`tamkin.pftools.ReactionAnalysis` for the details.

The ultimate purpose of a ``KineticModel`` object is to estimate the kinetic
parameters :math:`A` and :math:`E_a` in the empirical Arrhenius law. This can
be accomplished with a ReactionAnalysis object. As soon as a reaction analysis
object is created, the kinetic parameters are fitted and stored internally for
post-processing or direct output:

    ra = ReactionAnalysis(km, temp_low, temp_high)

The first argument is a kinetic model, which may include a tunneling correction.
The two following arguments are the boundaries for the grid on the temperature
axis used to fit the kinetic parameters. A fourth optional argument fixes the
spacing of the temperature grid. The default spacing is 10K.

The kinetic parameters and some related data are accessible as attributes:

* ``A`` and ``Ea`` -- The kinetic parameters in atomic units.
* ``R2`` -- The Pearson R^2 of the fit.
* ``temps`` -- An array with the temperature grid in Kelvin
* ``rate_consts`` -- the rate constants at the grid points in atomic units
* ``ln_rate_consts`` -- the logarithm of `the rate constants in atomic units`

The pre-exponential factor can be printed in SI units, knowing that it has the
same unit as the rate constant itself::

    print "Pre-exponential factor [%s] = %.5e" % (km.unit_name, ra.A/km.unit)
    print "The activation energy [kJ/mol] = %.2f" % (ra.Ea/kjmol)

When performing a high-throughput level-of-theory study, it is very convenient
to use the attributes of the ``ReactionAnalysis`` in a post-processing program.

One can write an overview of the reaction analysis to a text file as follows::

    ra.write_to_file("reaction_analysis.txt")

The computed rate constants and the linear fit can be plotted in a typical
Arrhenius plot as follows::

    ra.plot_arrhenius("arrhenius.txt")
