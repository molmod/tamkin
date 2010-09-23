Chemical Physics with TAMkin
============================

TODO: a chapter prior to this one, discussing the math of all supported
contributions to the partition function in TAMkin.


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

The steady state limit of a chemical reaction is completely characterized by
the equilibrium constant. It is one of the most important quantities that can
be derived from the partition functions in TAMkin.

In the case of ideal gases, this quantity only depends on the temperature, not
on the total pressure. For this reason, it is practically never necessary to set
the pressure in the translational contribution to the partition function.

Definition of the equilibrium constant
--------------------------------------

McQuarry
^^^^^^^^

It is instructive to review to the definition of the equilibrium constant given
in `Physical chemistry, a molecular approach`, by McQuarry and Simon
[McQuarry1997]_ (page 981). For a chemical reaction of the form

.. math:: \nu_A A(g) + \nu_b B(g) \rightleftharpoons \nu_C C(g) + \nu_D D(g)

the equilibrium constant in terms of concentrations is defined as

.. math:: K_c(T) = \frac{(Z_{1,C}/V)^{\nu_C}(Z_{1,D}/V)^{\nu_D}}
                        {(Z_{1,A}/V)^{\nu_A}(Z_{1,V}/V)^{\nu_B}},

where :math:`Z_{1,X}` is the single-particle partition function of species `X`
and V is the total volume of the system. One can derive the equilibrium constant
in terms of partial pressures using the ideal-gas law:

.. math:: K_p(T) = K_c(T) \left(\frac{c^0k_BT}{P_0}\right)^{\nu_C+\nu_D-\nu_A-\nu_B}.

Although this expressions for :math:`K_c` and :math:`K_p` are perfectly valid, they
are only applicable to the case where all reactants and products are 3D gas phase
particles sitting in the same reactor volume, :math:`V`. TAMkin also supports
partition functions for gases in other dimensions, or even for systems that have
no translational degrees of freedom at all. Moreover, for some applications, one
needs to find the equilibrium between systems that are physically discjunct
instead of sharing the same volume. Therefore we derive a more general
expression in the following section that coincides with the form of McQuarry in
the case of 3D gases.

General form
^^^^^^^^^^^^

Consider again the same chemical balance,

.. math:: \nu_A A + \nu_b B \rightleftharpoons \nu_C C + \nu_D D,

where we dropped the labels :math:`(g)`  as we do no longer consider the
only conventional gas phase systems. An extension with more reactions and
products is trivial.

The grand canonical partition function of this system is written as

.. math:: \mathcal{Z} = \sum_{N_A} \sum_{N_B} \sum_{N_C} \sum_{N_D}
                        Z(N_A, N_B, N_C, N_D, \ldots)

where :math:`Z` is the partition function for a fixed number of particles of
each species. We now introduce the first approximation, i.e. that the
interaction between the particles of different species can be neglected. This
means that the partition function :math:`Z` can be factorized into contributions
from partition functions per species:

.. math:: \mathcal{Z} = \sum_{N_A} \sum_{N_B} \sum_{N_C} \sum_{N_D}
                        Z_A(N_A, \ldots) Z_B(N_B, \ldots)
                        Z_C(N_C, \ldots) Z_D(N_D, \ldots)

where :math:`Z_X(N_X, \ldots)` is the parition function of a system with
:math:`N_X` reactants of species `X`. We do not need to know in detail what
kind of partition function :math:`Z_X` represents. It may be an NVT, NpT or any
other ensemble with a fixed number of particles.

The probability of a certain mixture of reactants is proportional to the product
of fixed particle partition functions:

.. math:: p(N_A, N_B, N_C, N_D) \propto Z_A(N_A, \ldots) Z_B(N_B, \ldots) Z_C(N_C, \ldots) Z_D(N_D, \ldots)

Now assume that we start from a reference state

.. math:: (N^0_A, N^0_B, N^0_C, N^0_D).

When we introduce a reaction coordinate :math:`\xi`, all other states
reachable through the chemical reaction can be written as

.. math:: (N^0_A - \xi\nu_A, N^0_B - \xi\nu_B, N^0_C + \xi\nu_C, N^0_D + \xi\nu_D)

To find the most probable system, the chemical equilibrium, we must find the
state that maximizes the probability :math:`p(N_A, N_B, N_C, N_D)`.
Mathematically, this means that we want to find a non-trivial solution to the
equation

.. math:: \frac{\partial p(N^0_A - \xi_{\text{eq}}\nu_A,
                           N^0_B - \xi_{\text{eq}}\nu_B,
                           N^0_C + \xi_{\text{eq}}\nu_C,
                           N^0_D + \xi_{\text{eq}}\nu_D)}
               {\partial \xi_{\text{eq}}} = 0.

To solve this problem, we rephrase it in terms of free energies, i.e. using
:math:`F_X = -k_Bt\ln(Z_X)` and the fact that the logarithmic function is
monotonous. The most probably state is therefore the state that minimizes the
total free energy.

.. math:: \frac{\partial (F_A(N^0_A - \xi_{\text{eq}}\nu_A, \ldots)
                         +F_B(N^0_B - \xi_{\text{eq}}\nu_B, \ldots)
                         +F_C(N^0_C + \xi_{\text{eq}}\nu_C, \ldots)
                         +F_D(N^0_D + \xi_{\text{eq}}\nu_D, \ldots)}
               {\partial \xi_{\text{eq}}} = 0

Using the the definition of the chemical potential, :math:`\mu(N_X, \ldots) =
\frac{\partial F_X(N_X, \ldots)}{\partial N_X}`, we end up with a very familiar
expression for the equilibrium condition:

.. math:: \nu_C \mu_C(N_{C,\text{eq}}, \ldots) + \nu_D \mu_D(N_{D,\text{eq}}, \ldots)
          - \nu_A \mu_A(N_{A,\text{eq}}, \ldots) - \nu_B \mu_B(N_{B,\text{eq}}, \ldots) = 0

where :math:`N_{X, \text{eq}}` is a shorthand for :math:`N^0_{X} +
\xi_{\text{eq}}\nu_X`. Now we rephrase these equations back in terms of the
partition functions. We rely on the classical gas limit of many-particle
partition function:

.. math::
    :nowrap:

    \begin{align*}
      \mu_X & = -k_BT \left(\frac{\partial \ln(Z_X(N_X, \ldots)}{\partial N_X}\right) \\
            & = -k_BT \left(\frac{\partial \ln\left(\frac{Z^{N_X}_X(1, \ldots)}{N_X!}\right)}{\partial N_X}\right) \\
            & = -k_BT \left(\frac{\partial (N_X\ln(Z_X(1, \ldots) - N_X\ln(N_X) + N_X)}{\partial N_X}\right) \\
            & = -k_BT \ln\left(\frac{Z_X(1, \ldots)}{N_X}\right)
    \end{align*}

**TODO:** This only valid when :math:`Z_X(1, \ldots)` does not explicitly depend
on the :math:`N_X`, which is only the case in constant volume ensembles. This
derivation should be generalized.

This expression for the chemical potential can be plugged back into the
equilibrium condition to get

.. math:: \frac{N_{C,\text{eq}}^{\nu_C}\,N_{D,\text{eq}}^{\nu_D}}
               {N_{A,\text{eq}}^{\nu_A}\,N_{B,\text{eq}}^{\nu_B}} =
          \frac{Z_C(1, \ldots)^{\nu_C}\,Z_D(1, \ldots)^{\nu_D}}
               {Z_A(1, \ldots)^{\nu_A}\,Z_B(1, \ldots)^{\nu_B}},

which is a standard text-book result. Now comes the hard part, where we have to
keep the derivation general enough to cover 3D gases, 2D gases, and systems
without translational freedom. In each case we must introduce a definition of
a density, which is required for a general expression of :math:`K_c`:

* **3D gas**: :math:`\rho_X = N_X/V_X`, where :math:`V_X` is the volume of the
  system containing particles of species X.

* **2D gas**: :math:`\rho_X = N_X/A_X`, where :math:`A_X` is the area of the
  system containing particles of species X.

* **Non-translational**: :math:`\rho_X = N_X`, which is simply the occupation
  number of the site X, or the probability that it is occupied. In the classical
  limit, this number is always well below unity.

In analogy, we must introduce different types of `dimensionless` partition
functions:

* **3D gas**: :math:`Z'_X(1, \ldots) = Z_X(1, \ldots)/V_X`, where :math:`V_X` is
  the volume of the system containing particles of species X.

* **2D gas**: :math:`Z'_X(1, \ldots) = Z_X(1, \ldots)/A_X`, where :math:`A_X` is
  the area of the system containing particles of species X.

* **Non-translational**: :math:`Z'_X(1, \ldots) = Z_X(1, \ldots)`.

We can finally write down the general form of :math:`K_c`:

.. math:: K_c(T) = \frac{\rho_{C,\text{eq}}^{\nu_C}\,\rho_{D,\text{eq}}^{\nu_D}}
                        {\rho_{A,\text{eq}}^{\nu_A}\,\rho_{B,\text{eq}}^{\nu_B}}
                 = \frac{Z'^{\nu_C}_C(1, \ldots)\,Z'^{\nu_D}_D(1, \ldots)}
                        {Z'^{\nu_A}_A(1, \ldots)\,Z'^{\nu_B}_B(1, \ldots)}


Implementation in TAMkin
^^^^^^^^^^^^^^^^^^^^^^^^

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


The unit of :math:`K_c`
^^^^^^^^^^^^^^^^^^^^^^^

By construction :math:`K_c` is no longer a dimensionless quantity. It's unit is
defined by the partition functions that go into the equilibrium constant.

- For each gas phase reactant, there is a factor bohr\ :sup:`d`, where `d` is
  the dimension of the gas.
- For a each gas phase product, there is a factor bohr\ :sup:`-d`, where `d` is
  the dimension of the gas.

In SI units, this becomes:

- For each gas phase reactant, there is a factor meter\ :sup:`d` mol\ :sup:`-1`,
  where `d` is the dimension of the gas.
- For a each gas phase product, there is a factor meter\ :sup:`-d` mol, where
  `d` is the dimension of the gas.

The standard change in free energy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

TODO


Computation of the equilibrium constant
---------------------------------------

Given a list of partition functions of reactants (``pfs_react``) and a list of
product partition functions (``pfs_prod``), the equilibrium constant is computed
at a certain temperature, ``temp``, as follows::

    K = compute_equilibrium_constant(pfs_react, pfs_prod, temp)

This function takes one optional argument: ``do_log``, which is by default
``False``. When set to True, the logarithm of the partition function is
returned.

Computation


Computation of the standard change in free energy
-------------------------------------------------

TODO


ThermodynamicModel objects
--------------------------

For later reference it is convenient to mention the ``ThermodynamicModel``
class. It is a simple object oriented representation of a thermodynamic
equilibrium.

Given a list of partition functions of reactants (``pfs_react``) and a list of
product partition functions (``pfs_prod``), a ``ThermodynamicModel`` object is
created as follows::

    tm = ThermodynamicModel(pfs_react, pfs_prod)

This can be used to compute the equilibrium constant as function of the
temperature::

    print "Equilibrium constant at 300K.", tm.compute_equilibrium_constant(300)

**TODO:** Add the standard change in free energy.

Reaction kinetics
~~~~~~~~~~~~~~~~~

Definition of the equilibrium constant
--------------------------------------

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
