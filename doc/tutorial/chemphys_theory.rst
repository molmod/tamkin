Chemical Physics -- Theoretical Background
==========================================


The partition function of gas-phase molecules
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Introduction
------------

The canonical partition function of a molecular gas in the classical and ideal
gas limit can be written in terms of the partition function of a single gas
molecule:

.. math:: Z_N = \frac{Z_1^N}{N!}.

The subscript :math:`N` and :math:`1` distinguish between the many- and
single-particle partition function. We will also use the notation :math:`Z(N,
\ldots)` and :math:`Z(1, \ldots)`, where the dots stand for the
other boundary conditions (or natural variables) that define the partition
function. This can be volume (:math:`V`) and temperature (:math:`T`), or
pressure (:math:`p`) and temperature (:math:`T`), or yet something else.

The second approximation is the separation of the single-particle partition
function into factors:

.. math:: Z_1 = Z_{1,\text{elec}} Z_{1,\text{trans}} Z_{1,\text{rot}} Z_{1,\text{vib}} \ldots

Only the four most conventional factors are shown here, but an extension with
internal rotors and other contributions is also possible. These factors stand
for:

* :math:`Z_{1,\text{elec}}`: the electronic contribution, i.e. electronic ground
  state energy of the molecule and the multiplicity.
* :math:`Z_{1,\text{trans}}`: the translational contribution, taking into
  account constant volume or constant pressure boundary conditions.
* :math:`Z_{1,\text{rot}}`: the rotational contribution.
* :math:`Z_{1,\text{vib}}`: the vibrational contribution.

All supported contributions in TAMkin will be discussed in detail below. This
separation into factors is an approximation because the corresponding terms in
the molecular Hamiltonian are not completely decoupled. In practice one
constructs a partition function for each subsystem under the assumption that all
other degrees of freedom are in their ground state.


Quantities of interest
----------------------

In the following subsections, the contribution to the following quantities of
each factor in the partition function will be discussed:

.. math::
    :nowrap:

    \begin{align*}
        \mathsf{log} & = \frac{\ln Z_N }{N} \\
        \mathsf{logt} & = \frac{1}{N}\frac{\partial \ln Z_N}{\partial T} \\
        \mathsf{logtt} & = \frac{1}{N}\frac{\partial^2 \ln Z_N}{\partial T^2} \\
        \mathsf{logn} & = \frac{\partial \ln Z_N}{\partial N} \\
        \mathsf{logv} & = \frac{\partial \ln Z_N}{\partial N} - \ln\left(\frac{V}{N}\right)
    \end{align*}


These are all intensive quantities. All other thermodynamic functions will be
derived from these basic quantities in a general fashion.

The first step is to rewrite these quantities in terms of the
single-particle partition function, using the classical gas limit and
Stirling's approximation:

.. math::
    :nowrap:

    \begin{align*}
        \mathsf{log} & = \frac{N\ln Z_1 - N\ln N + N}{N} = 1 + \ln\left(\frac{Z_1}{N}\right) \\
        \mathsf{logt} & = \frac{\partial \ln Z_1 }{\partial T} \\
        \mathsf{logtt} & = \frac{\partial^2 \ln Z_1 }{\partial T^2} \\
        \mathsf{logn} & = \frac{\partial [N \ln Z_1 - N\ln N + N]}{\partial N} = \ln\left(\frac{Z_1}{N}\right) + N\frac{\partial \ln Z_1}{\partial N} \\
        \mathsf{logv} & = \ln\left(\frac{Z_1}{N}\right) + N\frac{\partial \ln Z_1}{\partial N} - \ln\left(\frac{V}{N}\right)
    \end{align*}

In these equations, one can fill in the factorization of the single-particle
partition function. We use the following strategy to make sure that the individual
contributions are all intensive quantities.

1. The translational contribution will include all the `effects` of the
   many-body nature of the total partition function.

    .. math::
        :nowrap:

        \begin{align*}
            \mathsf{log}_{\text{trans}} & = 1 + \ln\left(\frac{Z_{1,\text{trans}}}{N}\right) \\
            \mathsf{logt}_{\text{trans}} & = \frac{\partial \ln Z_{1,\text{trans}} }{\partial T} \\
            \mathsf{logtt}_{\text{trans}} & = \frac{\partial^2 \ln Z_{1,\text{trans}} }{\partial T^2} \\
            \mathsf{logn}_{\text{trans}} & = \ln\left(\frac{Z_{1,\text{trans}}}{N}\right) + N\frac{\partial \ln Z_{1,\text{trans}}}{\partial N} \\
            \mathsf{logv}_{\text{trans}} & = \ln\left(\frac{Z_{1,\text{trans}}}{N}\right) + N\frac{\partial \ln Z_{1,\text{trans}}}{\partial N} - \ln\left(\frac{V}{N}\right)
        \end{align*}

2. All other contributions are treated as if they are simple single-particle
   contributions.

    .. math::
        :nowrap:

        \begin{align*}
            \mathsf{log}_{\text{other}} & = \ln Z_{1,\text{other}} \\
            \mathsf{logt}_{\text{other}} & = \frac{\partial \ln Z_{1,\text{other}} }{\partial T} \\
            \mathsf{logtt}_{\text{other}} & = \frac{\partial^2 \ln Z_{1,\text{other}} }{\partial T^2} \\
            \mathsf{logn}_{\text{other}} & = \ln Z_{1,\text{other}} + N\frac{\partial \ln Z_{1,\text{other}}}{\partial N} \\
            \mathsf{logv}_{\text{other}} & = \ln Z_{1,\text{other}} + N\frac{\partial \ln Z_{1,\text{other}}}{\partial N}
        \end{align*}

This strategy has the additional advantage that particles without translational
degrees of freedom can be treated within the same framework. For such systems,
the classical gas limit does not apply either and one has :math:`Z_N = Z_1^N`.
One can simply drop the translational contributions to :math:`\mathsf{log*}`.

Electronic contribution
-----------------------


Translational contribution
--------------------------


Rotational contribution
-----------------------


Vibrational contribution
------------------------


Rigid free rotor correction
---------------------------


Hindered free rotor correction
------------------------------


Quantities derived from one partition function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Free energy
-----------


Internal energy
---------------


Heat capacity
-------------


Entropy
-------


Chemical potential
------------------


Quantities derived from multiple partition functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


The equilibrium constant
------------------------

The steady state limit of a chemical reaction is completely characterized by
the equilibrium constant. It is one of the most important quantities that can
be derived from the partition functions in TAMkin.

In the case of ideal gases, this quantity only depends on the temperature, not
on the total pressure. For this reason, it is practically never necessary to set
the pressure in the translational contribution to the partition function.


McQuarrie
^^^^^^^^^

It is instructive to review to the definition of the equilibrium constant given
in `Physical chemistry, a molecular approach`, by McQuarrie and Simon
[McQuarrie1997]_ (page 981). For a chemical reaction of the form

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
needs to find the equilibrium between systems that are physically disjunct
instead of sharing the same volume. Therefore we derive a more general
expression in the following section that coincides with the form of McQuarrie in
the case of 3D gases.

General form
^^^^^^^^^^^^

Consider again the same chemical balance,

.. math:: \nu_A A + \nu_b B \rightleftharpoons \nu_C C + \nu_D D,

where we dropped the labels :math:`(g)` as we do no longer consider the
only conventional gas phase systems. An extension with more reactions and
products is trivial. We assume that this reaction takes place in a closed
system, e.g. a reactor vessel, where the number of particles of each species may
only change through the chemical reaction. All possible states in the closed
system are known once we assume a reference state

.. math:: (N^0_A, N^0_B, N^0_C, N^0_D, \ldots).

where :math:`N^0_X` is the reference number of particles of species `X`. The
dots stand for all other natural variables of the closed system, .e.g. total
volume or external pressure, which remain constant during the course of the
reaction.

When we introduce a reaction coordinate :math:`\xi`, all other states
reachable through the chemical reaction can be written as

.. math:: (N^0_A - \xi\nu_A, N^0_B - \xi\nu_B, N^0_C + \xi\nu_C, N^0_D + \xi\nu_D, \ldots)

The grand partition function for all states of the mixture is written as:

.. math:: \mathcal{Z} = \sum_{\xi = \xi_{\text{min}}}^{\xi_{\text{max}}}
                Z(N^0_A - \xi\nu_A, N^0_B - \xi\nu_B, N^0_C + \xi\nu_C, N^0_D + \xi\nu_D, \ldots)

where :math:`Z` is the partition function of the mixture at a fixed reaction
coordinate. Assuming that the interactions between particles of different
species can be neglected, the grand partition function becomes:

.. math:: \mathcal{Z} = \sum_{\xi = \xi_{\text{min}}}^{\xi_{\text{max}}}
                Z_A(N^0_A - \xi\nu_A, \ldots)
                Z_B(N^0_B - \xi\nu_B, \ldots)
                Z_C(N^0_C + \xi\nu_C, \ldots)
                Z_D(N^0_D + \xi\nu_D, \ldots)

where :math:`Z_X(N_X, \ldots)` is the partition function of a system with
:math:`N_X` reactants of species `X`. We do not need to know in detail
what kind of partition function :math:`Z_X` represents. It may be an NVT, NpT or
any other ensemble with a fixed number of particles.

The probability of a mixture of reactants and products is proportional to the
product of fixed particle partition functions:

.. math:: p(N_A, N_B, N_C, N_D) \propto Z_A(N_A, \ldots) Z_B(N_B, \ldots) Z_C(N_C, \ldots) Z_D(N_D, \ldots)

where :math:`N_X` is a shorthand for :math:`N^0_{X} + \xi\nu_X`. To find the
most probable state of the system, the chemical equilibrium, we must find the
:math:`\xi` that maximizes the probability :math:`p(N_A, N_B, N_C, N_D)`.
Mathematically, this means that we want to find a non-trivial solution to the
equation

.. math:: \frac{\partial p(N_A, N_B, N_C, N_D + \xi_{\text{eq}}\nu_D)}
               {\partial \xi_{\text{eq}}} = 0.

To solve this problem, we rephrase it in terms of free energies, i.e. using
:math:`F_X = -k_Bt\ln(Z_X)` and the fact that the logarithmic function is
monotonous. The most probable state is therefore the state that minimizes the
total free energy.

.. math:: \frac{\partial [F_A(N^0_A - \xi_{\text{eq}}\nu_A, \ldots)
                         +F_B(N^0_B - \xi_{\text{eq}}\nu_B, \ldots)
                         +F_C(N^0_C + \xi_{\text{eq}}\nu_C, \ldots)
                         +F_D(N^0_D + \xi_{\text{eq}}\nu_D, \ldots)]}
               {\partial \xi_{\text{eq}}} = 0

Using the the definition of the chemical potential, :math:`\mu_X(N_X, \ldots) =
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
            & = -k_BT \left(\frac{\partial [N_X\ln(Z_X(1, \ldots)) - N_X\ln(N_X) + N_X]}{\partial N_X}\right) \\
            & = -k_BT \ln\left(\frac{Z_X(1, \ldots)}{N_X}\right)
    \end{align*}

The last step is only valid when :math:`Z_X(1, \ldots)` does not
explicitly depend on the :math:`N_X`, which is only true for ideal gases.

This expression for the chemical potential can be plugged back into the
equilibrium condition to get

.. math:: \frac{N_{C,\text{eq}}^{\nu_C}\,N_{D,\text{eq}}^{\nu_D}}
               {N_{A,\text{eq}}^{\nu_A}\,N_{B,\text{eq}}^{\nu_B}} =
          \frac{Z^{\nu_C}_C(1, \ldots)\,Z^{\nu_D}_D(1, \ldots)}
               {Z^{\nu_A}_A(1, \ldots)\,Z^{\nu_B}_B(1, \ldots)},

which is a standard text-book equation, but now derived in a much more general
context. Now comes the hard part, where we have to keep the derivation general
enough to cover 3D gases, 2D gases, and systems without translational freedom.
In each case we must introduce a definition of a density, which is required for
a general expression of :math:`K_c`:

* **3D gas**: :math:`\rho_X = N_X/V_X`, where :math:`V_X` is the volume of the
  system containing particles of species X.

* **2D gas**: :math:`\rho_X = N_X/A_X`, where :math:`A_X` is the area of the
  system containing particles of species X.

* **Non-translational**: :math:`\rho_X = N_X`, which is simply the occupation
  number of the site X, or the probability that it is occupied. In the classical
  limit, this number is always well below unity.

In analogy, we must introduce different types of `dimensionless auxiliary
partition functions`:

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

The unit of K\ :sub:`c`
^^^^^^^^^^^^^^^^^^^^^^^

By construction :math:`K_c` is no longer a dimensionless quantity. (This is
different from the approach followed by McQuarrie, where :math:`K_c` is made
dimensionless by assuming some reference concentration for each quantity.)
The unit of :math:`K_c` is defined by the partition functions that go into the
equilibrium constant.

- For each gas phase reactant, there is a factor :math:`\text{bohr}^d`, where `d` is
  the dimension of the gas.
- For a each gas phase product, there is a factor :math:`\text{bohr}^{-d}`, where `d` is
  the dimension of the gas.

In SI units, this becomes:

- For each gas phase reactant, there is a factor :math:`m^d\,mol^{-1}`,
  where `d` is the dimension of the gas.
- For a each gas phase product, there is a factor :math:`mol\,m^{-d}`, where
  `d` is the dimension of the gas.


The change in free energy
-------------------------

The change in free energy associated with a reaction, :math:`\Delta_r G`, is
defined as the chemical potential of the products minus the chemical potential
of the reactants

.. math:: \Delta_r G(T) = \nu_C \mu_C(N_C, \ldots) + \nu_D \mu_D(N_D, \ldots)
                        - \nu_A \mu_A(N_A, \ldots) - \nu_B \mu_B(N_B, \ldots)

where the chemical potentials are all computed at a certain well-defined state
of the ensemble. For example, for 3D gases, :math:`\Delta_r G` depends on the
pressure and the temperature. When the number is expressed in Hartree/particle,
it is the free energy required to transform :math:`\nu_A` molecules of reactant
A and :math:`\nu_B` molecules of reactant B into :math:`\nu_C` molecules of
product C and :math:`\nu_D` molecules of product B, at a certain reference
state.

Let us now use the relation

.. math:: \mu_X = -k_BT\ln\left(\frac{Z_X(1,\ldots)}{N_X}\right)

to rewrite the change in free energy in terms of partition functions.

.. math:: \Delta_r G(T) = -k_BT \ln\left(
                \frac{Z^{\nu_C}_C(1,\ldots) Z^{\nu_D}_D(1,\ldots)}
                     {Z^{\nu_A}_A(1,\ldots) Z^{\nu_B}_B(1,\ldots)}
                \frac{N^{\nu_A}_A N^{\nu_B}_B}{N^{\nu_C}_C N^{\nu_D}_D}
            \right)

We now assume a reference state for each partition function that leads to a
reference `density`, :math:`\rho_{X,0}`, for each subsystem. The meaning the term
`density` may depend on the dimension of the gas, as discussed previously. We
can further rewrite the change in free energy as:

.. math:: \Delta_r G(T) = -k_BT \ln\left(
                \frac{Z'^{\nu_C}_C(1,\ldots) Z'^{\nu_D}_D(1,\ldots)}
                     {Z'^{\nu_A}_A(1,\ldots) Z'^{\nu_B}_B(1,\ldots)}
                \frac{\rho^{\nu_A}_{A,0} \rho^{\nu_B}_{B,0}}{\rho^{\nu_C}_{C,0} \rho^{\nu_D}_{D,0}}
            \right).

The first factor in the logarithm is the equilibrium constant, so we get:

.. math:: \Delta_r G(T) = -k_BT \ln\left(
                K_c \frac{\rho^{\nu_A}_{A,0} \rho^{\nu_B}_{B,0}}
                         {\rho^{\nu_C}_{C,0} \rho^{\nu_D}_{D,0}}
            \right).


:math:`K_c` is (for ideal gases) independent of the density or pressure of each
component. It still depends on the temperature. The second factor does not depend on
temperature, and bundles all the density or pressure information of the
reference state at which the change in free energy is computed.


The rate constant
-----------------


Kinetic parameters (A and E\ :sub:`a`)
--------------------------------------
