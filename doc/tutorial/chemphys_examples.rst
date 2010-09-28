Chemical Physics -- Practical examples
======================================

**TODO**: add some more introductory text.


Conformational Equilibrium
~~~~~~~~~~~~~~~~~~~~~~~~~~

**TODO**


Chemical Equilibrium
~~~~~~~~~~~~~~~~~~~~

**TODO**


Heat of formation
~~~~~~~~~~~~~~~~~

**TODO**


Reaction Kinetics (bimolecular)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this example we show how one estimates kinetic parameters for the addition of
ethene to ethyl in the gas phase at constant pressure. The reaction balance is

.. math::
    :nowrap:

    \newcommand{\rad}{\raisebox{0.4ex}{\scriptsize{\ensuremath{\bullet}}}}

    CH$_2$=CH$_2$ (gas) + \rad{C}H$_2-$CH$_3$ (gas) $\rightarrow$ \rad{C}H$_2$-CH$_2$-CH$_2$-CH$_3$ (gas)

For this example we prepared three frequency computations:

- One for each ground state geometry of the reactants (ethene and ethyl). The
  formatted checkpoint files of the frequency jobs are ``ethene.fchk`` and
  ``ethyl.fchk``.

- One for the transition state where ethene performs a `trans` attack on ethyl.
  The geometry of the transition state is optimized towards the saddle point
  in the potential energy surface. The formatted checkpoint file of the
  frequency job is ``ts_trans.fchk``.

The frequency computations are carried out with Gaussian03. The level of theory
is B3LYP/6-31G(d). The following script computes the kinetic parameters (A and
E\ :sub:`a`) through a linear fit of :math:`\ln(k)` versus :math:`T` in the
temperature range 300K-600K.

File: ``examples/019_ethyl_ethene_simple/kinetic.py``

.. literalinclude:: ../../examples/019_ethyl_ethene_simple/kinetic.py
   :lines: 37-
   :linenos:

The scripts writes several output files discussed in the subsections below.

CSV Files with energetic analysis
---------------------------------

CSV files are created for different temperatures: 300K, 400K, 500K and 600K.
The file at 300 K contains the following data:

.. csv-table::

    Temperature [K],300,,,
    ,,,,
    **Quantity**,**Ethyl**,**Ethene**,**Transition state**,**Linear combination (always in kJ/mol)**
    Signed stoichiometry,-1,-1,1,
    **Values in a.u.**,,,,
    Electronic energy,-79.1579,-78.5875,-157.7371,22
    Zero-point energy,-79.0982,-78.5362,-157.6231,30
    Internal heat (300.00K),-79.0933,-78.5322,-157.6157,26
    Chemical potential (300.00K),-79.1225,-78.5573,-157.6536,69
    **Corrections in kJ/mol**,,,,
    Zero-point energy,157,134,299,8
    Internal heat (300.00K),170,145,319,4
    Chemical potential (300.00K),93,79,219,47
    ,,,,
    **Other quantities**,Unit,Value,,
    Rate constant,m**3*mol**-1/second,4.52253403913e-11,,

The numbers in this table are rounded at some precision to improve the
readability, but the actual CSV file contains all numbers in full machine
precision. The linear combination of the chemical potentials is also known as
the `change in free energy` associated with the reaction.

An Arrhenius plot
-----------------

This plot can be used for a visual check of the linear regression analysis to
estimate the kinetic parameters

.. image:: arrhenius_bimolecular_gas_phase.png

A log file with an analysis of the kinetic parameters
-----------------------------------------------------

This file is written to the file ``kinetic_parameters.txt``. It contains the
following data::

    Summary
    A [m**3*mol**-1/second] = 8.93643e+04
    ln(A [a.u.]) = -10.63
    Ea [kJ/mol] = 33.16
    R2 (Pearson) = 99.94%

    Temperature grid
    T_low [K] = 300.0
    T_high [K] = 600.0
    T_step [K] = 10.0
    Number of temperatures = 31

    Reaction rate constants
        T [K]    Delta_r F [kJ/mol]      k(T) [m**3*mol**-1/second]
        300.00          68.7              1.66848e-01
        310.00          70.1              2.48896e-01
        320.00          71.6              3.62870e-01
        330.00          73.0              5.18109e-01
        340.00          74.4              7.25806e-01
        350.00          75.9              9.99197e-01
        360.00          77.3              1.35374e+00
        370.00          78.7              1.80731e+00
        380.00          80.1              2.38032e+00
        390.00          81.6              3.09593e+00
        400.00          83.0              3.98015e+00
        410.00          84.4              5.06200e+00
        420.00          85.8              6.37365e+00
        430.00          87.3              7.95046e+00
        440.00          88.7              9.83115e+00
        450.00          90.1              1.20579e+01
        460.00          91.5              1.46762e+01
        470.00          92.9              1.77354e+01
        480.00          94.4              2.12883e+01
        490.00          95.8              2.53912e+01
        500.00          97.2              3.01043e+01
        510.00          98.6              3.54916e+01
        520.00         100.0              4.16205e+01
        530.00         101.4              4.85623e+01
        540.00         102.8              5.63922e+01
        550.00         104.2              6.51888e+01
        560.00         105.6              7.50346e+01
        570.00         107.1              8.60158e+01
        580.00         108.5              9.82222e+01
        590.00         109.9              1.11747e+02
        600.00         111.3              1.26688e+02

    Electronic energy barrier [kJ/mol] = 21.6
    Zero-point energy barrier [kJ/mol] = 29.7

    Reactant 0 partition function
    Title: Ethyl
    Electronic energy [au]: -79.15787
    Zero-point contribution [kJ/mol]: 156.6101213
    Zero-point energy [au]: -79.09822
    Contributions to the partition function:
      ELECTRONIC
        Multiplicity: 2
        Electronic energy: -79.1578683
      ROTATIONAL
        Rotational symmetry number: 1
        Moments of inertia [amu*bohr**2]: 17.468474  79.684868 85.942337
        Threshold for non-zero moments of inertia [amu*bohr**2]: 5.485799e-04
        Non-zero moments of inertia: 3
      TRANSLATIONAL
        Dimension: 3
        Constant pressure: True
        Pressure [bar]: 1.01325
          BIG FAT WARNING!!!
          This is an NpT partition function.
          Internal energy contains a PV term (and is therefore the enthalpy).
          Free energy contains a PV term (and is therefore the Gibbs free energy).
          The heat capacity is computed at constant pressure.
        Mass [amu]: 29.039125
      VIBRATIONAL
        Number of zero wavenumbers: 0
        Number of real wavenumbers: 15
        Number of imaginary wavenumbers: 0
        Frequency scaling factor: 1.0000
        Zero-point scaling factor: 1.0000
        Real Wavenumbers [1/cm]:
           123.7   457.5   817.9   995.0  1074.2  1207.6  1430.1  1492.4
          1510.8  1514.8  2965.5  3058.3  3102.3  3168.1  3264.7
        Zero-point contribution [kJ/mol]: 156.6101213

    Reactant 1 partition function
    Title: Ethene
    Electronic energy [au]: -78.58746
    Zero-point contribution [kJ/mol]: 134.4868825
    Zero-point energy [au]: -78.53624
    Contributions to the partition function:
      ELECTRONIC
        Multiplicity: 1
        Electronic energy: -78.5874587
      ROTATIONAL
        Rotational symmetry number: 4
        Moments of inertia [amu*bohr**2]: 12.280076  60.075552 72.355628
        Threshold for non-zero moments of inertia [amu*bohr**2]: 5.485799e-04
        Non-zero moments of inertia: 3
      TRANSLATIONAL
        Dimension: 3
        Constant pressure: True
        Pressure [bar]: 1.01325
          BIG FAT WARNING!!!
          This is an NpT partition function.
          Internal energy contains a PV term (and is therefore the enthalpy).
          Free energy contains a PV term (and is therefore the Gibbs free energy).
          The heat capacity is computed at constant pressure.
        Mass [amu]: 28.031300
      VIBRATIONAL
        Number of zero wavenumbers: 0
        Number of real wavenumbers: 12
        Number of imaginary wavenumbers: 0
        Frequency scaling factor: 1.0000
        Zero-point scaling factor: 1.0000
        Real Wavenumbers [1/cm]:
           834.8   956.1   976.1  1070.1  1248.0  1395.8  1494.3  1720.2
          3151.9  3167.3  3222.2  3247.7
        Zero-point contribution [kJ/mol]: 134.4868825

    Transition state partition function
    Title: Transition state
    Electronic energy [au]: -157.73711
    Zero-point contribution [kJ/mol]: 299.2533370
    Zero-point energy [au]: -157.62313
    Contributions to the partition function:
      ELECTRONIC
        Multiplicity: 2
        Electronic energy: -157.7371095
      ROTATIONAL
        Rotational symmetry number: 1
        Moments of inertia [amu*bohr**2]: 92.846631  597.569081 642.613097
        Threshold for non-zero moments of inertia [amu*bohr**2]: 5.485799e-04
        Non-zero moments of inertia: 3
      TRANSLATIONAL
        Dimension: 3
        Constant pressure: True
        Pressure [bar]: 1.01325
          BIG FAT WARNING!!!
          This is an NpT partition function.
          Internal energy contains a PV term (and is therefore the enthalpy).
          Free energy contains a PV term (and is therefore the Gibbs free energy).
          The heat capacity is computed at constant pressure.
        Mass [amu]: 57.070425
      VIBRATIONAL
        Number of zero wavenumbers: 0
        Number of real wavenumbers: 32
        Number of imaginary wavenumbers: 1
        Frequency scaling factor: 1.0000
        Zero-point scaling factor: 1.0000
        Real Wavenumbers [1/cm]:
            48.5   154.6   157.1   247.0   370.8   547.0   765.2   823.5
           831.9   848.8   917.8  1024.8  1035.8  1075.2  1228.4  1247.8
          1317.9  1432.1  1487.2  1498.5  1514.5  1518.4  1609.4  2985.9
          3061.7  3100.7  3149.1  3153.7  3163.8  3225.6  3237.2  3251.4
        Imaginary Wavenumbers [1/cm]:
          -383.6
        Zero-point contribution [kJ/mol]: 299.2533370
