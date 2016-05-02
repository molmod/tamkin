.. _releases:

Release history
###############


* **May 2, 2016. Version 1.0.9**

    - Improved VASP IO routine, option to pass OUTCAR file of single point calculation
      used for a refined energy calculation.
    - Use two digits in temperature filename of ``tamkin-driver.py``.

* **February 5, 2016. Version 1.0.8**

    - The script tamkin-driver.py now reads the fixed atoms from the Gaussian fchk file
      and carries out a PHVA analysis if fixed atoms are present.
    - The ``load_indices`` now also handles ranges of atoms.

* **October 30, 2015. Version 1.0.7**

    - Improved VASP IO (partial Hessians and compatibility with different VASP versions)
    - Improved CP2K IO (partial Hessians and array with fixed atoms)
    - More output is written by the script tamkin-driver.py
    - Fixed mistakes in the installation instructions.

* **September 27, 2014. Version 1.0.6**

    - More robust internal rotor code

* **September 9, 2014. Version 1.0.5**

    - Several bug fixes and cleanups

* **August 8, 2014. Version 1.0.4**

    - Improved ``tamkin-driver.py`` with improved documentation.
    - Added missing test file.

* **August 7, 2014. Version 1.0.3**

    - Fix bug in load_rotscan_g03log. This function now also works with G09.
    - More output (rotors and tunneling correction)
    - Initial version of the TAMkin driver script

* **August 5, 2014. Version 1.0.2**

    - Add some statistics to several plots
    - Make the output of dump_modes_gaussian compatible with Avogadro.

* **April 2, 2014. Version 1.0.1**

    - Bug fix related to reading CP2K output files.
    - Bug fix for linear molecules when option ``treatment=Full()`` is used in
      the normal mode analysis.

* **December 6, 2013. Version 1.0.0**

    - Initial release through github.
