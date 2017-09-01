..
    : TAMkin is a post-processing toolkit for normal mode analysis, thermochemistry
    : and reaction kinetics.
    : Copyright (C) 2008-2012 Toon Verstraelen <Toon.Verstraelen@UGent.be>, An Ghysels
    : <An.Ghysels@UGent.be> and Matthias Vandichel <Matthias.Vandichel@UGent.be>
    : Center for Molecular Modeling (CMM), Ghent University, Ghent, Belgium; all
    : rights reserved unless otherwise stated.
    :
    : This file is part of TAMkin.
    :
    : TAMkin is free software; you can redistribute it and/or
    : modify it under the terms of the GNU General Public License
    : as published by the Free Software Foundation; either version 3
    : of the License, or (at your option) any later version.
    :
    : In addition to the regulations of the GNU General Public License,
    : publications and communications based in parts on this program or on
    : parts of this program are required to cite the following article:
    :
    : "TAMkin: A Versatile Package for Vibrational Analysis and Chemical Kinetics",
    : An Ghysels, Toon Verstraelen, Karen Hemelsoet, Michel Waroquier and Veronique
    : Van Speybroeck, Journal of Chemical Information and Modeling, 2010, 50,
    : 1736-1750W
    : http://dx.doi.org/10.1021/ci100099g
    :
    : TAMkin is distributed in the hope that it will be useful,
    : but WITHOUT ANY WARRANTY; without even the implied warranty of
    : MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    : GNU General Public License for more details.
    :
    : You should have received a copy of the GNU General Public License
    : along with this program; if not, see <http://www.gnu.org/licenses/>
    :
    : --

.. _releases:

Release history
###############

* **Version 1.2.3** September 1, 2017

    - First release in molmod channel on Anaconda.

* **Version 1.2.2** August 10, 2017

    - Fix typos in the docs.

* **Version 1.2.1** August 09, 2017

    - Git is no longer needed when installing from pip.

* **Version 1.2.0** August 04, 2017

    - Restructured package, such that it can be installed with setuptools
    - Tests are installed as well
    - Packaging for PyPI, Github source, anaconda.org using Travis Deployment
    - Deployment of website with Travis
    - Lowercase package filenames

* **Version 1.1.0** July 9, 2017

    - Fix bug in CP2K fixed atom loader.
    - Fixed link to CP2K website.
    - Improved format for atom indices: comments are allowed and shift can be set in the
      file.
    - tamkin-driver support for CP2K and VASP
    - tesing on Travis-CI

* **Version 1.0.9** May 2, 2016

    - Improved VASP IO routine, option to pass OUTCAR file of single point calculation
      used for a refined energy calculation.
    - Use two digits in temperature filename of ``tamkin-driver.py``.

* **Version 1.0.8** February 5, 2016

    - The script tamkin-driver.py now reads the fixed atoms from the Gaussian fchk file
      and carries out a PHVA analysis if fixed atoms are present.
    - The ``load_indices`` now also handles ranges of atoms.

* **Version 1.0.7** October 30, 2015

    - Improved VASP IO (partial Hessians and compatibility with different VASP versions)
    - Improved CP2K IO (partial Hessians and array with fixed atoms)
    - More output is written by the script tamkin-driver.py
    - Fixed mistakes in the installation instructions.

* **Version 1.0.6** September 27, 2014

    - More robust internal rotor code

* **Version 1.0.5** September 9, 2014

    - Several bug fixes and cleanups

* **Version 1.0.4** August 8, 2014

    - Improved ``tamkin-driver.py`` with improved documentation.
    - Added missing test file.

* **Version 1.0.3** August 7, 2014

    - Fix bug in load_rotscan_g03log. This function now also works with G09.
    - More output (rotors and tunneling correction)
    - Initial version of the TAMkin driver script

* **Version 1.0.2** August 5, 2014

    - Add some statistics to several plots
    - Make the output of dump_modes_gaussian compatible with Avogadro.

* **Version 1.0.1** April 2, 2014

    - Bug fix related to reading CP2K output files.
    - Bug fix for linear molecules when option ``treatment=Full()`` is used in
      the normal mode analysis.

* **Version 1.0.0** December 6, 2013

    - Initial release through github.
