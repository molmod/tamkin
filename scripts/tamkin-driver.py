#!/usr/bin/env python
# -*- coding: utf-8 -*-
# TAMkin is a post-processing toolkit for normal mode analysis, thermochemistry
# and reaction kinetics.
# Copyright (C) 2008-2012 Toon Verstraelen <Toon.Verstraelen@UGent.be>, An Ghysels
# <An.Ghysels@UGent.be> and Matthias Vandichel <Matthias.Vandichel@UGent.be>
# Center for Molecular Modeling (CMM), Ghent University, Ghent, Belgium; all
# rights reserved unless otherwise stated.
#
# This file is part of TAMkin.
#
# TAMkin is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# In addition to the regulations of the GNU General Public License,
# publications and communications based in parts on this program or on
# parts of this program are required to cite the following article:
#
# "TAMkin: A Versatile Package for Vibrational Analysis and Chemical Kinetics",
# An Ghysels, Toon Verstraelen, Karen Hemelsoet, Michel Waroquier and Veronique
# Van Speybroeck, Journal of Chemical Information and Modeling, 2010, 50,
# 1736-1750W
# http://dx.doi.org/10.1021/ci100099g
#
# TAMkin is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
#--
#!/usr/bin/env python
u'''\

The TAMkin driver script
########################

The TAMkin driver script (``tamkin-driver.py``) runs a standardized TAMkin
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

The following config files are read by the ``tamkin-driver.py`` script:

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
'''


import sys, os, numpy as np
from glob import glob
from tamkin import *
from molmod.periodic import periodic

# A few variables for debugging purposes only
gradient_threshold = 1e-2
v_threshold = 0.01


def load_cfg(fn):
    '''Read a config file and return contents as a dictionary.'''
    result = {}

    if not os.path.isfile(fn):
        return {}

    def parse(value):
        '''Convert a string into a useful datatype, depending on the value'''
        if value.lower() == 'true':
            return True
        if value.lower() == 'false':
            return False
        try:
            return int(value)
        except ValueError:
            pass
        try:
            return float(value)
        except ValueError:
            pass
        return value

    with open(fn) as f:
        for line in f:
            line = line[:line.find('#')].strip()
            if len(line) == 0:
                continue
            words = line.split()
            if len(words) < 2:
                raise ValueError('Error reading %s. Each non-empty line must contain at least two words: key value1 [value2 ...].')
            key = words[0]
            values = [parse(word) for word in words[1:]]
            if len(values) == 1:
                result[key] = values[0]
            else:
                result[key] = values
    return result


def get_atom_indexes(cfg, key):
    result = cfg.get(key)
    if result is not None:
        result = np.array(result)
        if cfg.get('fortran'):
            result -= 1
    return result


def get_pf(dn, temps):
    '''Construct a partition function from a molecule directory.

       **Arguments:**

       dn
            A molecule directory name

       temps
            A list of temperatures at which the properties of the molecule must be
            evaluated.
    '''
    print '  Loading partition function from', dn

    # A1) load the molecule
    mol_cfg = load_cfg('%s/molecule.cfg' % dn)
    molecule = None

    fn_fchk_freq = '%s/freq/gaussian.fchk' % dn
    if os.path.isfile(fn_fchk_freq):
        fn_fchk_sp = '%s/sp/gaussian.fchk' % dn
        if not os.path.isfile(fn_fchk_sp):
            fn_fchk_sp = None
        molecule = load_molecule_g03fchk(fn_fchk_freq, fn_fchk_sp)

    fn_outcar_freq = '%s/freq/OUTCAR' % dn
    if os.path.isfile(fn_outcar_freq):
        if molecule is not None:
            raise IOError('Cannot mix outputs from different codes.')
        fn_poscar_freq = '%s/freq/POSCAR' % dn
        if not os.path.isfile(fn_poscar_freq):
            raise IOError('Missing POSCAR file in %s/freq.' % dn)
        molecule = load_molecule_vasp(fn_poscar_freq, fn_outcar_freq,
                                      is_periodic=mol_cfg.get('periodic', True))
        fixed = load_fixed_vasp(fn_outcar_freq)
        if len(fixed) > 0:
            molecule = molecule.copy_with(fixed=fixed)

    fn_cp2k_freq = '%s/freq/cp2k.out' % dn
    if os.path.isfile(fn_cp2k_freq):
        if molecule is not None:
            raise IOError('Cannot mix outputs from different codes.')
        fn_cp2k_sp = '%s/sp/cp2k.out' % dn
        if not os.path.isfile(fn_cp2k_sp):
            raise IOError('Missing file %s/sp/cp2k.out' % dn)
        molecule = load_molecule_cp2k(fn_cp2k_sp, fn_cp2k_freq,
                                      is_periodic=mol_cfg.get('periodic', True))
        fixed = load_fixed_cp2k(fn_cp2k_freq)
        if len(fixed) > 0:
            molecule = molecule.copy_with(fixed=fixed)

    # A2) if present add a grimme correction
    fn_dftd3 = '%s/dftd3/dftd3.out' % dn
    if os.path.isfile(fn_dftd3):
        print '    Loading DFT-D3 corrections'
        disp_energy = load_dftd3(fn_dftd3)
        molecule = molecule.copy_with(energy=molecule.energy + disp_energy)

    # B) Load all rotors and keep directories of each rotor (dns_rotor).
    rotors = []
    dns_rotor = []

    # B1) load all rotors computed with Gaussian
    for dn_rotor in glob('%s/rotor_g_*' % dn):
        print '    Loading rotor', dn_rotor
        dns_rotor.append(dn_rotor)
        fn_log = '%s/gaussian.log' % dn_rotor
        # Load the config file
        rotor_cfg = load_cfg(os.path.join(dn_rotor, 'rotor.cfg'))
        # Load the rotational scan data from the Gaussian log file.
        rotor_scan = load_rotscan_g03log(fn_log, get_atom_indexes(rotor_cfg, 'top'))
        # Construct a Rotor object (solves Schrodinger equation etc.)
        rotor = Rotor(rotor_scan, molecule,
                      suffix=os.path.basename(dn_rotor)[6:],
                      rotsym=rotor_cfg.get('rotsym', 1),
                      even=rotor_cfg.get('even', False),
                      num_levels=rotor_cfg.get('num_levels', 50),
                      dofmax=rotor_cfg.get('dofmax', 5),
                      v_threshold=v_threshold)
        rotors.append(rotor)

    # B2) load all free rotors
    for dn_rotor in glob('%s/rotor_f_*' % dn):
        print '    Loading rotor', dn_rotor
        dns_rotor.append(dn_rotor)
        fn_cfg = '%s/rotor.cfg' % dn_rotor
        # Load the config file
        rotor_cfg = load_cfg(fn_cfg)
        # Construct a Rotor object (solves Schrodinger equation etc.)
        rotor_scan = RotScan(get_atom_indexes(rotor_cfg, 'dihed'),
                             molecule,
                             get_atom_indexes(rotor_cfg, 'top'))
        rotor = Rotor(rotor_scan, molecule,
                      suffix=os.path.basename(dn_rotor)[6:],
                      rotsym=rotor_cfg.get('rotsym', 1),
                      num_levels=rotor_cfg.get('num_levels', 50),
                      dofmax=rotor_cfg.get('dofmax', 5),
                      v_threshold=v_threshold)
        rotors.append(rotor)

    # B3) load all custom rotors
    for dn_rotor in glob('%s/rotor_c_*' % dn):
        print '    Loading rotor', dn_rotor
        dns_rotor.append(dn_rotor)
        fn_dat = '%s/rotor.dat' % dn_rotor
        # Load the potential data
        potential = np.loadtxt(fn_dat).T
        # Load the config file
        rotor_cfg = load_cfg(os.path.join(dn_rotor, 'rotor.cfg'))
        # Construct a Rotor object (solves Schrodinger equation etc.)
        rotor_scan = RotScan(get_atom_indexes(rotor_cfg, 'dihed'),
                             molecule,
                             get_atom_indexes(rotor_cfg, 'top'),
                             potential)
        rotor = Rotor(rotor_scan, molecule,
                      suffix=os.path.basename(dn_rotor)[6:],
                      rotsym=rotor_cfg.get('rotsym', 1),
                      even=rotor_cfg.get('even', False),
                      num_levels=rotor_cfg.get('num_levels', 50),
                      dofmax=rotor_cfg.get('dofmax', 5),
                      v_threshold=v_threshold)
        rotors.append(rotor)

    # C) Load fixed atoms and rigid blocks. The fixed atoms loaded here may override fixed
    #    atom settings in the output files of the quantum chemistry codes.
    fn_fixed = '%s/fixed.txt' % dn
    if os.path.isfile(fn_fixed):
        print '    Loading fixed atoms from %s/fixed.txt' % dn
        fixed = load_indices(fn_fixed)
        molecule = molecule.copy_with(fixed=fixed)

    blocks = None
    fn_blocks = '%s/blocks.txt' % dn
    if os.path.isfile(fn_blocks):
        print '    Loading rigid blocks from %s/blocks.txt' % dn
        blocks = load_indices(fn_blocks, groups=True)

    # D) Perform a normal mode analysis
    ext_dof = False
    if molecule.fixed is None:
        if blocks is None:
            if molecule.periodic:
                print '    Performing NMA on a periodic system, without fixed atoms'
                nma = NMA(molecule)
            else:
                print '    Performing NMA with fixed external degrees of freedom'
                nma = NMA(molecule, ConstrainExt(gradient_threshold=gradient_threshold))
                ext_dof = True
        else:
            print '    Performing NMA without %i mobile blocks' % len(blocks)
            nma = NMA(molecule, MBH(blocks))
    else:
        if blocks is None:
            print '    Performing NMA with %i fixed atoms (out of %i)' % (
                len(molecule.fixed), molecule.size)
            nma = NMA(molecule, PHVA(molecule.fixed))
        else:
            print '    Performing NMA with %i fixed atoms (out of %i) and %i mobile blocks' % (
                len(molecule.fixed), molecule.size, len(blocks))
            nma = NMA(molecule, PHVA_MBH(molecule.fixed, blocks))

    # E) Define the partition function
    pf_contributions = [Vibrations(freq_scaling=mol_cfg.get('freq_scaling', 1.0),
                                   zp_scaling=mol_cfg.get('zp_scaling', 1.0))]
    if ext_dof:
        pf_contributions.append(ExtTrans())
        pf_contributions.append(ExtRot(mol_cfg.get('symnum', None)))
    pf_contributions.extend(rotors)
    pf = PartFun(nma, pf_contributions)
    # Store the atomic numbers as an attribute of the partition function (used
    # to check that no atoms get lost in a reaction.
    pf.numbers = molecule.numbers

    # F) Make plots of the rotor
    for dn_rotor, rotor in zip(dns_rotor, rotors):
        rotor.plot_levels('%s/levels.png' % dn_rotor, 300)

    # G) Write a partf.txt file
    pf.write_to_file('%s/partf.txt' % dn)

    # H) Write a thermo analysis if temperatures are provided
    if len(temps) > 0:
        ta = ThermoAnalysis(pf, temps)
        ta.write_to_file('%s/thermo.csv' % dn)

    return pf


def get_chemical_formula(numbers):
    tmp = dict((number, (numbers == number).sum()) for number in np.unique(numbers))
    return ' '.join(['%s_%i' % (periodic[number].symbol, count) for number, count in tmp.iteritems()])


def check_mass_balance(pfs1, pfs2, name1, name2):
    cf1 = get_chemical_formula(np.concatenate([pf.numbers for pf in pfs1]))
    cf2 = get_chemical_formula(np.concatenate([pf.numbers for pf in pfs2]))
    if cf1 != cf2:
        print 'Error: mass balance mismatch!'
        print 'Chemical formula of the %s: %s' % (name1, cf1)
        for pf in pfs1:
            print '%30s:    %s' % (pf.title, get_chemical_formula(pf.numbers))
        print 'Chemical formula of the %s: %s' % (name1, cf2)
        for pf in pfs2:
            print '%30s:    %s' % (pf.title, get_chemical_formula(pf.numbers))
        print
        raise RuntimeError('Mass balance mismatch')


def main_kinetics(res, ts, prs):
    cfg = load_cfg('kinetics.cfg')

    if cfg.get('tunneling', False):
        # Define the tunneling correction
        tunneling = Eckart(res, ts, prs)
    else:
        tunneling = None

    km_trans = KineticModel(res, ts, tunneling)
    ra_trans = ReactionAnalysis(km_trans, cfg['temp_low'], cfg['temp_high'], cfg.get('temp_step', 10))
    ra_trans.plot_arrhenius("arrhenius.png")
    ra_trans.write_to_file("kinetics.txt")


def main_equilibrium(res, prs):
    cfg = load_cfg('equilibrium.cfg')

    tm = ThermodynamicModel(res, prs)
    tm.write_to_file("equilibrium.txt")
    for temp in cfg.get('temps', []):
        tm.write_table(temp, "thermo_%07.2f.csv" % temp)


def main():
    if len(sys.argv) > 1:
        raise ValueError('The TAMKin driver script takes no command line arguments.')

    temps = load_cfg('thermo.cfg').get('temps', [])

    reactants = []
    tss = []
    products = []
    for dn in glob('*/'):
        # chop trailing backslash
        dn = dn[:-1]
        if dn.startswith('re_'):
            reactants.append(get_pf(dn, temps))
        elif dn.startswith('ts_'):
            tss.append(get_pf(dn, temps))
        elif dn.startswith('pr_'):
            products.append(get_pf(dn, temps))
        else:
            print 'WARNING: skipping directory %s' % dn

    if len(reactants) == 0:
        raise RuntimeError('At least one reactant must be present.')
    if len(tss) > 1:
        raise RuntimeError('At most one transition state may be present.')

    if len(tss) == 1:
        check_mass_balance(reactants, tss, 'reactants', 'transition state')
        if len(products) > 0:
            check_mass_balance(reactants, products, 'reactants', 'products')
        main_kinetics(reactants, tss[0], products)
    if len(products) > 0:
        check_mass_balance(reactants, products, 'reactants', 'products')
        main_equilibrium(reactants, products)


if __name__ == '__main__':
    if len(sys.argv) > 1:
        print __doc__
    else:
        try:
            main()
        except:
            print 'An Error Occured! For detailed documentation, run ``tamkin-driver.py -h``.'
            print
            raise
