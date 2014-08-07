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


import os, numpy as np
from glob import glob
from tamkin import *


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


def get_pf(dn):
    '''Construct a partition function from a standard directory layout.

       **Arguments:**

       dn
            A directory name

       The following file must be present in the directory::

            freq/gaussian.fchk

       The following files may be present in the directory::

            molecule.cfg           # to specify additional parameters of the
                                   # molecule. (Details given below.)

            dftd3/dftd3.out        # to add an a posteriori Grimme Dispersion
                                   # correction.

            sp/gaussian.fchk       # to replace the energy from the frequency
                                   # job by a more accurate one.

            rotor_g_*/gaussian.log # to add a rotor based on a Gaussian relaxed
            rotor_g_*/rotor.cfg    # scan. (More than one allowed.)

            rotor_f_*/rotor.cfg    # to specify a free rotor. (Details given
                                   # below. More than one allowed.)

            rotor_c_*/rotor.dat    # to specify a custom hindered rotor.
            rotor_c_*/rotor.cfg    # (Details given below. More than one
                                   # allowed.)

       All other files are simply ignored.

       All configuration files (*.cfg) have the following format. Every
       non-empty lines consists of a key and one or more values. Comments can
       be added with a #, just as in Python source code. Each file has its
       specific keys that are processed. Unknown keys are ignored.

       **molecule.cfg:**

       ``symnum``
            This keywords can be used to assign the rotational symmetry number.
            For molecules with less than 10 atoms, this number is estimated
            automatically when not given. For larger molecules, the default
            value is 1.

       **rotor_g_*/rotor.cfg**

       ``rotsym``
            The rotational symmetry of the internal rotor. [default=1]

       ``even``
            A boolean (True of False) to indicate that the torsional potential
            is even. [default=False]

       ``num_levels``
            The number of energy levels considered in the QM treatment of the
            rotor. [default=50]

       ``dofmax``
            The maximum number of cosines used to represent the torsional
            potential. if the potential is not even, the same number of sines is
            also used. [default=5]

       **rotor_f_*/rotor.cfg**

       ``rotsym``
            The rotational symmetry of the internal rotor. [default=1]

       ``num_levels``
            The number of energy levels considered in the QM treatment of the
            rotor. [default=50]

       ``dofmax``
            The maximum number of cosines used to represent the torsional
            potential. if the potential is not even, the same number of sines is
            also used. [default=5]

       **rotor_c_*/rotor.dat**

       This file just contains two colums of data, angles (rad) and energies
       (hartree), that specify the custom potential. It does not follow the cfg
       format.

       **rotor_c_*/rotor.cfg**

       ``dihed``
            A list of four atom indexes that define the dihedral angle,
            separated by whitespace.

       ``fortran``
            A boolean (True of False) to indicate that the atom indexes are
            given in Fortran convention (starting by one). [default=False]

       ``rotsym``
            The rotational symmetry of the internal rotor. [default=1]

       ``even``
            A boolean (True of False) to indicate that the torsional potential
            is even. [default=False]

       ``num_levels``
            The number of energy levels considered in the QM treatment of the
            rotor. [default=50]

       ``dofmax``
            The maximum number of cosines used to represent the torsional
            potential. if the potential is not even, the same number of sines is
            also used. [default=5]

    '''
    # A1) load the molecule
    molecule = load_molecule_g03fchk('%s/freq/gaussian.fchk' % dn)

    # A2) if present, load a more accurate energy
    fn_sp = '%s/sp/gaussian.fchk' % dn
    if os.path.isfile(fn_sp):
        molecule_energy = load_molecule_g03fchk(fn_sp)
        molecule = molecule.copywith(energy=molecule_energy.energy)

    # A3) if present add a grimme correction
    fn_dftd3 = '%s/dftd3/dftd3.out' % dn
    if os.path.isfile(fn_dftd3):
        disp_energy = load_dftd3(fn_dftd3)
        molecule = molecule.copywith(energy=molecule.energy + disp_energy)

    # B) Load all rotors and keep directories of each rotor (dns_rotor).
    rotors = []
    dns_rotor = []

    # B1) load all rotors computed with Gaussian
    for fn_log in glob('%s/rotor_g_*/gaussian.log' % dn):
        dn_rotor = os.path.dirname(fn_log)
        dns_rotor.append(dn_rotor)
        # Load the rotational scan data from the Gaussian log file.
        rotor_scan = load_rotscan_g03log(fn_log)
        # Load the config file
        rotor_cfg = load_cfg(os.path.join(dn_rotor, 'rotor.cfg'))
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
    for fn_cfg in glob('%s/rotor_f_*/rotor.cfg' % dn):
        dn_rotor = os.path.dirname(fn_cfg)
        dns_rotor.append(dn_rotor)
        # Load the config file
        rotor_cfg = load_cfg(fn_cfg)
        # Construct a Rotor object (solves Schrodinger equation etc.)
        rotor_scan = RotScan(rotor_cfg['atoms'], molecule)
        rotor = Rotor(rotor_scan, molecule,
                      suffix=os.path.basename(dn_rotor)[6:],
                      rotsym=rotor_cfg.get('rotsym', 1),
                      num_levels=rotor_cfg.get('num_levels', 50),
                      dofmax=rotor_cfg.get('dofmax', 5),
                      v_threshold=v_threshold)
        rotors.append(rotor)

    # B3) load all custom rotors
    for fn_dat in glob('%s/rotor_c_*/rotor.dat' % dn):
        dn_rotor = os.path.dirname(fn_dat)
        dns_rotor.append(dn_rotor)
        # Load the potential data
        potential = np.loadtxt(fn_dat).T
        # Load the config file
        rotor_cfg = load_cfg(os.path.join(dn_rotor, 'rotor.cfg'))
        # Construct a Rotor object (solves Schrodinger equation etc.)
        dihed = np.array(rotor_cfg['dihed']) - rotor_cfg.get('fortran', False)
        rotor_scan = RotScan(dihed, molecule, potential=potential)
        rotor = Rotor(rotor_scan, molecule,
                      suffix=os.path.basename(dn_rotor)[6:],
                      rotsym=rotor_cfg.get('rotsym', 1),
                      even=rotor_cfg.get('even', False),
                      num_levels=rotor_cfg.get('num_levels', 50),
                      dofmax=rotor_cfg.get('dofmax', 5),
                      v_threshold=v_threshold)
        rotors.append(rotor)

    # C) Perform a normal mode analysis
    nma = NMA(molecule, ConstrainExt(gradient_threshold=gradient_threshold))

    # D) Define the partition function
    mol_cfg = load_cfg('%s/molecule.cfg')
    pf = PartFun(nma, [ExtTrans(), ExtRot(mol_cfg.get('symnum', 1))] + rotors)

    # E) Make plots of the rotor
    for dn_rotor, rotor in zip(dns_rotor, rotors):
        rotor.plot_levels('%s/levels.png' % dn_rotor, 300)

    return pf


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
    tm = ThermodynamicModel(res, prs)
    tm.write_to_file("equilibrium.txt")


def main():
    '''Runs a standardized TAMkin computation.

       This script assumes that the current directory has subdirectories as
       follows::

            re_*/   # reactant subdirectories (at least one)
            ts_*/   # on transition state subdirectory (optional, at most one)
            pr_*/   # reaction product subdirectories (optional)

       The layout of each directory is documented in the function get_pf. At
       least one transition state or product must be present.

       When a transition state is present, the kinetic parameters are computed.
       In this case, a config file ``kinetics.cfg`` must be present. It contains
       the following keys:

       ``temp_low``, ``temp_high`` (mandatory)
            These specifiy the minimum and maximum temperature for the Arrhenius
            plot.

       ``temp_int``
            The temperature interval for the datapoints for the Arrhenius plot.
            [default=10]

       When reaction products are present, the equilibrium constant is computed.
    '''
    reactants = []
    tss = []
    products = []
    for dn in glob('*/'):
        # chop trailing backslash
        dn = dn[:-1]
        if dn.startswith('re_'):
            reactants.append(get_pf(dn))
        elif dn.startswith('ts_'):
            tss.append(get_pf(dn))
        elif dn.startswith('pr_'):
            products.append(get_pf(dn))
        else:
            print 'Warning: skipping directory %s' % dn

    if len(reactants) == 0:
        raise RuntimeError('At least one reactant must be present.')
    if len(tss) > 1:
        raise RuntimeError('At most one transition state may be present.')
    if len(tss) == 0 and len(products) == 0:
        raise RuntimeError('At least a transition state or one product must be present.')

    if len(tss) == 1:
        main_kinetics(reactants, tss[0], products)
    if len(products) > 0:
        main_equilibrium(reactants, products)


if __name__ == '__main__':
    main()
