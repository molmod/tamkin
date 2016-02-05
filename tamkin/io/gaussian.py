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


from tamkin.data import Molecule, RotScan

from molmod.io import FCHKFile
from molmod import dihed_angle, amu, angstrom

import numpy as np


__all__ = [
    "load_fixed_g03com", "load_molecule_g03fchk", "load_molecule_g98fchk",
    "load_rotscan_g03log",
]


def load_fixed_g03com(filename):
    """Load fixed atoms from a gaussian input file.

       Argument:
        | ``filename`` -- The gaussian input file.

       A fixed atom is recognized by the '-1' after the atom symbol in the
       molecule specification. The '-1' is the second word in the line.
    """
    f = file(filename)
    for line in f:
        # iterate until we reach the title line
        line = line.strip()
        if not (len(line)==0 or line[0] == "#" or line[0] == "%"):
            break
    for line in f:
        # skip lines until charge/multiplicty is reached
        words = line.split()
        if len(words) == 2:
            try:
                int(words[0])
                int(words[1])
                break
            except ValueError:
                pass
    counter = 0
    fixed_atoms = []
    for line in f:
        # go trough the list of atoms and store the fixed ones.
        words = line.split()
        if len(words) == 0:
            break
        if words[1] == "-1":
            fixed_atoms.append(counter)
        counter += 1
    return fixed_atoms


def iter_floats_file(fn):
    f = file(fn)
    for line in f:
        words = line.split()
        for w in words:
            yield float(w.replace('D', 'E'))
    f.close()


def load_molecule_g03fchk(fn_freq, fn_ener=None, fn_vdw=None, energy=None, fn_punch=None):
    """Load a molecule from Gaussian03 formatted checkpoint files.

       Arguments:
         | ``fn_freq`` -- the formatted checkpoint file of the frequency job

       Optional arguments:
         | ``fn_ener`` -- the formatted checkpoint file of a single point
                          computation for the energy. When not given, the energy
                          is taken from the frequency job.
         | ``fn_vdw`` -- An orca output file containing a Van der Waals
                         correction for the energy.
         | ``energy`` -- Override the energy from the formatted checkpoint file
                         with the given value.
         | ``punch`` -- A Gaussian derivatives punch file. When given, the
                        gradient and the Hessian are read from this file
                        instead.
    """

    fchk_freq = FCHKFile(fn_freq, ignore_errors=True, field_labels=[
        "Cartesian Force Constants", "Real atomic weights", "Total Energy",
        "Multiplicity", "Cartesian Gradient", "MicOpt",
    ])
    if fn_ener is None:
        fchk_ener = fchk_freq
    else:
        fchk_ener = FCHKFile(fn_ener, ignore_errors=True, field_labels=[
            "Total Energy"
        ])
    if energy is None:
        energy = fchk_ener.fields["Total Energy"]
    vdw = 0
    if fn_vdw is not None:
        from tamkin.io.dispersion import load_dftd_orca
        vdw = load_dftd_orca(fn_vdw)

    natom = fchk_freq.molecule.size
    if fchk_freq.molecule.size == 1 and \
       "Cartesian Force Constants" not in fchk_freq.fields:
        gradient = np.zeros((1,3), float)
        hessian = np.zeros((3,3), float)
    elif fn_punch is None:
        gradient = fchk_freq.fields["Cartesian Gradient"].copy()
        gradient.shape = (natom, 3)
        hessian = fchk_freq.get_hessian()
    else:
        gradient = np.zeros((natom, 3), float)
        hessian = np.zeros((3*natom, 3*natom), float)
        iterator = iter_floats_file(fn_punch)
        for i in xrange(natom):
            for j in xrange(3):
                gradient[i,j] = iterator.next()
        for i in xrange(3*natom):
            for j in xrange(i+1):
                v = iterator.next()
                hessian[i,j] = v
                hessian[j,i] = v

    if "MicOpt" in fchk_freq.fields:
        fixed = (fchk_freq.fields["MicOpt"] == -2).nonzero()[0]
        if len(fixed) == 0:
            fixed = None
    else:
        fixed = None

    return Molecule(
        fchk_freq.molecule.numbers,
        fchk_freq.molecule.coordinates,
        fchk_freq.fields["Real atomic weights"]*amu,
        energy+vdw,
        gradient,
        hessian,
        fchk_freq.fields["Multiplicity"],
        None, # gaussian is very poor at computing the rotational symmetry number
        False,
        title=fchk_freq.title,
        fixed=fixed,
    )


g98_masses = np.array([
    1.0079, 4.0026, 6.94, 9.01218, 10.81, 12.011, 14.0067, 15.9994, 18.9984,
    20.179, 22.98977, 24.305, 26.98154, 28.0855, 30.97376, 32.06, 35.453,
    39.948, 39.0983, 40.08, 44.9559, 47.9, 50.9415, 51.996, 54.938, 55.847,
    58.9332, 58.71, 63.546, 65.38, 69.735, 72.59, 74.9216, 78.96, 79.904, 83.8,
    85.4678, 87.62, 88.9059, 91.22, 92.9064, 95.94, 98.9062, 101.07, 102.9055,
    106.4, 107.868, 112.41, 114.82, 118.69, 121.75, 127.6, 126.9045, 131.3,
    132.9054, 137.33, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 178.49, 180.9479, 183.85, 186.207, 190.2, 192.22,
    195.09, 196.9665, 200.59, 204.37, 207.2, 208.9804, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0,
])*amu


def load_molecule_g98fchk(fn_freq, fn_ener=None, energy=None):
    """Load a molecule from Gaussian98 formatted checkpoint files.

       Arguments:
         | ``fn_freq`` -- The formatted checkpoint file of the frequency job.

       Optional arguments:
         | ``fn_ener`` -- The formatted checkpoint file of a single point
                          computation for the energy. When not given, the energy
                          is taken from the frequency job.
         | ``energy`` -- Override the energy from the formatted checkpoint file
                         with the given value.
    """

    fchk_freq = FCHKFile(fn_freq, ignore_errors=True, field_labels=[
        "Cartesian Force Constants", "Total Energy",
        "Multiplicity", "Cartesian Gradient"
    ])
    if fn_ener is None:
        fchk_ener = fchk_freq
    else:
        fchk_ener = FCHKFile(fn_ener, ignore_errors=True, field_labels=[
            "Total Energy"
        ])
    masses = np.array([g98_masses[n-1] for n in fchk_freq.molecule.numbers])
    if energy is None:
        energy = fchk_ener.fields["Total Energy"]

    return Molecule(
        fchk_freq.molecule.numbers,
        fchk_freq.molecule.coordinates,
        masses,
        energy,
        np.reshape(np.array(fchk_freq.fields["Cartesian Gradient"]), (len(fchk_freq.molecule.numbers),3)),
        fchk_freq.get_hessian(),
        fchk_freq.fields["Multiplicity"],
        None, # gaussian is very poor at computing the rotational symmetry number
        False,
    )


def load_rotscan_g03log(fn_log, top_indexes=None):
    """Load the torsional potential from a Gaussian 03 log/output file.

       Argument:
        | ``fn_log`` -- The filename of the gaussian output.

       Optional argument:
        | ``top_indexes`` -- The atom indexes that define the rotor. These do
                             not have to include the atoms that define the
                             rotational axis. When not given, an attempt is made
                             to derive this information from the dihedral angle
                             used for the scan.
    """
    # find the line that specifies the dihedral angle
    with open(fn_log) as f:
        found_modredundant = False
        for line in f:
            if line.startswith(" The following ModRedundant input section has been read:"):
                found_modredundant = True
                break
        if not found_modredundant:
            raise IOError("Could not find the ModRedundant section in the log file.")

        dihedral = None
        num_geoms = 0
        for line in f:
            line = line.strip()
            if len(line) == 0:
                break
            else:
                words = line.split()
                if words[0] == "D" and (words[5] == "S" or words[5] == "F"):
                    if dihedral is None:
                        dihedral = list(int(word)-1 for word in words[1:5])
                    else:
                        raise IOError("Found multiple dihedral angle scan, which is not supported.")

        if dihedral is None:
            raise IOError("Could not find the dihedral angle of the rotational scan.")

        # load all the energies and compute the corresponding angles
        energies = []
        angles = []
        geometries = []
        for line in f:
            if line.startswith("                          Input orientation:") or \
               line.startswith("                         Standard orientation:"):
                # read the molecule
                numbers = []
                last_coordinates = []
                ## skip four lines
                for i in xrange(4):
                    f.next()
                # read atoms
                for line in f:
                    if line.startswith(" -----"):
                        break
                    words = line.split()
                    numbers.append(int(words[1]))
                    last_coordinates.append((float(words[3]), float(words[4]), float(words[5])))
            if line.startswith(" SCF Done:"):
                # read the energy
                last_energy = float(line[line.find("=")+1:].split()[0])
            if line.startswith("    -- Stationary point found."):
                # store last emergy and geometry in list
                energies.append(last_energy)
                last_coordinates = np.array(last_coordinates)*angstrom
                geometries.append(last_coordinates)
                angles.append(dihed_angle(last_coordinates[dihedral])[0])

        if len(energies) == 0:
            raise IOError("Could not find any stationary point")

        if top_indexes is None:
            # Define the molecular geometry that is used in the constructor of
            # RotScan to detect the top.
            from molmod.molecules import Molecule as BaseMolecule
            molecule = BaseMolecule(numbers, geometries[0])
        else:
            molecule = None
        result = RotScan(
            dihedral, molecule, top_indexes,
            np.array([angles, energies])
        )
        result.geometries = np.array(geometries)
        return result
