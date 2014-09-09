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
"""Convenient interfaces to define thermodynamic and kinetic models"""


import numpy as np, csv

from molmod import boltzmann, kjmol, second, meter, mol, planck

from tamkin.partf import PartFun


__all__ = [
    "BaseModel", "ThermodynamicModel", "BaseKineticModel",
    "KineticModel", "ActivationKineticModel"
]


class BaseModel(object):
    """Base class for all physico-chemical models."""
    def __init__(self):
        """
           Useful attribute:
            | ``pfs_all`` -- A dictionary with all partition functions involved
                             in the model. The values are the partition function
                             objects. The keys are the corresponding signed
                             stoichiometries.

           The methods in the base class are mainly used by the Monte Carlo
           routine in the ReactionAnalysis class.
        """
        self.pfs_all = {}
        self.pfs_list = []

    def _add_pfs(self, pfs, factor=1.0):
        """Add partition functions to the global list.

           Argument:
            | ``pfs`` -- A list of partition functions.

           Optional argument
            | ``factor`` -- The stoichiometry of all partition functions will be
                            multiplied with this factor.
        """
        for pf in pfs:
            if isinstance(pf, PartFun):
                add_st = 1.0
            else:
                # It may also be a tuple of a partition function and a
                # stoichiometry.
                pf, add_st = pf
            st = self.pfs_all.get(pf)
            if st is None:
                # This is used to maintain an ordered list of partition functions
                self.pfs_list.append(pf)
                st = 0.0
            self.pfs_all[pf] = st + factor*add_st

    def _set_unit(self, per_second=False):
        """Set the unit of the equilibrium constant or pre-exponential factor.

           Optional argument:
            | ``per_second`` -- Boolean, when True the unit is divided by
                                second.

           This functions is called by the constructor after all partition
           functions are added. (See :meth:`_add_pfs`.)
        """
        meter_power = 0
        mol_power = 0
        for pf, st in self._iter_pfs():
            if hasattr(pf, "translational"):
                meter_power -= pf.translational.dim*st
            mol_power += st
        unit = mol**mol_power*meter**meter_power
        unit_name = ""
        if meter_power != 0:
            if meter_power == 1:
                unit_name += "m"
            else:
                unit_name += "m^%i" % meter_power
        if mol_power != 0:
            if mol_power == 1:
                unit_name += " mol"
            else:
                unit_name += " mol^%i" % mol_power
        if per_second:
            unit /= second
            if len(unit_name) != 0:
                unit_name += " "
            unit_name += "s^-1"
        self.unit = unit
        self.unit_name = unit_name

    def _iter_pfs(self):
        """Iterate over all partition functions with their stoichiometries."""
        for pf in self.pfs_list:
            st = self.pfs_all[pf]
            if abs(st) > 0:
                yield pf, st

    def backup_freqs(self):
        """Keep a backup copy of the frequencies and the energy of each partition function."""
        for pf in self.pfs_all:
            pf.vibrational.positive_freqs_orig = pf.vibrational.positive_freqs.copy()
            pf.vibrational.negative_freqs_orig = pf.vibrational.negative_freqs.copy()
            pf.electronic.energy_backup = pf.electronic.energy

    def alter_freqs(self, freq_error, scale_energy):
        """Randomly distort the frequencies and energies.

           Arguments:
            | ``freq_error`` -- The absolute error to be introduced in the
                                frequencies.
            | ``scale_energy`` -- The relative error to be introduced in the
                                  (electronic) energies.
        """
        for pf in self.pfs_all:
            N = len(pf.vibrational.positive_freqs)
            freq_shift = np.random.normal(0, freq_error, N)
            pf.vibrational.positive_freqs = pf.vibrational.positive_freqs_orig + freq_shift
            pf.vibrational.positive_freqs[pf.vibrational.positive_freqs<=0] = 0.01
            N = len(pf.vibrational.negative_freqs)
            freq_shift = np.random.normal(0, freq_error, N)
            pf.vibrational.negative_freqs = pf.vibrational.negative_freqs_orig + freq_shift
            pf.vibrational.negative_freqs[pf.vibrational.negative_freqs>=0] = -0.01
            pf.electronic.energy = pf.electronic.energy_backup*scale_energy

    def restore_freqs(self):
        """Restore the backup of the frequencies and the energy of each partition function."""
        for pf in self.pfs_all:
            pf.vibrational.positive_freqs = pf.vibrational.positive_freqs_orig
            pf.vibrational.negative_freqs = pf.vibrational.negative_freqs_orig
            pf.electronic.energy = pf.electronic.energy_backup
            del pf.vibrational.positive_freqs_orig
            del pf.vibrational.negative_freqs_orig
            del pf.electronic.energy_backup

    def free_energy_change(self, temp):
        """Compute the change in free energy.

           The change in free energy depends on the the pressure in the
           ``ExtTrans`` contribution to the partition function, if such a
           contribution would be present.

           Argument:
            | ``temp``  -- the temperature
        """
        return sum(pf.chemical_potential(temp)*st for pf, st in self._iter_pfs())

    def energy_difference(self):
        """Compute the electronic energy difference between (+) products and (-) reactants."""
        return sum(pf.electronic.energy*st for pf, st in self._iter_pfs())

    def zero_point_energy_difference(self):
        """Compute the zero-point energy difference between (+) products and (-) reactants."""
        return sum(pf.zero_point_energy()*st for pf, st in self._iter_pfs())

    def internal_heat_difference(self, temp):
        """Compute the internal_heat difference between (+) products and (-) reactants.

           Argument:
            | ``temp`` -- The temperature.
        """
        return sum(pf.internal_heat(temp)*st for pf, st in self._iter_pfs())

    def equilibrium_constant(self, temp, do_log=False):
        """Compute the equilibrium constant at the given temperature.

           Argument:
            | ``temp`` -- The temperature.

           Optional argument:
            | ``do_log`` -- When True, the logarithm of the equilibrium constant
                            is returned instead of just the equilibrium constant
                            itself. [default=False]
        """
        log_K = sum(pf.logv(temp)*st for pf, st in self._iter_pfs())
        if do_log:
            return log_K
        else:
            return np.exp(log_K)

    def write_table(self, temp, filename):
        """Write a CSV file with the principal energies to a file.

           Arguments:
            | ``temp`` -- The temperature to use for the temperature-dependent
                          quantities.
            | ``filename`` -- The name of the CSV file.
        """
        f = file(filename, "w")
        c = csv.writer(f)
        self.dump_table(temp, c)
        f.close()


    def dump_table(self, temp, c):
        """Write a CSV file with the principal energies to a stream .

           Arguments:
            | ``temp`` -- The temperature to use for the temperature-dependent
                          quantities.
            | ``c`` -- A csv.writer object from the built-in Python csv module.
        """
        c.writerow(["Temperature [K]", temp])
        c.writerow([])
        c.writerow(["Quantity"] +
                   [pf.title for pf, st in self._iter_pfs()] +
                   ["Linear combination (always in kJ/mol)"])
        c.writerow(["Signed stoichiometry"] +
                   [st for pf, st in self._iter_pfs()])
        c.writerow(["**Values in a.u.**"])
        c.writerow(["Electronic energy"] +
                   [pf.electronic.energy for pf, st in self._iter_pfs()] +
                   [self.energy_difference()/kjmol])
        c.writerow(["Zero-point energy"] +
                   [pf.zero_point_energy() for pf, st in self._iter_pfs()] +
                   [self.zero_point_energy_difference()/kjmol])
        c.writerow(["Internal heat (%.2fK)" % temp] +
                   [pf.internal_heat(temp) for pf, st in self._iter_pfs()] +
                   [self.internal_heat_difference(temp)/kjmol])
        c.writerow(["Chemical potential (%.2fK)" % temp] +
                   [pf.chemical_potential(temp) for pf, st in self._iter_pfs()] +
                   [self.free_energy_change(temp)/kjmol])
        c.writerow(["**Corrections in kJ/mol**"])
        c.writerow(["Zero-point energy"] +
                   [(pf.zero_point_energy() - pf.electronic.energy)/kjmol for pf, st in self._iter_pfs()] +
                   [(self.zero_point_energy_difference() - self.energy_difference())/kjmol])
        c.writerow(["Internal heat (%.2fK)" % temp] +
                   [(pf.internal_heat(temp) - pf.electronic.energy)/kjmol for pf, st in self._iter_pfs()] +
                   [(self.internal_heat_difference(temp) - self.energy_difference())/kjmol])
        c.writerow(["Chemical potential (%.2fK)" % temp] +
                   [(pf.chemical_potential(temp) - pf.electronic.energy)/kjmol for pf, st in self._iter_pfs()] +
                   [(self.free_energy_change(temp) - self.energy_difference())/kjmol])
        c.writerow([])
        c.writerow(["**Other quantities**", "Unit", "Value"])

    def write_to_file(self, filename):
        """Write the model to a text file.

           One argument:
            | ``filename`` -- The file to write the output.
        """
        f = file(filename, "w")
        self.dump(f)
        f.close()

    def dump(self, f):
        """Write all info about the model to a file."""
        print >> f, "The chemical balance:"
        print >> f, "  ", " + ".join("%s*(\"%s\")" % (-st, pf.title) for pf, st in self._iter_pfs() if st < 0),
        print >> f, " <--> ",
        print >> f, " + ".join("%s*(\"%s\")" % (st, pf.title) for pf, st in self._iter_pfs() if st > 0)
        print >> f
        for counter, (pf, st) in enumerate(self._iter_pfs()):
            print >> f, "Partition function %i" % counter
            print >> f, "Signed stoichiometry: %i" % st
            pf.dump(f)
            print >> f


class ThermodynamicModel(BaseModel):
    """A model for a thermodynamic equilibrium."""
    def __init__(self, pfs_react, pfs_prod):
        """
           Arguments:
            | ``pfs_react`` -- A list with reactant partition functions.
            | ``pfs_prod`` -- A list with product partition functions.

           Both arguments are lists whose items should be PartFun objects.
           One may also replace an item by a ``(pf, st)`` tuple, where ``pf`` is
           the partition function and ``st`` is the stoichiometry.

           Useful attributes:
            | ``unit_name`` -- A string with the SI unit of the equilibrium
                               constant
            | ``unit`` -- The conversion factor to transform the equilibrium
                          constant into SI units
                          (``equilibrium_const/self.unit``)
        """
        BaseModel.__init__(self)
        self._add_pfs(pfs_react, -1)
        self._add_pfs(pfs_prod, +1)
        self._set_unit()

    def dump_table(self, temp, c):
        """Write a CSV file with the principal energies to a stream .

           Arguments:
            | ``temp`` -- The temperature to use for the temperature-dependent
                          quantities.
            | ``c`` -- A csv.writer object from the built-in Python csv module.
        """
        BaseModel.dump_table(self, temp, c)
        c.writerow(["Equilibrium constant", self.unit_name, self.equilibrium_constant(temp)/self.unit])

    def dump(self, f):
        """Write all info about the thermodynamic model to a file."""
        delta_E = self.energy_difference()
        print >> f, "Electronic energy difference [kJ/mol] = %.1f" % (delta_E/kjmol)
        delta_ZPE = self.zero_point_energy_difference()
        print >> f, "Zero-point energy difference [kJ/mol] = %.1f" % (delta_ZPE/kjmol)
        print >> f
        BaseModel.dump(self, f)


class BaseKineticModel(BaseModel):
    """A generic model for the rate constant of a chemical reaction."""

    def rate_constant(self, temp, do_log=False):
        """Compute the rate constant of the reaction in this analysis

           Arguments:
            | ``temp`` -- The temperature.

           Optional argument:
            | ``do_log`` -- When True, the logarithm of the rate constant is
                            returned instead of just the rate constant itself.
                            [default=False]
        """
        raise NotImplementedError


class KineticModel(BaseKineticModel):
    """A model for the rate constant of a single-step chemical reaction."""
    def __init__(self, pfs_react, pf_trans, tunneling=None):
        """
           Arguments:
            | ``pfs_react`` -- A list of partition functions for the reactants.
            | ``pf_trans`` -- The partition function of the transition state.

           The first argument is a list whose items should be PartFun objects.
           One may also replace an item by a ``(pf, st)`` tuple, where ``pf`` is
           the partition function and ``st`` is the stoichiometry.

           Optional argument:
            | ``tunneling`` -- A tunneling correction object. If not given, no
                               tunneling correction is applied.


           Useful attributes:
            | ``unit_name`` -- A string containing the SI unit of the rate
                               constant.
            | ``unit`` -- The conversion factor to transform the rate constant
                          to SI units (``rate_const/self.unit``)
        """
        BaseKineticModel.__init__(self)
        if len(pfs_react) == 0:
            raise ValueError("At least one reactant must be given.")
        self._add_pfs(pfs_react, -1)
        self._add_pfs([pf_trans], +1)
        self.tunneling = tunneling
        self._set_unit(per_second=True)

    def rate_constant(self, temp, do_log=False):
        """See :meth:`BaseKineticModel.rate_constant`

           The implementation is based on transition state theory.
        """
        result = self.equilibrium_constant(temp, do_log)
        if do_log:
            result = np.log(boltzmann*temp/planck) + result
            if self.tunneling is not None:
                result += np.log(self.tunneling(temp))
        else:
            result = boltzmann*temp/planck*result
            if self.tunneling is not None:
                result *= self.tunneling(temp)
        return result

    def dump_table(self, temp, c):
        """Write a CSV file with the principal energies to a stream .

           Arguments:
            | ``temp`` -- The temperature to use for the temperature-dependent
                          quantities.
            | ``c`` -- A csv.writer object from the built-in Python csv module.
        """
        BaseKineticModel.dump_table(self, temp, c)
        c.writerow(["Rate constant", self.unit_name, self.rate_constant(temp)/self.unit])

    def dump(self, f):
        """Write all info about the kinetic model to a file."""
        delta_E = self.energy_difference()
        print >> f, "Electronic energy barrier [kJ/mol] = %.1f" % (delta_E/kjmol)
        delta_ZPE = self.zero_point_energy_difference()
        print >> f, "Zero-point energy barrier [kJ/mol] = %.1f" % (delta_ZPE/kjmol)
        print >> f
        BaseKineticModel.dump(self, f)


class ActivationKineticModel(BaseKineticModel):
    """A model for the rate constant of a single-step chemical reaction with a pre-reactive complex."""
    def __init__(self, tm, km):
        """
           Arguments:
            | ``tm`` -- The thermodynamic model for the pre-reactive complex.
            | ``km`` -- The kinetic model for the single-step reaction.
        """
        BaseKineticModel.__init__(self)
        self.tm = tm
        self.km = km
        self._add_pfs(tm.pfs_all.iteritems())
        self._add_pfs(km.pfs_all.iteritems())
        self._set_unit(per_second=True)

    def rate_constant(self, temp, do_log=False):
        """See :meth:`BaseKineticModel.rate_constant`"""
        if do_log:
            return self.tm.equilibrium_constant(temp, True) + \
                   self.km.rate_constant(temp, True)
        else:
            return self.tm.equilibrium_constant(temp, False)* \
                   self.km.rate_constant(temp, False)

    def dump(self, f):
        """Write all info about the kinetic model to a file."""
        delta_E = self.energy_difference()
        print >> f, "Electronic energy barrier [kJ/mol] = %.1f" % (delta_E/kjmol)
        delta_ZPE = self.zero_point_energy_difference()
        print >> f, "Zero-point energy barrier [kJ/mol] = %.1f" % (delta_ZPE/kjmol)
        print >> f
        print >> f, "Thermodynamic submodel"
        self.tm.dump(f)
        print >> f
        print >> f, "Kinetic submodel"
        self.km.dump(f)
