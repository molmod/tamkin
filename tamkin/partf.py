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
"""Partition functions based on the harmonic oscillator approximation.

   The workhorse of this module is the **PartFun** class. A PartFun object
   represents a partition function with a solid interface. PartFun objects
   can be used to study chemical equilibrium, rate coefficients and various
   thermodynamic properties.

   These are the other classes and functions present in this module:

   * **Abstract classes:**
       * Info
       * StatFys
       * StatFysTerms
   * **Contributions to the partition function:**
       * Electronic
       * ExtTrans
       * ExtRot
       * Vibrations
       * Rotor (see rotor.py)
   * **Helper functions:**
       * helper_levels, helpert_levels, helpertt_levels,
       * helper_vibrations, helpert_vibrations, helpertt_vibrations

   **Important**: Partition functions can be constructed for NpT gases, NVT
   gases and many other systems. The return values of methods such as
   ``free_energy``, ``internal_heat`` and ``heat_capacity`` are often given
   specialized names in the context of different partition functions. For
   example, chemists tend to use the following names:

   * 3D NVT gas:
       - ``PartFun.free_energy`` -> the Helmholtz free energy
       - ``PartFun.internal_heat`` -> the internal energy
       - ``PartFun.heat_capacity`` -> the heat capacity at constant volume

   * 3D NpT gas:
       - ``PartFun.free_energy`` -> the Gibbs free energy
       - ``PartFun.internal_heat`` -> the enthalpy
       - ``PartFun.heat_capacity`` -> the heat capacity at constant pressure

   Don't say we did not warn you. Terminology can be very confusing.

   All the extensive thermodynamic quantities computed here are in atomic units
   per molecule. If you want to express the Gibbs free energy at 300 Kelvin of a
   system in kJ/mol, use the molmod module to perform unit conversions. For
   example::

     >>> pf = PartFun(..., [ExtTrans(cp=True)])
     >>> print pf.free_energy(300)/kjmol
"""


from molmod import boltzmann, lightspeed, atm, bar, amu, centimeter, kjmol, \
    planck, mol, meter, newton

import numpy as np


__all__ = [
    "Info", "StatFys", "StatFysTerms",
    "helper_levels", "helpert_levels", "helpertt_levels",
    "Electronic", "ExtTrans", "ExtRot", "PCMCorrection",
    "Vibrations",
    "helper_vibrations", "helpert_vibrations", "helpertt_vibrations",
    "PartFun",
]


class Info(object):
    """An object that has a name and that can dump info to a file."""
    def __init__(self, name):
        """
           Arguments:
            | ``name`` -- the name used for this object in the output
        """
        self.name = name

    def dump(self, f):
        """Write a description to file.

           Arguments:
            | ``f`` -- the file object to write to
        """
        print >> f, "  %s" % self.name.upper()

    def dump_values(self, f, label, values, format, num_col=8):
        """Write a nicely formatted array of numbers to file.

           Arguments:
            | ``f`` -- the file object to write to
            | ``label`` -- a label that explains the meaning of the numbers
            | ``values`` -- the array with numbers
            | ``format`` -- a Python format string for one number

           Optional argumet:
            | ``num_col`` -- the number of columns [default=8]
        """
        parts = ["    "]
        for counter, value in enumerate(values):
            parts.append(format % value)
            if counter % num_col == num_col-1:
                parts.append("\n    ")
        if len(parts) > 1:
            if parts[-1] == "\n    ":
                parts.pop()
            parts.insert(0, "    %s:\n" % label)
            print >> f, "".join(parts)


class StatFys(object):
    """Abstract class for (contributions to) the parition function.

       The constructor (__init__) and four methods (init_part_fun, helper,
       helpert, helpertt) must be implemented in derived classes.
    """
    def init_part_fun(self, nma, partf):
        """Compute parameters that depend on nma and partition function.

           Arguments:
            | ``nma`` -- an NMA object
            | ``partf`` -- A PartFun object

           This method is called by the PartFun object and should not be called
           by the user. When the PartFun object is initialized, each
           contribution can further initialize its parameters based on the
           information available in the nma and partf objects.
        """
        pass

    def helper(self, temp, n):
        r"""Helper function.

           Returns

           .. math:: T^n \frac{\ln(Z_N)}{N}

           where :math:`Z_N` is (the contribution to) the many body partition
           function and :math:`N` is the total number of particles. In all
           cases, except for the translational contribution, this comes down to:

           .. math:: T^n \ln(Z_1)

           where :math:`\ln(Z_1)` is the single-particle contribution to the
           partition function.

           Arguments:
            | ``temp`` -- the temperature
            | ``n`` -- the power for the temperature factor
        """
        raise NotImplementedError

    def helpert(self, temp, n):
        r"""Helper function T.

           Returns

           .. math:: T^n \frac{d \left(\frac{\ln(Z_N)}{N}\right)}{dT}

           where :math:`Z_N` is (the contribution to) the many body partition
           function and :math:`N` is the total number of particles.  In all
           cases, except for the translational contribution, this comes down to:

           .. math:: T^n \frac{\partial \ln(Z_1)}{\partial T}

           where :math:`\ln(Z_1)` is the single-particle contribution to the
           partition function.

           Arguments:
            | ``temp`` -- the temperature
            | ``n`` -- the power for the temperature factor
        """
        raise NotImplementedError

    def helpertt(self, temp, n):
        r"""Helper function TT.

           Returns

           .. math:: T^n \frac{d^2 \left(\frac{\ln(Z_N)}{N}\right)}{dT^2}

           where :math:`Z_N` is (the contribution to) the many body partition
           function and :math:`N` is the total number of particles. In all
           cases, except for the translational contribution, this comes down to:

           .. math:: T^n \frac{\partial^2 \ln(Z_1)}{\partial T^2}

           where :math:`\ln(Z_1)` is the single-particle contribution to the
           partition function.

           Arguments:
            | ``temp`` -- the temperature
            | ``n`` -- the power for the temperature factor
        """
        raise NotImplementedError

    def helpern(self, temp, n):
        r"""Helper function N.

           Returns

           .. math:: T^n \frac{\partial \ln(Z_N)}{\partial N}

           where :math:`Z_N` is (the contribution to) the many body partition
           function and :math:`N` is the total number of particles. This is used
           to compute the chemical potential and reaction free energies. In all
           cases, except for the total partition function and translational
           contribution, this comes down to the method ``helper``.

           Arguments:
            | ``temp`` -- the temperature
            | ``n`` -- the power for the temperature factor
        """
        # By default return helper
        return self.helper(temp, n)

    def helperv(self, temp, n):
        r"""Helper function V

           This always the same as the method ``helpern``, except for the
           total partition function and the translational contribution to the
           partition function. Then it returns

           .. math:: T^n \frac{\partial \ln(Z_N)}{\partial N} - \ln\left(\frac{V}{N}\right)

           This is used for the computation of equilibrium constants, and rate
           constants.

           Arguments:
            | ``temp`` -- temperature
            | ``n`` -- the power for the temperature factor
        """
        # By default return helper
        return self.helper(temp, n)

    def log(self, temp, helper=None):
        r"""Log function

           The logarithm of the N-particle partition function divided by the
           number of particles:

           .. math:: \frac{\ln(Z_N)}{N}

           For all contributions, except for the translational, this comes down
           to:

           .. math:: \ln(Z_1)

           Argument:
            | ``temp`` -- the temperature

           Optional argument:
            | ``helper`` -- an alternative implementation of helper
                            [default=self.helper]
        """
        if helper is None:
            helper = self.helper
        return helper(temp, 0)

    def logt(self, temp, helpert=None):
        r"""Log function T

           The derivative towards temperature of the logarithm of the N-particle
           partition function divided by the number of particles:

           .. math:: \frac{d \left(\frac{\ln(Z_N)}{N}\right)}{dT}

           For all contributions, except for the translational, this comes down
           to:

           .. math:: \frac{\partial \ln(Z_1)}{\partial T}

           Argument:
            | ``temp`` -- the temperature

           Optional arguments:
            | ``helpert`` -- an alternative implementation of helpert
                             [default=self.helpert]
        """
        if helpert is None:
            helpert = self.helpert
        return helpert(temp, 0)

    def logtt(self, temp, helpertt=None):
        r"""Log function TT

           The second derivative towards temperature of the logarithm of the
           N-particle partition function divided by the number of particles.

           .. math:: \frac{d^2 \left(\frac{\ln(Z_N)}{N}\right)}{dT^2}

           For all contributions, except for the translational, this comes down
           to:

           .. math:: \frac{\partial^2 \ln(Z_1)}{\partial T^2}

           Argument:
            | ``temp`` -- the temperature

           Optional arguments:
            | ``helpertt`` -- an alternative implementation of helpertt
                              [default=self.helpertt]
        """
        if helpertt is None:
            helpertt = self.helpertt
        return helpertt(temp, 0)

    def logn(self, temp, helpern=None):
        r"""Log function N

           The derivative of the logarithm of the many-particle partition
           function towards to number of particles.

           .. math:: \frac{\partial \ln(Z_N)}{\partial N}

           In most cases, this is the same as the logarithm of the
           single-particle partition divided by the number of particles. (see
           method ``log``.) The only exception is the translational contribution
           to the partition function, and hence also the total partition
           function.

           Argument:
            | ``temp`` -- the temperature

           Optional argument:
            | ``helpern`` -- an alternative implementation of helpern
                             [default=self.helpern]
        """
        if helpern is None:
            helpern = self.helpern
        return helpern(temp, 0)

    def logv(self, temp, helperv=None):
        r"""Log function V

           This always the same as the method ``logn``, except for the
           total partition function and the translational contribution to the
           partition function. Then it returns

           .. math:: \frac{\partial \ln(Z_N)}{\partial N} - \ln\left(\frac{V}{N}\right)

           This is used for the computation of equilibrium constants, and rate
           constants.

           Argument:
            | ``temp`` -- the temperature

           Optional argument:
            | ``helperv`` -- an alternative implementation of helperv
                             [default=self.helperv]
        """
        if helperv is None:
            helperv = self.helperv
        return helperv(temp, 0)

    def internal_heat(self, temp, helpert=None):
        """Computes the internal heat per molecule.

           Argument:
            | ``temp`` -- the temperature

           Optional argument:
            | ``helpert`` -- an alternative implementation of helpert
                             [default=self.helpert]
        """
        if helpert is None:
            helpert = self.helpert
        return boltzmann*helpert(temp, 2)

    def heat_capacity(self, temp, helpert=None, helpertt=None):
        """Computes the heat capacity per molecule.

           Argument:
            | ``temp`` -- the temperature

           Optional arguments:
            | ``helpert`` -- an alternative implementation of helpert
                             [default=self.helpert]
            | ``helpertt`` -- an alternative implementation of helpertt
                              [default=self.helpertt]
        """
        if helpert is None:
            helpert = self.helpert
        if helpertt is None:
            helpertt = self.helpertt
        return boltzmann*(2*helpert(temp, 1) + helpertt(temp, 2))

    def entropy(self, temp, helper=None, helpert=None):
        """Computes the entropy contribution per molecule.

           Argument:
            | ``temp`` -- the temperature

           Optional arguments:
            | ``helper`` -- an alternative implementation of helper
                            [default=self.helper]
            | ``helpert`` -- an alternative implementation of helpert
                             [default=self.helpert]
        """
        if helper is None:
            helper = self.helper
        if helpert is None:
            helpert = self.helpert
        return boltzmann*(helper(temp, 0) + helpert(temp, 1))

    def free_energy(self, temp, helper=None):
        """Computes the free energy per molecule.

           Argument:
            | ``temp`` -- the temperature

           Optional argument:
            | ``helper`` -- an alternative implementation of helper
                            [default=self.helper]
        """
        if helper is None:
            helper = self.helper
        return -boltzmann*helper(temp, 1)

    def zero_point_energy(self, helpern=None):
        """Return the zero-point energy.

           Optional argument:
            | ``helper`` -- an alternative implementation of helpern
                            [default=self.helpern]

           In TAMkin the zero-point energy is defined as the limit of the
           chemical potential for the temperature going towards zero.
        """
        return self.chemical_potential(0, helpern)

    def chemical_potential(self, temp, helpern=None):
        """Computes the chemical potential.

           Argument:
            | ``temp`` -- the temperature

           Optional argument:
            | ``helper`` -- an alternative implementation of helpern
                            [default=self.helpern]
        """
        if helpern is None:
            helpern = self.helpern
        return -boltzmann*helpern(temp, 1)


class StatFysTerms(StatFys):
    """Abstract class for (contributions to) the parition function with multiple terms.

       The different terms (or factors if you like) are of the same mathematical
       structure.

       The constructor (__init__) and the four methods (init_part_fun,
       helper_terms, helpert_terms, helpertt_terms) must be implemented in
       derived classes.
    """
    def __init__(self, num_terms):
        """
           Arguments:
             num_terms`` -- the number of terms present in this contribution.
        """
        self.num_terms = num_terms

    def helper(self, temp, n):
        """See :meth:`StatFys.helper`."""
        return self.helper_terms(temp, n).sum()

    def helpert(self, temp, n):
        """See :meth:`StatFys.helpert`."""
        return self.helpert_terms(temp, n).sum()

    def helpertt(self, temp, n):
        """See :meth:`StatFys.helpertt`."""
        return self.helpertt_terms(temp, n).sum()

    def helpern(self, temp, n):
        """See :meth:`StatFys.helpern`."""
        return self.helpern_terms(temp, n).sum()

    def helperv(self, temp, n):
        """See :meth:`StatFys.helperv`."""
        return self.helperv_terms(temp, n).sum()

    def helper_terms(self, temp, n):
        """Returns an array with all the helper results for the distinct terms.

           This is just an array version of :meth:`StatFys.helper`.
        """
        raise NotImplementedError

    def helpert_terms(self, temp, n):
        """Returns an array with all the helpert results for the distinct terms.

           This is just an array version of :meth:`StatFys.helpert`.
        """
        raise NotImplementedError

    def helpertt_terms(self, temp, n):
        """Returns an array with all the helpertt results for the distinct terms.

           This is just an array version of :meth:`StatFys.helpertt`.
        """
        raise NotImplementedError

    def helpern_terms(self, temp, n):
        """Returns an array with all the helpern results for the distinct terms.

           This is just an array version of :meth:`StatFys.helpern`.
        """
        # by default, this is the same as helper_terms
        return self.helper_terms(temp, n)

    def helperv_terms(self, temp, n):
        """Returns an array with all the helperv results for the distinct terms.

           This is just an array version of :meth:`StatFys.helperv`.
        """
        # by default, this is the same as helper_terms
        return self.helper_terms(temp, n)

    def log_terms(self, temp):
        """Returns an array with log results for the distinct terms.

           This is just an array version of :meth:`StatFys.log`.
        """
        return self.log(temp, self.helper_terms)

    def logt_terms(self, temp):
        """Returns an array with logt results for the distinct terms.

           This is just an array version of :meth:`StatFys.logt`.
        """
        return self.logt(temp, self.helpert_terms)

    def logtt_terms(self, temp):
        """Returns an array with logtt results for the distinct terms.

           This is just an array version of :meth:`StatFys.logtt`.
        """
        return self.logtt(temp, self.helpertt_terms)

    def logv_terms(self, temp):
        """Returns an array with logv results for the distinct terms.

           This is just an array version of :meth:`StatFys.logv`.
        """
        return self.logv(temp, self.helperv_terms)

    def logn_terms(self, temp):
        """Returns an array with logn results for the distinct terms.

           This is just an array version of :meth:`StatFys.logn`.
        """
        return self.logn(temp, self.helperv_terms)

    def internal_heat_terms(self, temp):
        """Returns an array with internal_heat results for the distinct terms.

           This is just an array version of :meth:`StatFys.internal_heat`.
        """
        return self.internal_heat(temp, self.helpert_terms)

    def heat_capacity_terms(self, temp):
        """Returns an array with heat_capacity results for the distinct terms.

           This is just an array version of :meth:`StatFys.heat_capacity`.
        """
        return self.heat_capacity(temp, self.helpert_terms, self.helpertt_terms)

    def entropy_terms(self, temp):
        """Returns an array with entropy results for the distinct terms.

           This is just an array version of :meth:`StatFys.entropy`.
        """
        return self.entropy(temp, self.helper_terms, self.helpert_terms)

    def free_energy_terms(self, temp):
        """Returns an array with free_energy results for the distinct terms.

           This is just an array version of :meth:`StatFys.free_energy`.
        """
        return self.free_energy(temp, self.helper_terms)

    def zero_point_energy_terms(self):
        """Returns an array with zero_point_energy results for the distinct terms.

           This is just an array version of :meth:`StatFys.chemical_potential`.
        """
        return self.zero_point_energy(self.helpern_terms)

    def chemical_potential_terms(self, temp):
        """Returns an array with chemical_potential results for the distinct terms.

           This is just an array version of :meth:`StatFys.chemical_potential`.
        """
        return self.chemical_potential(temp, self.helpern_terms)


def _check_levels(temp, bfs, Z):
    if bfs[bfs.argmin()]/Z > 0.01:
        raise ValueError('The highest energy level is occupied by more than 1%.')


def helper_levels(temp, n, energy_levels, check=False):
    """Helper 0 function for a system with the given energy levels.

       Returns T^n ln(Z), where Z is the partition function

       Arguments:
        | ``temp`` -- the temperature
        | ``n`` -- the power for the temperature factor
        | ``energy_levels`` -- an array with energy levels

       Optional argument:
        | ``check`` -- when set to True, an error is raise when the highest
                       energy level is occupied by more than 1% the the current
                       temperature.
    """
    # this is defined as a function because multiple classes need it
    if temp == 0:
        energy = energy_levels[0]
        degeneracy = 1
        while (degeneracy<len(energy_levels) and energy_levels[0]==energy_levels[degeneracy]):
            degeneracy += 1
        return float(temp)**n*np.log(degeneracy) - float(temp)**(n-1)*energy
    else:
        bfs = np.exp(-energy_levels/(boltzmann*temp))
        Z = bfs.sum()
        if check:
            _check_levels(temp, bfs, Z)
        return float(temp)**n*np.log(Z)

def helpert_levels(temp, n, energy_levels, check=False):
    """Helper 1 function for a system with the given energy levels.

       Returns T^n (d ln(Z) / dT), where Z is the partition function

       Arguments:
        | ``temp`` -- the temperature
        | ``n`` -- the power for the temperature factor
        | ``energy_levels`` -- an array with energy levels

       Optional argument:
        | ``check`` -- when set to True, an error is raise when the highest
                       energy level is occupied by more than 1% the the current
                       temperature.
    """
    # this is defined as a function because multiple classes need it
    if temp == 0:
        raise NotImplementedError
    else:
        bfs = np.exp(-energy_levels/(boltzmann*temp))
        Z = bfs.sum()
        if check:
            _check_levels(temp, bfs, Z)
        return float(temp)**(n-2)*(bfs*energy_levels).sum()/Z/boltzmann

def helpertt_levels(temp, n, energy_levels, check=False):
    """Helper 2 function for a system with the given energy levels.

       Returns T^n (d^2 ln(Z) / dT^2), where Z is the partition function

       Arguments:
        | ``temp`` -- the temperature
        | ``n`` -- the power for the temperature factor
        | ``energy_levels`` -- an array with energy levels

       Optional argument:
        | ``check`` -- when set to True, an error is raise when the highest
                       energy level is occupied by more than 1% the the current
                       temperature.
    """
    # this is defined as a function because multiple classes need it
    if temp == 0:
        raise NotImplementedError
    else:
        bfs = np.exp(-energy_levels/(boltzmann*temp))
        Z = bfs.sum()
        if check:
            _check_levels(temp, bfs, Z)
        return float(temp)**(n-4)/boltzmann**2*((bfs*energy_levels**2).sum()/Z) \
               -2*float(temp)**(n-3)/boltzmann*((bfs*energy_levels).sum()/Z) \
               -float(temp)**(n-4)/boltzmann**2*((bfs*energy_levels).sum()/Z)**2


class Electronic(Info, StatFys):
    """The electronic contribution to the partition function."""
    # TODO: this should also include the potential energy from the ab initio
    # computation. This is now added in the PartFun object.
    def __init__(self, multiplicity=None):
        """
           Optional argument:
            | ``multiplicity`` -- the spin multiplicity of the electronic system

           When the optional argument is not given, it is determined when the
           PartFun object is constructed.
        """
        self.multiplicity = multiplicity
        self.energy = None
        Info.__init__(self, "electronic")

    def init_part_fun(self, nma, partf):
        """See :meth:`StatFys.init_part_fun`."""
        if self.multiplicity is None:
            self.multiplicity = nma.multiplicity
            if self.multiplicity is None:
                raise ValueError("Spin multiplicity is not defined.")
        self.energy = nma.energy

    def dump(self, f):
        """See :meth:`Info.dump`."""
        Info.dump(self, f)
        print >> f, "    Multiplicity: %i" % self.multiplicity
        print >> f, "    Electronic energy: %.7f" % self.energy

    def helper(self, temp, n):
        """See :meth:`StatFys.helper`."""
        result = float(temp)**n*np.log(self.multiplicity)
        if temp == 0.0:
            if n < 1:
                raise NotImplementedError
            else:
                result -= self.energy/boltzmann
        else:
            result -= float(temp)**(n-1)*self.energy/boltzmann
        return result

    def helpert(self, temp, n):
        """See :meth:`StatFys.helpert`."""
        if temp == 0.0:
            if n < 2:
                raise NotImplementedError
            else:
                return self.energy/boltzmann
        else:
            return float(temp)**(n-2)*self.energy/boltzmann

    def helpertt(self, temp, n):
        """See :meth:`StatFys.helpertt`."""
        if temp == 0.0:
            if n < 3:
                raise NotImplementedError
            else:
                return -2.0*self.energy/boltzmann
        else:
            return -2.0*float(temp)**(n-3)*self.energy/boltzmann


class ExtTrans(Info, StatFys):
    """The contribution from the external translation.

       This contribution includes many body terms and optional constant pressure
       corrections. It is based on the classical ideal gas approximation.
    """

    def __init__(self, cp=True, pressure=None, density=None, dim=3, mobile=None):
        """
           Optional arguments:
            | ``cp`` -- When True, an additional factor is included in the
                        partition function to model a constant pressure (or
                        constant surface tension) ensemble instead of a constant
                        volume (or constant surface) ensemble.
            | ``pressure`` -- (only allowed when cp==True)
                              The external pressure exerted on the system in
                              case of the NpT ensemble. The default is 1 atm
                              for 3D gases. The default for 2D systems is 75.64
                              miliNewton per meter, i.e. the surface tension of
                              water. For other dimensions, the default is 1.0.
            | ``density`` -- (only allowed when cp==False)
                             The density of the system in case of the NVT
                             ensemble. The default is 1.0 mol/meter**dim.
            | ``dim`` -- The dimension of the ideal gas.
            | ``mobile`` -- A list of atom indexes that are free to translate. In
                            case of a mobile molecule adsorbed on a surface, only
                            include atom indexes of the adsorbate. The default is
                            that all atoms are mobile.
        """
        if dim <= 0:
            raise ValueError("The dimension of the gas must be strictly positive.")
        if dim == 3:
            self.pressure_unit = bar
            self.pressure_unit_name = "bar"
        elif dim == 2:
            self.pressure_unit = 1e-3*newton/meter
            self.pressure_unit_name = "mN/m"
        else:
            self.pressure_unit = 1.0
            self.pressure_unit_name = "a.u."
        self.density_unit = mol/meter**dim
        self.density_unit_name = "mol*meter**%i" % (-dim)

        if cp:
            if density is not None:
                raise ValueError("The density can not be fixed in the NpT ensemble, i.e. it depends on the temperature.")
            if pressure is None:
                if dim == 3:
                    pressure = 1*atm
                elif dim == 2:
                    # approximately the surface tension of water:
                    pressure = self.pressure_unit*75.64
                else:
                    # whatever...
                    pressure = 1.0
            self._pressure = pressure
        else:
            if pressure is not None:
                raise ValueError("The pressure can not be fixed in the NVT ensemble, i.e. it depends on the temperature.")
            if density is None:
                density = self.density_unit
            self._density = density

        self._cp = cp
        self.dim = dim
        self._mobile = mobile
        Info.__init__(self, "translational")

    cp = property(lambda self: self._cp)
    mobile = property(lambda self: self._mobile)

    def _get_pressure(self):
        """The pressure in case of an NpT ensemble."""
        if self.cp:
            return self._pressure
        else:
            raise ValueError("The pressure is not a known constant in the NVT ensemble, i.e. it depends on the temperature.")

    def _set_pressure(self, pressure):
        """Set the pressure in case of an NpT ensemble.

           This will raise a ValueError in case this is a constant volume ensemble.
        """
        if self.cp:
            self._pressure = pressure
        else:
            raise ValueError("The pressure can not be fixed in the NVT ensemble, i.e. it depends on the temperature.")

    pressure = property(_get_pressure, _set_pressure)

    def _get_density(self):
        """The density in case of an NVT ensemble."""
        if not self.cp:
            return self._density
        else:
            raise ValueError("The density is not a known constant in the NpT ensemble, i.e. it depends on the temperature.")

    def _set_density(self, density):
        """Set the density in case of an NVT ensemble.

           This will raise a ValueError in case this is a constant pressure ensemble.
        """
        if not self.cp:
            self._density = density
        else:
            raise ValueError("The density can not be fixed in the NpT ensemble, i.e. it depends on the temperature.")

    density = property(_get_density, _set_density)

    def init_part_fun(self, nma, partf):
        """See :meth:`StatFys.init_part_fun`."""
        if self.mobile is None:
            self.mass = nma.mass
        else:
            self.mass = nma.masses[self.mobile].sum()

    def dump(self, f):
        """See :meth:`Info.dump`."""
        Info.dump(self, f)
        print >> f, "    Dimension: %i" % self.dim
        print >> f, "    Constant pressure: %s" % self.cp
        if self.cp:
            print >> f, "    Pressure [%s]: %.5f" % (self.pressure_unit_name, self._pressure/self.pressure_unit)
        else:
            print >> f, "    Density [%s]: %.5f" % (self.density_unit_name, self._density/self.density_unit)
        if self.cp:
            print >> f, "      BIG FAT WARNING!!!"
            print >> f, "      This is an NpT partition function."
            print >> f, "      Internal heat contains a PV term (and is therefore the enthalpy)."
            print >> f, "      Free energy contains a PV term (and is therefore the Gibbs free energy)."
            print >> f, "      The heat capacity is computed at constant pressure."
        else:
            print >> f, "      BIG FAT WARNING!!!"
            print >> f, "      This is an NVT partition function."
            print >> f, "      Internal heat does NOT contain a PV term."
            print >> f, "      Free energy does NOT contain a PV term (and is therefore the Helmholtz free energy)."
            print >> f, "      The heat capacity is computed at constant volume."
        print >> f, "    Mass [amu]: %f" % (self.mass/amu)

    def _z1(self, temp):
        return 0.5*self.dim*np.log(2*np.pi*self.mass*boltzmann*temp/planck**2)

    def helper(self, temp, n):
        """See :meth:`StatFys.helper`."""
        if temp == 0:
            if n > 0:
                return 0.0
            else:
                raise NotImplementedError
        else:
            result = self._z1(temp)
            if self.cp:
                result += np.log(boltzmann*temp/self._pressure)
            else:
                result += 1.0 - np.log(self.density)
            return result*float(temp)**n

    def helpert(self, temp, n):
        """See :meth:`StatFys.helpert`."""
        if temp == 0:
            raise NotImplementedError
        else:
            result = 0.5*self.dim
            if self.cp:
                result += 1
            return result*float(temp)**(n-1)

    def helpertt(self, temp, n):
        """See :meth:`StatFys.helpertt`."""
        if temp == 0:
            raise NotImplementedError
        else:
            result = -0.5*self.dim
            if self.cp:
                result -= 1
            return result*float(temp)**(n-2)

    def helpern(self, temp, n):
        """See :meth:`StatFys.helpern`."""
        if temp == 0:
            if n > 0:
                return 0.0
            else:
                raise NotImplementedError
        else:
            result = self._z1(temp)
            if self.cp:
                result += np.log(boltzmann*temp/self._pressure)
            else:
                result += -np.log(self._density)
            return result*float(temp)**n

    def helperv(self, temp, n):
        r"""See :meth:`StatFys.helperv`."""
        if temp == 0:
            if n > 0:
                return 0.0
            else:
                raise NotImplementedError
        else:
            return self._z1(temp)*float(temp)**n


class ExtRot(Info, StatFys):
    """The contribution from the external rotation.

       This class is based on the integral approximation of the partition
       function.
    """
    def __init__(self, symmetry_number=None, im_threshold=1.0):
        """
           Optional arguments:
            | ``symmetry_number`` -- The rotational symmetry number of the
                                     molecule.
            | ``im_threshold``  --  When a moment of inertia drops below this
                                    threshold, it is discarded, which matters
                                    for linear molecules.
        """
        self.symmetry_number = symmetry_number
        self.im_threshold = im_threshold
        Info.__init__(self, "rotational")

    def init_part_fun(self, nma, partf):
        """See :meth:`StatFys.init_part_fun`."""
        if nma.periodic:
            raise ValueError("There is no external rotation in periodic systems.")
        self.inertia_tensor = nma.inertia_tensor
        self.moments = np.linalg.eigvalsh(nma.inertia_tensor)
        if self.symmetry_number == None:
            self.symmetry_number = nma.symmetry_number
            natom = len(nma.numbers)
            if self.symmetry_number == None and natom < 10:
                from molmod import Molecule
                # compute the rotational symmetry number
                tmp_mol = Molecule(nma.numbers, nma.coordinates)
                self.symmetry_number = tmp_mol.compute_rotsym()
            else:
                self.symmetry_number = 1
                print 'WARNING: molecule is too large (%i atoms > 10) to quickly estimate the rotational symmetry number.' % natom
        self.factor = np.sqrt(np.product([
            2*np.pi*m*boltzmann for m in self.moments if m > self.im_threshold
        ]))/self.symmetry_number/np.pi
        self.count = (self.moments > self.im_threshold).sum()

    def dump(self, f):
        """See :meth:`Info.dump`."""
        Info.dump(self, f)
        print >> f, "    Rotational symmetry number: %i" % self.symmetry_number
        print >> f, "    Moments of inertia [amu*bohr**2]: %f  %f %f" % tuple(self.moments/amu)
        print >> f, "    Threshold for non-zero moments of inertia [amu*bohr**2]: %e" % (self.im_threshold/amu)
        print >> f, "    Non-zero moments of inertia: %i" % self.count

    def helper(self, temp, n):
        """See :meth:`StatFys.helper`."""
        if temp == 0:
            if n > 0:
                return 0.0
            else:
                raise NotImplementedError
        else:
            return float(temp)**n*(np.log(temp)*0.5*self.count + np.log(self.factor))

    def helpert(self, temp, n):
        """See :meth:`StatFys.helpert`."""
        return float(temp)**(n-1)*0.5*self.count

    def helpertt(self, temp, n):
        """See :meth:`StatFys.helpertt`."""
        return -float(temp)**(n-2)*0.5*self.count


class PCMCorrection(Info, StatFys):
    """A correction to the free energy as function of the temperature.

       The correction can be a constant shift of the free energy or a linear
       shift of the free energy as function of the temperature.
    """

    def __init__(self, point1, point2=None):
        """
           Argument:
            | ``point1`` -- A 2-tuple with free energy and a temperature. A
                            correction for the free energy at the given
                            temperature. (If no second point is given, the same
                            correction is applied to all temperatures.)

           Optional argument:
            | ``point2`` -- A 2-tuple with free energy and a temperature. In
                            combination with point1, a linear free energy
                            correction as function of the temperature is added.
        """
        if (not hasattr(point1, "__len__")) or len(point1) != 2:
            raise ValueError("The first argument must be a (delta_G, temp) pair.")
        if point2 is not None and ((not hasattr(point2, "__len__")) or len(point2)) != 2:
            raise ValueError("The second argument must be None or a (delta_G, temp) pair.")
        self.point1 = point1
        self.point2 = point2
        Info.__init__(self, "pcm_correction")

    def dump(self, f):
        """See :meth:`Info.dump`."""
        Info.dump(self, f)
        print >> f, "    Point 1:"
        print >> f, "       Delta G [kJ/mol]: %.2f" % (self.point1[0]/kjmol)
        print >> f, "       Temperature [K]: %.2f" % (self.point1[1])
        print >> f, "    Point 2:"
        if self.point2 is not None:
            print >> f, "       Delta G [kJ/mol]: %.2f" % (self.point2[0]/kjmol)
            print >> f, "       Temperature [K]: %.2f" % (self.point2[1])
        else:
            print >> f, "       Not Defined!! Only rely on computations on temperature of point 1!!"
        print >> f, "    Zero-point contribution [kJ/mol]: %.7f" % (self.zero_point_energy()/kjmol)

    def _eval_free(self, temp):
        if self.point2 is None:
            return self.point1[0], 0.0, 0.0
        else:
            slope = (self.point2[0]-self.point1[0])/(self.point2[1]-self.point1[1])
            return (
                self.point1[0] + slope*(temp-self.point1[1]),
                slope,
                0.0
            )

    def helper(self, temp, n):
        """See :meth:`StatFys.helper`."""
        F, Fp, Fpp = self._eval_free(temp)
        return -F*float(temp)**(n-1)/boltzmann

    def helpert(self, temp, n):
        """See :meth:`StatFys.helpert`."""
        F, Fp, Fpp = self._eval_free(temp)
        return (F*float(temp)**(n-2) - Fp*float(temp)**(n-1))/boltzmann

    def helpertt(self, temp, n):
        """See :meth:`StatFys.helpertt`."""
        F, Fp, Fpp = self._eval_free(temp)
        return (-Fpp*float(temp)**(n-1) + 2*(Fp*float(temp)**(n-2) - F*float(temp)**(n-3)))/boltzmann


def helper_vibrations(temp, n, freqs, classical=False, freq_scaling=1, zp_scaling=1):
    """Helper 0 function for a set of harmonic oscillators.

       Returns T^n ln(Z), where Z is the partition function.

       Arguments:
        | ``temp`` -- the temperature
        | ``n`` -- the power for the temperature factor
        | ``freqs`` -- an array with frequencies

       Optional arguments:
        | ``classical`` -- When True, the classical partition function is used.
                           [default=False]
        | ``freq_scaling`` -- Scale the frequencies with the given factor.
                              [default=1]
        | ``freq_zp`` -- Scale the zero-point energy correction with the given.
                         factor [default=1]
    """
    # this is defined as a function because multiple classes need it
    if classical:
        if temp == 0:
            if n >= 1:
                return np.zeros(len(freqs))
            else:
                raise NotImplementedError
        else:
            Af = planck*freqs*freq_scaling/(boltzmann*temp)
            return -float(temp)**n*np.log(Af)
    else:
        # The zero point correction is included in the vibrational partition
        # function.
        if temp == 0:
            Abis = freqs*(0.5*planck*zp_scaling/boltzmann)
            if n >= 1:
                return -Abis*float(temp)**(n-1)
            else:
                raise NotImplementedError
        else:
            A = freqs*(planck/(boltzmann*temp))
            B = np.exp(-freq_scaling*A)
            return -((0.5*zp_scaling)*A + np.log(1 - B))*float(temp)**n

def helpert_vibrations(temp, n, freqs, classical=False, freq_scaling=1, zp_scaling=1):
    """Helper 1 function for a set of harmonic oscillators.

       Returns T^n (d ln(Z) / dT), where Z is the partition function.

       Arguments:
        | ``temp`` -- the temperature
        | ``n`` -- the power for the temperature factor
        | ``energy_levels`` -- an array with energy levels

       Optional arguments:
        | ``classical`` -- When True, the classical partition function is used
                           [default=False]
        | ``freq_scaling`` -- Scale the frequencies with the given factor
                              [default=1]
        | ``freq_zp`` -- Scale the zero-point energy correction with the given
                         factor [default=1]
    """
    # this is defined as a function because multiple classes need it
    if classical:
        if temp == 0:
            raise NotImplementedError
        else:
            result = float(temp)**(n-1)
            if hasattr(freqs, "__len__"):
                result *= np.ones(len(freqs))
            return result
    else:
        if temp == 0:
            raise NotImplementedError
        else:
            A = freqs*(planck/(boltzmann*temp))
            B = np.exp(-freq_scaling*A)
            C = B/(1 - B)
            return A*float(temp)**(n-1)*(0.5*zp_scaling + freq_scaling*C)

def helpertt_vibrations(temp, n, freqs, classical=False, freq_scaling=1, zp_scaling=1):
    """Helper 2 function for a set of harmonic oscillators.

       Returns T^n (d^2 ln(Z) / dT^2), where Z is the partition function.

       Arguments:
        | ``temp`` -- the temperature
        | ``n`` -- the power for the temperature factor
        | ``energy_levels`` -- an array with energy levels

       Optional arguments:
        | ``classical`` -- When True, the classical partition function is used
                           [default=False]
        | ``freq_scaling`` -- Scale the frequencies with the given factor
                              [default=1]
        | ``freq_zp`` -- Scale the zero-point energy correction with the given
                         factor [default=1]
    """
    # this is defined as a function because multiple classes need it
    if classical:
        if temp == 0:
            raise NotImplementedError
        else:
            result = -float(temp)**(n-2)
            if hasattr(freqs, "__len__"):
                result *= np.ones(len(freqs))
            return result
    else:
        if temp == 0:
            raise NotImplementedError
        else:
            A = freqs*(planck/(boltzmann*temp))
            Af = freq_scaling*A
            B = np.exp(-Af)
            C = B/(1.0 - B)
            return -A*float(temp)**(n-2)*(zp_scaling + freq_scaling*C*(2 - Af/(1-B)))


class Vibrations(Info, StatFysTerms):
    """The vibrational contribution to the partition function."""
    def __init__(self, classical=False, freq_scaling=1, zp_scaling=1, freq_threshold=None):
        """
           Optional arguments:
            | ``classical`` -- When True, the vibrations are treated classically
                               [default=False]
            | ``freq_scaling`` -- Scale factor for the frequencies [default=1]
            | ``zp_scaling`` -- Scale factor for the zero-point energy
                                correction [default=1]
            | ``freq_threshold`` -- Frequencies whose absolute value is below
                                    this threshold will be left out of the
                                    partition function, in addition to those
                                    already indicated as 'almost' zero by the
                                    NMA. This option is only needed to fix some
                                    pathological cases.
        """
        self.classical = classical
        self.freq_scaling = freq_scaling
        self.zp_scaling = zp_scaling
        self.freq_threshold = freq_threshold
        Info.__init__(self, "vibrational")

    def init_part_fun(self, nma, partf):
        """See :meth:`StatFys.init_part_fun`."""
        zero_indexes = nma.zeros
        nonzero_mask = np.ones(len(nma.freqs), dtype=bool)
        nonzero_mask[zero_indexes] = False
        if self.freq_threshold is not None:
            nonzero_mask[abs(nma.freqs) < self.freq_threshold] = False

        self.freqs = nma.freqs[nonzero_mask]
        self.zero_freqs = nma.freqs[~nonzero_mask]
        self.positive_freqs = nma.freqs[(nma.freqs > 0) & nonzero_mask]
        self.negative_freqs = nma.freqs[(nma.freqs < 0) & nonzero_mask]

        StatFysTerms.__init__(self, len(self.positive_freqs))

    def dump(self, f):
        """See :meth:`Info.dump`."""
        Info.dump(self, f)
        print >> f, "    Number of zero wavenumbers: %i " % (len(self.zero_freqs))
        print >> f, "    Number of real wavenumbers: %i " % (len(self.positive_freqs))
        print >> f, "    Number of imaginary wavenumbers: %i" % (len(self.negative_freqs))
        print >> f, "    Frequency scaling factor: %.4f" % self.freq_scaling
        print >> f, "    Zero-point scaling factor: %.4f" % self.zp_scaling
        self.dump_values(f, "Zero Wavenumbers [1/cm]", self.zero_freqs/(lightspeed/centimeter), "% 8.1f", 8)
        self.dump_values(f, "Real Wavenumbers [1/cm]", self.positive_freqs/(lightspeed/centimeter), "% 8.1f", 8)
        self.dump_values(f, "Imaginary Wavenumbers [1/cm]", self.negative_freqs/(lightspeed/centimeter), "% 8.1f", 8)
        print >> f, "    Zero-point contribution [kJ/mol]: %.7f" % (self.zero_point_energy()/kjmol)

    def helper_terms(self, temp, n):
        """See :meth:`StatFysTerms.helper_terms`."""
        return helper_vibrations(
            temp, n, self.positive_freqs, self.classical, self.freq_scaling,
            self.zp_scaling
        )

    def helpert_terms(self, temp, n):
        """See :meth:`StatFysTerms.helpert_terms`."""
        return helpert_vibrations(
            temp, n, self.positive_freqs, self.classical, self.freq_scaling,
            self.zp_scaling
        )

    def helpertt_terms(self, temp, n):
        """See :meth:`StatFysTerms.helpertt_terms`."""
        return helpertt_vibrations(
            temp, n, self.positive_freqs, self.classical, self.freq_scaling,
            self.zp_scaling
        )


class PartFun(Info, StatFys):
    """The partition function.

       This object contains all contributions to the partition function in
       self.terms and makes sure they are properly initialized. It also
       implements all the methods defined in StatFys, e.g. it can compute
       the entropy, the free energy and so on.
    """
    __reserved_names__ = set(["terms"])

    def __init__(self, nma, terms=None):
        """
           Arguments:
            | ``nma`` -- NMA object
           Optional arguments:
            | ``terms`` -- list to select the contributions to the partition
                           function e.g. [Vibrations(classical=True), ExtRot(1)]
        """
        if terms is None:
            terms = []
        self.terms = terms
        # perform a sanity check on the names of the contributions:
        for term in self.terms:
            if term.name in self.__reserved_names__:
                raise ValueError("An additional partition function term can not have the name '%s'" % term.name)
        # done testing, start initialization
        self.vibrational = None
        self.electronic = None
        for term in self.terms:
            self.__dict__[term.name] = term

        if self.vibrational is None:
            self.vibrational = Vibrations()
            self.terms.append(self.vibrational)

        if self.electronic is None:
            self.electronic = Electronic()
            self.terms.append(self.electronic)

        self.terms.sort(key=(lambda t: t.name))

        for term in self.terms:
            try:
                term.init_part_fun(nma, self)
            except Exception as err:
                if not err.args:
                   err.args=('',)
                err.args = ('%s: %s' % (term.name, err.message),)+err.args[1:]
                raise

        self.title = nma.title
        self.chemical_formula = nma.chemical_formula
        Info.__init__(self, "total")

    def helper(self, temp, n):
        """See :meth:`StatFys.helper`."""
        return sum(term.helper(temp, n) for term in self.terms)

    def helpert(self, temp, n):
        """See :meth:`StatFys.helpert`."""
        return sum(term.helpert(temp, n) for term in self.terms)

    def helpertt(self, temp, n):
        """See :meth:`StatFys.helpertt`."""
        return sum(term.helpertt(temp, n) for term in self.terms)

    def helpern(self, temp, n):
        """See :meth:`StatFys.helpern`."""
        return sum(term.helpern(temp, n) for term in self.terms)

    def helperv(self, temp, n):
        """See :meth:`StatFys.helperv`."""
        return sum(term.helperv(temp, n) for term in self.terms)

    def dump(self, f):
        """See :meth:`Info.dump`."""
        print >> f, "Title:", self.title
        print >> f, "Chemical formula:", self.chemical_formula
        print >> f, "Electronic energy [au]: %.5f" % self.electronic.energy
        print >> f, "Zero-point contribution [kJ/mol]: %.7f" % ((self.zero_point_energy() - self.electronic.energy)/kjmol)
        print >> f, "Zero-point energy [au]: %.5f" % self.zero_point_energy()
        print >> f, "Contributions to the partition function:"
        for term in self.terms:
            term.dump(f)

    def write_to_file(self, filename):
        """Write an extensive description of the parition function to a file.

           Argument:
            | ``filename`` -- The name of the file to write to.
        """
        f = file(filename, 'w')
        self.dump(f)
        f.close()
