# -*- coding: utf-8 -*-
# TAMkin is a post-processing toolkit for thermochemistry and kinetics analysis.
# Copyright (C) 2008-2010 Toon Verstraelen <Toon.Verstraelen@UGent.be>,
# Matthias Vandichel <Matthias.Vandichel@UGent.be> and
# An Ghysels <An.Ghysels@UGent.be>, Center for Molecular Modeling (CMM), Ghent
# University, Ghent, Belgium; all rights reserved unless otherwise stated.
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
# parts of this program are required to cite the following five articles:
#
# "TAMkin: A Versatile Package for Vibrational Analysis and Chemical Kinetics",
# An Ghysels, Toon Verstraelen, Karen Hemelsoet, Michel Waroquier and Veronique
# Van Speybroeck, Journal of Chemical Information and Modeling, Articles ASAP
# (As Soon As Publishable)
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
# --
"""Partition functions based on the harmonic oscillator approximation.

   The workhorse of this module is the **PartFun** class. A PartFun object
   represents a partition function with a solid interface. PartFun objects
   can be used to study chemical equilibrium, rate coefficients and various
   thermodynamic properties.

   These are the other classes and functions present in this module:

   * **Gas Laws:**
       * IdealGasLaw
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
   * **Auxiliary routines:**
       * compute_rate_coeff
       * compute_equilibrium_constant

   **Important**: Partition functions can be constructed for NpT gases, NVT
   gases and many other systems. The return values of methods such as
   ``free_energy``, ``internal_energy`` and ``heat_capacity`` are often given
   specialized names in the context of different partition functions. For
   example, chemists tend to use the following names:

   * 3D NVT gas:
       - ``PartFun.free_energy`` -> the Helmholtz free energy
       - ``PartFun.internal_energy`` -> the internal energy
       - ``PartFun.heat_capacity`` -> the heat capacity at constant volume

   * 3D NpT gas:
       - ``PartFun.free_energy`` -> the Gibbs free energy
       - ``PartFun.internal_energy`` -> the enthalpy
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
    planck

import numpy


__all__ = [
    "IdealGasLaw", "Info", "StatFys", "StatFysTerms",
    "helper_levels", "helpert_levels", "helpertt_levels",
    "Electronic", "ExtTrans", "ExtRot", "PCMCorrection",
    "Vibrations",
    "helper_vibrations", "helpert_vibrations", "helpertt_vibrations",
    "PartFun", "compute_rate_coeff", "compute_equilibrium_constant"
]



class IdealGasLaw(object):
    """Bundles several functions related to the ideal gas law."""

    def __init__(self, pressure=None, dim=3):
        """
           Optional argument:
             | ``pressure`` -- the external pressure of the system. The default
                               is 1 atm for 3D gases. The default for 2D systems
                               is 75.64 mili Newton per meter, i.e. the surface
                               tension of water. For other dimensions, the
                               default is 1.0.
             | ``dim`` -- The dimensionality of the gas.
        """
        if pressure is None:
            if dim == 3:
                pressure = 1*atm
            elif dim == 2:
                # approximately the surface tension of water in atomic units:
                pressure = 4.86e-05
            else:
                # whatever...
                pressure = 1.0
        self.pressure = pressure
        self.dim = dim
        # decide on the units
        if dim == 3:
            self.p_unit = bar
            self.p_unit_name = "bar"
        else:
            self.p_unit = 1.0
            self.p_unit_name = "a.u."

    def pv(self, temp, n):
        """PV function.

           Arguments:
            | ``temp`` -- The temperature
            | ``n`` -- A power for the additional temperature factor.

           This is an auxiliary function for the translational partition
           function. It returns the product of pressure and volume multiplied by
           a power of the temperature.
        """
        if temp == 0:
            if n > -1:
                return 0.0
            elif n == -1:
                return boltzmann
            else:
                raise NotImplementedError
        else:
            return boltzmann*temp**(n+1)

    def pvt(self, temp, n):
        """PV function T.

           Arguments:
            | ``temp`` -- The temperature
            | ``n`` -- A power for the additional temperature factor.

           This is an auxiliary function for the translational partition
           function. It returns the derivative towards the temperature product
           of pressure and volume, multiplied by a power of the temperature.
           (Derivation is performed prior to multiplication with T)
        """
        return 0.0

    def pvtt(self, temp, n):
        """PV function TT.

           Arguments:
            | ``temp`` -- The temperature
            | ``n`` -- A power for the additional temperature factor.

           This is an auxiliary function for the translational partition
           function. It returns the second derivative towards the temperature of
           the product of pressure and volume, multiplied by a power of the
           temperature. (Derivation is performed prior to multiplication with T)
        """
        return 0.0

    def helper(self, temp, n):
        r"""Helper function.

           Returns

           .. math:: T^n \ln(V(T))

           where :math:`V` is the volume per molecule.

           Arguments:
            | ``temp`` -- the temperature
            | ``n`` -- the power for the temperature factor
        """
        if temp == 0:
            if n > 0:
                return 0.0
            else:
                raise NotImplementedError
        else:
            return temp**n*numpy.log(boltzmann*temp/self.pressure)

    def helpert(self, temp, n):
        r"""Helper function T.

           Returns

           .. math:: T^n \frac{d \ln(V(T))}{d T}

           where :math:`V` is the volume per molecule and the derivative is
           taken at constant pressure.

           Arguments:
            | ``temp`` -- the temperature
            | ``n`` -- the power for the temperature factor
        """
        if temp == 0:
            raise NotImplementedError
        else:
            return temp**(n-1)

    def helpertt(self, temp, n):
        r"""Helper function TT.

           Returns

           .. math:: T^n \frac{d^2 \ln(V(T))}{d T^2}

           where :math:`V` is the volume per molecule and the derivative is
           taken at constant pressure.

           Arguments:
            | ``temp`` -- the temperature
            | ``n`` -- the power for the temperature factor
        """
        if temp == 0:
            raise NotImplementedError
        else:
            return -temp**(n-2)

    def helpern(self, temp, n):
        r"""Helper function N.

           Returns

           .. math:: T^n \left(\frac{N}{V}-\frac{P}{kT}\right)\frac{dV(T)}{dN}

           where :math:`V` is the volume per molecule, :math:`N` is the number
           of particles, :math:`P` is the pressure and the derivative is taken
           at constant pressure. This is part of the computation of the chemical
           potential at constant pressure and is non-zero for non-ideal gases.
        """
        return 0.0

    def _get_description(self):
        """A one-line summary of the gas law."""
        return "Ideal gas law, dimension = %i, pressure [%s] = %.5f" % (
            self.dim, self.p_unit_name, self.pressure/self.p_unit
        )

    description = property(_get_description)


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
           function and :math:`N` is the total number of particles.

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
           function and :math:`N` is the total number of particles.

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
           function and :math:`N` is the total number of particles.

           Arguments:
            | ``temp`` -- the temperature
            | ``n`` -- the power for the temperature factor
        """
        raise NotImplementedError

    def helpern(self, temp, n):
        r"""Helper function N.

           Returns

           .. math:: T^n N \frac{d \left(\frac{\ln(Z_N)}{N}\right)}{dN}

           where :math:`Z_N` is (the contribution to) the many body partition
           function and :math:`N` is the total number of particles. This is used to
           compute the chemical potential.

           Arguments:
            | ``temp`` -- the temperature
            | ``n`` -- the power for the temperature factor
        """
        # Return zero by default as this is the proper behavior for most
        # constributions to the partition function. The external translation is
        # (for now) the only exception.
        return 0.0

    def log(self, temp, helper=None):
        """The logarithm of the partition function.

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
        """The derivative towards temperature of the logarithm of the partition function.

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
        """The second derivative towards temperature of the logarithm of the partition function.

           Argument:
            | ``temp`` -- the temperature

           Optional arguments:
            | ``helpertt`` -- an alternative implementation of helpertt
                              [default=self.helpertt]
        """
        if helpertt is None:
            helpertt = self.helpertt
        return helpertt(temp, 0)

    def internal_energy(self, temp, helpert=None):
        """Computes the internal energy per molecule.

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

    def chemical_potential(self, temp, helper=None, helpern=None):
        """Computes the chemical potential.

           Argument:
            | ``temp`` -- the temperature

           Optional argument:
            | ``helper`` -- an alternative implementation of helper
                            [default=self.helper]
            | ``helper`` -- an alternative implementation of helpern
                            [default=self.helpern]

           Note: as opposed to most other methods, this is an intensive
           function!
        """
        if helper is None:
            helper = self.helper
        if helpern is None:
            helpern = self.helpern
        return -boltzmann*(helper(temp, 1) + helpern(temp, 1))


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

    def internal_energy_terms(self, temp):
        """Returns an array with internal_energy results for the distinct terms.

           This is just an array version of :meth:`StatFys.internal_energy`.
        """
        return self.internal_energy(temp, self.helpert_terms)

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

    def chemical_potential_terms(self, temp):
        """Returns an array with chemical_potential results for the distinct terms.

           This is just an array version of :meth:`StatFys.chemical_potential`.
        """
        return self.chemical_potential(temp, self.helper_terms)


def helper_levels(temp, n, energy_levels):
    """Helper 0 function for a system with the given energy levels.

       Returns T^n ln(Z), where Z is the partition function

       Arguments:
        | ``temp`` -- the temperature
        | ``n`` -- the power for the temperature factor
        | ``energy_levels`` -- an array with energy levels
    """
    # this is defined as a function because multiple classes need it
    if temp == 0:
        energy = energy_levels[0]
        degeneracy = 1
        while (degeneracy<len(energy_levels) and energy_levels[0]==energy_levels[degeneracy]):
            degeneracy += 1
        return temp**n*numpy.log(degeneracy) - temp**(n-1)*energy
    else:
        Z = numpy.exp(-energy_levels/(boltzmann*temp)).sum()
        return temp**n*numpy.log(Z)

def helpert_levels(temp, n, energy_levels):
    """Helper 1 function for a system with the given energy levels.

       Returns T^n (d ln(Z) / dT), where Z is the partition function

       Arguments:
        | ``temp`` -- the temperature
        | ``n`` -- the power for the temperature factor
        | ``energy_levels`` -- an array with energy levels
    """
    # this is defined as a function because multiple classes need it
    if temp == 0:
        raise NotImplementedError
    else:
        es = energy_levels
        bfs = numpy.exp(-es/(boltzmann*temp))
        Z = bfs.sum()
        return temp**(n-2)*(bfs*es).sum()/Z/boltzmann

def helpertt_levels(temp, n, energy_levels):
    """Helper 2 function for a system with the given energy levels.

       Returns T^n (d^2 ln(Z) / dT^2), where Z is the partition function

       Arguments:
        | ``temp`` -- the temperature
        | ``n`` -- the power for the temperature factor
        | ``energy_levels`` -- an array with energy levels
    """
    # this is defined as a function because multiple classes need it
    if temp == 0:
        raise NotImplementedError
    else:
        es = energy_levels
        bfs = numpy.exp(-es/(boltzmann*temp))
        Z = bfs.sum()
        return temp**(n-4)/boltzmann**2*((bfs*es**2).sum()/Z) \
               -2*temp**(n-3)/boltzmann*((bfs*es).sum()/Z) \
               -temp**(n-4)/boltzmann**2*((bfs*es).sum()/Z)**2


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
        Info.__init__(self, "electronic")

    def init_part_fun(self, nma, partf):
        """See :meth:`StatFys.init_part_fun`."""
        if self.multiplicity is None:
            self.multiplicity = nma.multiplicity
            if self.multiplicity is None:
                raise ValueError("Spin multiplicity is not defined.")

    def dump(self, f):
        """See :meth:`Info.dump`."""
        Info.dump(self, f)
        print >> f, "    Multiplicity: %i" % self.multiplicity

    def helper(self, temp, n):
        """See :meth:`StatFys.helper`."""
        return temp**n*numpy.log(self.multiplicity)

    def helpert(self, temp, n):
        """See :meth:`StatFys.helpert`."""
        return 0.0

    def helpertt(self, temp, n):
        """See :meth:`StatFys.helpertt`."""
        return 0.0


class ExtTrans(Info, StatFys):
    r"""The contribution from the external translation.

       In the translational contribution, we take into account the terms that
       are typical for the classical limit of the many body partition function.
       Strictly speaking, these additions are not due to the fact that there is
       translational freedom, so this is to some extent an ugly hack, but a very
       common and convenient one.

       ExtTrans has a second feature that goes beyond the limits of just a
       simple translational partition function. All corrections due to the NpT
       ensemble are included by default. Leaving out these corrections is
       optional (cp=False), which leads to a partition function for the NVT
       ensemble.

       The translational partition function of a single d-dimensional particle
       reads

       .. math:: Z_{1,\text{trans}} = \left(\frac{2\pi m k_B T}{h^2}\right)^{\frac{d}{2}}V,

       and the logarithm is

       .. math:: \ln(Z_{1,\text{trans}}) = \frac{d}{2}\ln\left(\frac{2\pi m k_B T}{h^2}\right) + \ln(V).

       ExtTrans models the logarithm of the many-body translational partition
       per particle, in the classical limit:

       .. math:: \frac{\ln(Z_{N,\text{trans}})}{N} = \frac{\ln(\frac{1}{N!}Z_{1,\text{trans}}^N)}{N},

       where N is the total number of particles. Using Stirlings approximation,
       this leads to:

       .. math:: \frac{\ln(Z_{N,\text{trans}})}{N} = \frac{-N \ln(N) + N}{N} + \ln(Z_{1,\text{trans}}).

       The first term is split into two terms, :math:`-\ln(N)` and :math:`1`.
       The former is pushed into the expression of the translational partition
       function, while the latter just remains where it is. The final
       expression is:

       .. math:: \frac{\ln(Z_{N,\text{trans}})}{N} = 1+\frac{d}{2}\ln\left(\frac{2\pi m k_B T}{h^2}\right) + \ln\left(\frac{V}{N}\right).

       From this derivation it is clear that the many-body effects and the
       translational part must be done together, because the separate
       contributions depend on the number of particles, which is annoying.

       Note that in the case of constant pressure, all derivatives towards
       temperature are taken at constant pressure, and there is an extra term

       .. math:: -\frac{PV}{k_BTN}

       which becomes -1 in the case of an ideal gas.
    """

    def __init__(self, cp=True, gaslaw=None, dim=3, mobile=None):
        """
           Optional arguments:
            | ``cp`` -- When True, an additional factor is included in the
                        partition function to model a constant pressure (or
                        constant surface tension) ensemble instead of a constant
                        volume (or constant surface) ensemble.
            | ``gaslaw`` -- The gas law that the system under study obeys. This
                            is used to evaluation the PV term, and also to
                            compute the derivative of the volume towards the
                            temperature under constant pressure. By default, the
                            ideal gas law is used.
            | ``dim`` -- The dimension of the ideal gas.
            | ``mobile`` -- A list of atom indexes that are free to translate. In
                            case of a mobile molecule adsorbed on a surface, only
                            include atom indexes of the adsorbate. The default is
                            that all atoms are mobile.

           The dimension of the gas determines the unit of the partition
           function as follows::

               unit = bohr**dim/particle

           In 3D it is a volume per particle, or an inverse concentration.
        """
        self.cp = cp
        if gaslaw is None:
            self.gaslaw = IdealGasLaw(dim=dim)
        else:
            self.gaslaw = gaslaw
            assert self.gaslaw.dim == dim
        self.dim = dim
        self.mobile = mobile
        Info.__init__(self, "translational")

    def init_part_fun(self, nma, partf):
        """See :meth:`StatFys.init_part_fun`."""
        if self.mobile is None:
            self.mass = nma.mass
        else:
            self.mass = nma.masses[self.mobile].sum()

    def set_pressure(self, pressure):
        """Update the pressure setting in the gaslaw.

           Argument:
            | ``pressure`` -- The new pressure.
        """
        self.gaslaw.pressure = pressure

    def dump(self, f):
        """See :meth:`Info.dump`."""
        Info.dump(self, f)
        print >> f, "    Gas law: %s" % self.gaslaw.description
        print >> f, "    Dimension: %i" % self.dim
        print >> f, "    Constant pressure: %s" % self.cp
        if self.cp:
            print >> f, "      BIG FAT WARNING!!!"
            print >> f, "      This is an NpT partition function."
            print >> f, "      Internal energy contains a PV term (and is therefore the enthalpy)."
            print >> f, "      Free energy contains a PV term (and is therefore the Gibbs free energy)."
            print >> f, "      The heat capacity is computed at constant pressure."
        else:
            print >> f, "      BIG FAT WARNING!!!"
            print >> f, "      This is an NVT partition function."
            print >> f, "      Internal energy does NOT contain a PV term."
            print >> f, "      Free energy does NOT contain a PV term (and is therefore the Helmholtz free energy)."
            print >> f, "      The heat capacity is computed at constant volume."
        print >> f, "    Mass [amu]: %f" % (self.mass/amu)

    def helper(self, temp, n):
        """See :meth:`StatFys.helper`."""
        if temp == 0:
            if n > 0:
                return 0.0
            else:
                raise NotImplementedError
        else:
            result = (
                temp**n + # This is due to the 1/N!
                temp**n*0.5*self.dim*numpy.log(2*numpy.pi*self.mass*boltzmann*temp/planck**2) +
                self.gaslaw.helper(temp, n) # this is the T^n*ln(V/N), the /N is due to 1/N!
            )
            if self.cp:
                result -= self.gaslaw.pv(temp, n-1)/boltzmann
            return result

    def helpert(self, temp, n):
        """See :meth:`StatFys.helpert`."""
        if temp == 0:
            raise NotImplementedError
        else:
            result = 0.5*self.dim*temp**(n-1)
            if self.cp:
                result += self.gaslaw.helpert(temp, n)
                result -= self.gaslaw.pvt(temp, n-1)/boltzmann
            return result

    def helpertt(self, temp, n):
        """See :meth:`StatFys.helpertt`."""
        if temp == 0:
            raise NotImplementedError
        else:
            result = -0.5*self.dim*temp**(n-2)
            if self.cp:
                result += self.gaslaw.helpertt(temp, n)
                result -= self.gaslaw.pvtt(temp, n-1)/boltzmann
            return result

    def helpern(self, temp, n):
        """See :meth:`StatFys.helpern`."""
        if temp==0:
            if n > 0:
                return 0.0
            else:
                raise NotImplementedError
        else:
            result = -temp**n
        if self.cp:
            result += self.gaslaw.helpern(temp, n)
        return result


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
        self.moments = numpy.linalg.eigvalsh(nma.inertia_tensor)
        if self.symmetry_number == None:
            self.symmetry_number = nma.symmetry_number
            if self.symmetry_number == None:
                from molmod import Molecule
                # compute the rotational symmetry number
                tmp_mol = Molecule(nma.numbers, nma.coordinates)
                self.symmetry_number = tmp_mol.compute_rotsym()
        self.factor = numpy.sqrt(numpy.product([
            2*numpy.pi*m*boltzmann for m in self.moments if m > self.im_threshold
        ]))/self.symmetry_number/numpy.pi
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
            return temp**n*(numpy.log(temp)*0.5*self.count + numpy.log(self.factor))

    def helpert(self, temp, n):
        """See :meth:`StatFys.helpert`."""
        return temp**(n-1)*0.5*self.count

    def helpertt(self, temp, n):
        """See :meth:`StatFys.helpertt`."""
        return -temp**(n-2)*0.5*self.count


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
        print >> f, "    Zero-point contribution [kJ/mol]: %.7f" % (self.free_energy(0.0)/kjmol)

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
        return -F*temp**(n-1)/boltzmann

    def helpert(self, temp, n):
        """See :meth:`StatFys.helpert`."""
        F, Fp, Fpp = self._eval_free(temp)
        return (F*temp**(n-2) - Fp*temp**(n-1))/boltzmann

    def helpertt(self, temp, n):
        """See :meth:`StatFys.helpertt`."""
        F, Fp, Fpp = self._eval_free(temp)
        return (-Fpp*temp**(n-1) + 2*(Fp*temp**(n-2) - F*temp**(n-3)))/boltzmann


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
            if n > 0:
                return numpy.zeros(len(freqs))
            else:
                raise NotImplementedError
        else:
            return temp**n*numpy.log(0.5*boltzmann*temp/numpy.pi/freqs*freq_scaling)
    else:
        # The zero point correction is included in the partition function and
        # should not be taken into account when computing the reaction barrier.
        pfb = numpy.pi*freqs/boltzmann
        if temp == 0:
            return -zp_scaling*pfb*temp**(n-1)
        else:
            return -zp_scaling*pfb*temp**(n-1) - numpy.log(1-numpy.exp(-2*freq_scaling*pfb/temp))*temp**n
        # This would be the version when the zero point energy corrections are
        # included in the energy difference when computing the reaction rate:
        #return -numpy.log(1-numpy.exp(exp_arg*freq_scaling))

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
            return temp**(n-1)*numpy.ones(len(freqs))
    else:
        if temp == 0:
            raise NotImplementedError
        else:
            pfb = numpy.pi*freqs/boltzmann
            return pfb*temp**(n-2)*(zp_scaling - 2*freq_scaling/(1 - numpy.exp(2*freq_scaling*pfb/temp)))

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
            return -temp**(n-2)*numpy.ones(len(freqs))
    else:
        if temp == 0:
            raise NotImplementedError
        else:
            pfb = numpy.pi*freqs/boltzmann
            return -2*pfb*temp**(n-3)*(zp_scaling - 2*freq_scaling/(1 - numpy.exp(2*freq_scaling*pfb/temp))) + \
                   +temp**(n-4)*(freq_scaling*pfb/numpy.sinh(freq_scaling*pfb/temp))**2


class Vibrations(Info, StatFysTerms):
    """The vibrational contribution to the partition function."""
    def __init__(self, classical=False, freq_scaling=1, zp_scaling=1):
        """
           Optional arguments:
            | ``classical`` -- When True, the vibrations are treated classically
                               [default=False]
            | ``freq_scaling`` -- Scale factor for the frequencies [default=1]
            | ``zp_scaling`` -- Scale factor for the zero-point energy
                                correction [default=1]
        """
        self.classical = classical
        self.freq_scaling = freq_scaling
        self.zp_scaling = zp_scaling
        Info.__init__(self, "vibrational")

    def init_part_fun(self, nma, partf):
        """See :meth:`StatFys.init_part_fun`."""
        zero_indexes = nma.zeros
        nonzero_mask = numpy.ones(len(nma.freqs), dtype=bool)
        nonzero_mask[zero_indexes] = False

        self.freqs = nma.freqs[nonzero_mask]
        self.zero_freqs = nma.freqs[zero_indexes]
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
        print >> f, "    Zero-point contribution [kJ/mol]: %.7f" % (self.free_energy(0.0)/kjmol)

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
                raise ValueError("An additional partition function term can not have the name '%s'" % mod.name)
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
            term.init_part_fun(nma, self)

        self.energy = nma.energy
        self.title = nma.title
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

    def internal_energy(self, temp):
        """Compute the internal energy.

           If self is a constant pressure ensemble of a regular 3D gas, the
           return value is the enthalpy. If self is a constant volume ensemble
           of a regular 3D gas, the return value is the internal energy.

           Arguments:
            | ``temp`` -- the temperature
        """
        # The molecular ground state energy is added here. It is tempting
        # to include it in the electronic part of partition function.
        return StatFys.internal_energy(self, temp) + self.energy

    def entropy(self, temp):
        """Compute the total entropy.

           Arguments:
            | ``temp`` -- the temperature
        """
        return StatFys.entropy(self, temp)

    def free_energy(self, temp):
        """Computes the free energy.

           If self is a constant pressure ensemble of a regular 3D gas, the
           return value is the Gibbs free energy. If self is a constant volume
           ensemble of a regular 3D gas, the return value is the Helmholtz free
           energy.

           Arguments:
            | ``temp`` -- the temperature
        """
        # The molecular ground state energy is added here. It is tempting
        # to include it in the electronic part of partition function.
        return StatFys.free_energy(self, temp) + self.energy

    def chemical_potential(self, temp):
        """Computes the chemical potential.

           Argument:
            | ``temp`` -- the temperature

           Note: as opposed to most other methods, this is an intensive
           function!
        """
        # The molecular ground state energy is added here. It is tempting
        # to include it in the electronic part of partition function.
        return StatFys.chemical_potential(self, temp) + self.energy

    def dump(self, f):
        """See :meth:`Info.dump`."""
        print >> f, "Title:", self.title
        print >> f, "Energy at T=0K [au]: %.5f" % self.energy
        print >> f, "Zero-point contribution [kJ/mol]: %.7f" % ((self.free_energy(0.0) - self.energy)/kjmol)
        print >> f, "Energy including zero-point contribution [au]: %.5f" % self.free_energy(0.0)
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


def compute_rate_coeff(pfs_react, pf_trans, temp, do_log=False):
    """Computes a (forward) rate coefficient.

       The implementation is based on transition state theory.

       Arguments:
         | ``pfs_react`` -- a list of partition functions objects, one for each
                            reactant
         | ``pf_trans`` -- the partition function of the transition state
         | ``temp`` -- the temperature

       Optional argument:
         | ``do_log`` -- Return the logarithm of the rate coefficient instead of
                         just the rate coefficient itself.
    """
    delta_A = pf_trans.free_energy(temp)
    delta_A -= sum(pf_react.free_energy(temp) for pf_react in pfs_react)
    log_result = -delta_A/(boltzmann*temp)
    for pf_react in pfs_react:
        if hasattr(pf_react, "translational"):
            log_result += pf_react.translational.gaslaw.helper(temp,0)
    if hasattr(pf_trans, "translational"):
        log_result -= pf_trans.translational.gaslaw.helper(temp,0)
    if do_log:
        return numpy.log(boltzmann*temp/planck) + log_result
    else:
        return boltzmann*temp/planck*numpy.exp(log_result)


def compute_equilibrium_constant(pfs_A, pfs_B, temp, do_log=False):
    """Computes the equilibrium constant between reactants and products.

       Arguments:
         | ``pfs_A`` -- a list of reactant partition functions
         | ``pfs_B`` -- a list of product partition functions
         | ``temp`` -- the temperature

       Optional argument:
         | ``do_log`` -- Return the logarithm of the equilibrium constant
                         instead of just the equilibrium constant itself.
    """
    delta_A = 0.0
    delta_A -= sum(pf_A.free_energy(temp) for pf_A in pfs_A)
    delta_A += sum(pf_B.free_energy(temp) for pf_B in pfs_B)
    log_K = -delta_A/(boltzmann*temp)

    for pf_A in pfs_A:
        if hasattr(pf_A, "translational"):
            log_K += pf_A.translational.gaslaw.helper(temp,0)
    for pf_B in pfs_B:
        if hasattr(pf_B, "translational"):
            log_K -= pf_B.translational.gaslaw.helper(temp,0)
    if do_log:
        return log_K
    else:
        return numpy.exp(log_K)
