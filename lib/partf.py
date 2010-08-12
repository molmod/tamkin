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
# "Vibrational Modes in partially optimized molecular systems.", An Ghysels,
# Dimitri Van Neck, Veronique Van Speybroeck, Toon Verstraelen and Michel
# Waroquier, Journal of Chemical Physics, Vol. 126 (22): Art. No. 224102, 2007
# DOI:10.1063/1.2737444
#
# "Cartesian formulation of the Mobile Block Hesian Approach to vibrational
# analysis in partially optimized systems", An Ghysels, Dimitri Van Neck and
# Michel Waroquier, Journal of Chemical Physics, Vol. 127 (16), Art. No. 164108,
# 2007
# DOI:10.1063/1.2789429
#
# "Calculating reaction rates with partial Hessians: validation of the MBH
# approach", An Ghysels, Veronique Van Speybroeck, Toon Verstraelen, Dimitri Van
# Neck and Michel Waroquier, Journal of Chemical Theory and Computation, Vol. 4
# (4), 614-625, 2008
# DOI:10.1021/ct7002836
#
# "Mobile Block Hessian approach with linked blocks: an efficient approach for
# the calculation of frequencies in macromolecules", An Ghysels, Veronique Van
# Speybroeck, Ewald Pauwels, Dimitri Van Neck, Bernard R. Brooks and Michel
# Waroquier, Journal of Chemical Theory and Computation, Vol. 5 (5), 1203-1215,
# 2009
# DOI:10.1021/ct800489r
#
# "Normal modes for large molecules with arbitrary link constraints in the
# mobile block Hessian approach", An Ghysels, Dimitri Van Neck, Bernard R.
# Brooks, Veronique Van Speybroeck and Michel Waroquier, Journal of Chemical
# Physics, Vol. 130 (18), Art. No. 084107, 2009
# DOI:10.1063/1.3071261
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
"""Partition functions based on the harmonic oscillator approximation & extensions

   The workhorse of this module is the **PartFun** class. A PartFun object
   represents a partition function with a interface that does not depend on the
   distinct contributions present in the partition function. PartFun objects
   can be used to study chemical equilibrium, rate coefficients and various
   thermodynamic properties.

   These are the other classes and functions present in this module:

   * **Gas Laws:** IdealGasLaw
   * **Abstract classes:** Info, StatFys, StatFysTerms
   * **Contributions:** Electronic, ExtTrans, ExtRot, Vibrations, Rotor
     (see rotor.py)
   * **Helper functions:** helper0_levels, helper1_levels, helper2_levels,
     helper0_vibrations, helper1_vibrations, helper2_vibrations
   * **Applications:** compute_rate_coeff, compute_equilibrium_constant (see
     pftools.py for a friendly interface to these functions)

   Note that all the extensive thermodynamic quantities computed here are in
   atomic units per molecule. If you want to express the Gibbs free energy at
   300 Kelvin of a system in kJ/mol, use the molmod module to perform unit
   conversions. For example::

     >>> pf = PartFun(...)
     >>> print pf.gibss_free_energy(300)/kjmol
"""


from molmod import boltzmann, lightspeed, atm, bar, amu, centimeter, kjmol

import numpy


__all__ = [
    "IdealGasLaw", "Info", "StatFys", "StatFysTerms",
    "helper0_levels", "helper1_levels", "helper2_levels",
    "Electronic", "ExtTrans", "ExtRot", "PCMCorrection",
    "Vibrations",
    "helper0_vibrations", "helper1_vibrations", "helper2_vibrations",
    "PartFun", "compute_rate_coeff", "compute_equilibrium_constant"
]



class IdealGasLaw(object):
    """Bundles several functions related to the ideal gas law."""

    def __init__(self, pressure=1*atm):
        """
           Optional argument:
             | pressure  --  the external pressure of the systen [default=1*atm]
        """
        self.pressure = pressure

    def pv(self, temp):
        """Returns the product of pressure and volume.

           Note that in SI units this is RT. The internal units of this module
           are atomic units and per molecule, so the return value becomes kT.
        """
        return boltzmann*temp

    def helper0(self, temp, n):
        """Helper function zero

           Returns T^n ln(V(T)), where V is the volume per molecule

           Arguments:
            | temp  --  the temperature
            | n  --  the power for the temperature factor
        """
        if temp == 0:
            if n > 0:
                return 0.0
            else:
                raise NotImplementedError
        else:
            return temp**n*numpy.log(boltzmann*temp/self.pressure)

    def helper1(self, temp, n):
        """Helper function one

           Returns T^n (d ln(V(T)))/(d T), where V is the volume per molecule

           Arguments:
            | temp  --  the temperature
            | n  --  the power for the temperature factor
        """
        if temp == 0:
            raise NotImplementedError
        else:
            return temp**(n-1)

    def helper2(self, temp, n):
        """Helper function two

           Returns T^n (d^2 ln(V(T)))/(d T^2), where V is the volume per molecule

           Arguments:
            | temp  --  the temperature
            | n  --  the power for the temperature factor
        """
        if temp == 0:
            raise NotImplementedError
        else:
            return -temp**(n-2)

    def _get_description(self):
        """A one-line summary of the gas law"""
        return "Ideal gas law, external pressure [bar] = %.5f" % (self.pressure/bar)

    description = property(_get_description)


class Info(object):
    """An object that has a name and that can dump info to a file."""
    def __init__(self, name):
        """
           Arguments:
            | name  --  the name used for this object in the output
        """
        self.name = name

    def dump(self, f):
        """Write a description to file

           Arguments:
            | f  --  the file object to write to
        """
        print >> f, "  %s" % self.name.upper()

    def dump_values(self, f, label, values, format, num_col=8):
        """Write a nicely formatted array of numbers to file

           Arguments:
            | f  --  the file object to write to
            | label  --  a label that explains the meaning of the numbers
            | values  --  the array with numbers
            | format  --  a Python format string for one number

           Optional argumet:
            | num_col  --  the number of columns [default=8]
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
    """Abstract class for (contributions to) the parition function

       The constructor (__init__) and four methods (init_part_fun, helper0,
       helper1, helper2) must be implemented in derived classes.
    """
    def init_part_fun(self, nma, partf):
        """Compute parameters that depend on nma and partition function

           Arguments:
            | nma  --  an NMA object
            | partf  --  A PartFun object

           This method is called by the PartFun object and should not be called
           by the user. When the PartFun object is initialized, each
           contribution can further initialize its parameters based on the
           information available in the nma and partf objects.
        """
        pass

    def helper0(self, temp, n):
        """Helper function zero

           Returns T^n ln(Z), where Z is the partition function

           Arguments:
            | temp  --  the temperature
            | n  --  the power for the temperature factor
        """
        raise NotImplementedError

    def helper1(self, temp, n, cp=False):
        """Helper function one

           Returns T^n (d ln(Z) / dT), where Z is the partition function

           Arguments:
            | temp  --  the temperature
            | n  --  the power for the temperature factor

           Optional argument:
            | cp  --  When True, the derivative is taken at constant pressure,
                      otherwise it is taken at constant volume. [default=False]
        """
        raise NotImplementedError

    def helper2(self, temp, n, cp=False):
        """Helper function two

           Returns T^n (d^2 ln(Z) / dT^2), where Z is the partition function

           Arguments:
            | temp  --  the temperature
            | n  --  the power for the temperature factor

           Optional argument:
            | cp  --  When True, the derivative is taken at constant pressure,
                      otherwise it is taken at constant volume. [default=False]
        """
        raise NotImplementedError

    def log(self, temp, helper0=None):
        """The logarithm of the partition function

           Argument:
            | temp  --  the temperature

           Optional argument:
            | helper0  --  an alternative implementation of helper0
                           [default=self.helper0]
        """
        if helper0 is None:
            helper0 = self.helper0
        return helper0(temp, 0)

    def dlog(self, temp, helper1=None, cp=False):
        """The derivative towards temperature of the logarithm of the partition function

           Argument:
            | temp  --  the temperature

           Optional arguments:
            | helper1  --  an alternative implementation of helper1
                           [default=self.helper1]
            | cp  --  When True, the derivative is taken at constant pressure,
                      otherwise it is taken at constant volume. [default=False]
        """
        if helper1 is None:
            helper1 = self.helper1
        return helper1(temp, 0, cp)

    def ddlog(self, temp, helper2=None, cp=False):
        """The second derivative towards temperature of the logarithm of the partition function

           Argument:
            | temp  --  the temperature

           Optional arguments:
            | helper2  --  an alternative implementation of helper2
                           [default=self.helper2]
            | cp  --  When True, the derivative is taken at constant pressure,
                      otherwise it is taken at constant volume. [default=False]
        """
        if helper2 is None:
            helper2 = self.helper2
        return helper2(temp, 0, cp)

    def internal_energy(self, temp, helper1=None):
        """Computes the internal energy per molecule

           Argument:
            | temp  --  the temperature

           Optional argument:
            | helper1  --  an alternative implementation of helper1
                           [default=self.helper1]
        """
        if helper1 is None:
            helper1 = self.helper1
        return boltzmann*helper1(temp, 2)

    def heat_capacity_v(self, temp, helper1=None, helper2=None):
        """Computes the heat capacity per molecule at constant pressure

           Argument:
            | temp  --  the temperature

           Optional arguments:
            | helper1  --  an alternative implementation of helper1
                           [default=self.helper1]
            | helper2  --  an alternative implementation of helper2
                           [default=self.helper2]
        """
        if helper1 is None:
            helper1 = self.helper1
        if helper2 is None:
            helper2 = self.helper2
        return boltzmann*(2*helper1(temp, 1) + helper2(temp, 2))

    def heat_capacity_p(self, temp, helper1=None, helper2=None):
        """Computes the heat capacity per molecule at constant pressure

           Argument:
            | temp  --  the temperature

           Optional arguments:
            | helper1  --  an alternative implementation of helper1
                           [default=self.helper1]
            | helper2  --  an alternative implementation of helper2
                           [default=self.helper2]
        """
        if helper1 is None:
            helper1 = self.helper1
        if helper2 is None:
            helper2 = self.helper2
        return boltzmann*(2*helper1(temp, 1, cp=True) + helper2(temp, 2, cp=True))

    def entropy(self, temp, helper0=None, helper1=None):
        """Computes the entropy contribution per molecule

           Argument:
            | temp  --  the temperature

           Optional arguments:
            | helper0  --  an alternative implementation of helper0
                           [default=self.helper0]
            | helper1  --  an alternative implementation of helper1
                           [default=self.helper1]
        """
        if helper0 is None:
            helper0 = self.helper0
        if helper1 is None:
            helper1 = self.helper1
        return boltzmann*(helper0(temp, 0) + helper1(temp, 1))

    def free_energy(self, temp, helper0=None):
        """Computes the free energy per molecule

           Argument:
            | temp  --  the temperature

           Optional argument:
            | helper0  --  an alternative implementation of helper0
                           [default=self.helper0]

           Note: at this point there is no real distinction between Helmoholtz
           and Gibbs free energy. This distinction is only introduced at the
           level of the PartFun object.
        """
        if helper0 is None:
            helper0 = self.helper0
        return -boltzmann*helper0(temp, 1)


class StatFysTerms(StatFys):
    """Abstract class for (contributions to) the parition function with multiple terms

       The different terms (or factors if you like) are of the same mathematical
       structure.

       The constructor (__init__) and the four methods (init_part_fun,
       helper0_terms, helper1_terms, helper2_terms) must be implemented in
       derived classes.
    """
    def __init__(self, num_terms):
        """
           Arguments:
             num_terms  --  the number of terms present in this contribution.
        """
        self.num_terms = num_terms

    def helper0(self, temp, n):
        """See :meth:`StatFys.helper0`"""
        return self.helper0_terms(temp, n).sum()

    def helper1(self, temp, n, cp=False):
        """See :meth:`StatFys.helper1`"""
        return self.helper1_terms(temp, n, cp).sum()

    def helper2(self, temp, n, cp=False):
        """See :meth:`StatFys.helper2`"""
        return self.helper2_terms(temp, n).sum()

    def helper0_terms(self, temp, n):
        """Returns an array with all the helper0 results for the distinct terms.

           This is just an array version of :meth:`StatFys.helper0`.
        """
        raise NotImplementedError

    def helper1_terms(self, temp, n, cp=False):
        """Returns an array with all the helper1 results for the distinct terms.

           This is just an array version of :meth:`StatFys.helper1`.
        """
        raise NotImplementedError

    def helper2_terms(self, temp, n, cp=False):
        """Returns an array with all the helper2 results for the distinct terms.

           This is just an array version of :meth:`StatFys.helper2`.
        """
        raise NotImplementedError

    def log_terms(self, temp):
        """Returns an array with log results for the distinct terms.

           This is just an array version of :meth:`StatFys.log`.
        """
        return self.log(temp, self.helper0_terms)

    def dlog_terms(self, temp, cp=False):
        """Returns an array with dlog results for the distinct terms.

           This is just an array version of :meth:`StatFys.dlog`.
        """
        return self.dlog(temp, self.helper1_terms, cp)

    def ddlog_terms(self, temp, cp=False):
        """Returns an array with ddlog results for the distinct terms.

           This is just an array version of :meth:`StatFys.ddlog`.
        """
        return self.ddlog(temp, self.helper2_terms, cp)

    def internal_energy_terms(self, temp):
        """Returns an array with internal_energy results for the distinct terms.

           This is just an array version of :meth:`StatFys.internal_energy`.
        """
        return self.internal_energy(temp, self.helper1_terms)

    def heat_capacity_v_terms(self, temp):
        """Returns an array with heat_capacity_v results for the distinct terms.

           This is just an array version of :meth:`StatFys.heat_capacity_v`.
        """
        return self.heat_capacity_v(temp, self.helper1_terms, self.helper2_terms)

    def heat_capacity_p_terms(self, temp):
        """Returns an array with heat_capacity_p results for the distinct terms.

           This is just an array version of :meth:`StatFys.heat_capacity_v`.
        """
        return self.heat_capacity_p(temp, self.helper1_terms, self.helper2_terms)

    def entropy_terms(self, temp):
        """Returns an array with entropy results for the distinct terms.

           This is just an array version of :meth:`StatFys.entropy`.
        """
        return self.entropy(temp, self.helper0_terms, self.helper1_terms)

    def free_energy_terms(self, temp):
        """Returns an array with free_energy results for the distinct terms.

           This is just an array version of :meth:`StatFys.free_energy`.
        """
        return self.free_energy(temp, self.helper0_terms)


def helper0_levels(temp, n, energy_levels):
    """Helper 0 function for a system with the given energy levels

       Returns T^n ln(Z), where Z is the partition function

       Arguments:
        | temp  --  the temperature
        | n  --  the power for the temperature factor
        | energy_levels  --  an array with energy levels
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

def helper1_levels(temp, n, energy_levels):
    """Helper 1 function for a system with the given energy levels

       Returns T^n (d ln(Z) / dT), where Z is the partition function

       Arguments:
        | temp  --  the temperature
        | n  --  the power for the temperature factor
        | energy_levels  --  an array with energy levels
    """
    # this is defined as a function because multiple classes need it
    if temp == 0:
        raise NotImplementedError
    else:
        es = energy_levels
        bfs = numpy.exp(-es/(boltzmann*temp))
        Z = bfs.sum()
        return temp**(n-2)*(bfs*es).sum()/Z/boltzmann

def helper2_levels(temp, n, energy_levels):
    """Helper 2 function for a system with the given energy levels

       Returns T^n (d^2 ln(Z) / dT^2), where Z is the partition function

       Arguments:
        | temp  --  the temperature
        | n  --  the power for the temperature factor
        | energy_levels  --  an array with energy levels
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
    """The electronic contribution to the partition function"""
    # TODO: this should also include the potential energy from the ab initio
    # computation. This is now added in the PartFun object.
    def __init__(self, multiplicity=None):
        """
           Optional argument:
            | multiplicity  --  the spin multiplicity of the electronic system

           When the optional argument is not given, it is determined when the
           PartFun object is constructed.
        """
        self.multiplicity = multiplicity
        Info.__init__(self, "electronic")

    def init_part_fun(self, nma, partf):
        """See :meth:`StatFys.init_part_fun`"""
        if self.multiplicity is None:
            self.multiplicity = nma.multiplicity
            if self.multiplicity is None:
                raise ValueError("Spin multiplicity is not defined.")

    def dump(self, f):
        """See :meth:`Info.dump`"""
        Info.dump(self, f)
        print >> f, "    Multiplicity: %i" % self.multiplicity

    def helper0(self, temp, n):
        """See :meth:`StatFys.helper0`"""
        return temp**n*numpy.log(self.multiplicity)

    def helper1(self, temp, n, cp=False):
        """See :meth:`StatFys.helper1`"""
        return 0.0

    def helper2(self, temp, n, cp=False):
        """See :meth:`StatFys.helper2`"""
        return 0.0


class ExtTrans(Info, StatFys):
    """The contribution from the external translation"""

    def __init__(self):
        self.gaslaw = None
        Info.__init__(self, "translational")

    def init_part_fun(self, nma, partf):
        """See :meth:`StatFys.init_part_fun`"""
        self.mass = nma.mass
        self.gaslaw = partf.gaslaw

    def dump(self, f):
        """See :meth:`Info.dump`"""
        Info.dump(self, f)
        print >> f, "    Mass [amu]: %f" % (self.mass/amu)

    def helper0(self, temp, n):
        """See :meth:`StatFys.helper0`"""
        if temp == 0:
            if n > 0:
                return 0.0
            else:
                raise NotImplementedError
        else:
            return (
                temp**n*1.5*numpy.log(0.5*self.mass*boltzmann*temp/numpy.pi) +
                self.gaslaw.helper0(temp, n)
            )

    def helper1(self, temp, n, cp=False):
        """See :meth:`StatFys.helper1`"""
        if temp == 0:
            raise NotImplementedError
        else:
            result = 1.5*temp**(n-1)
            if cp:
                result += self.gaslaw.helper1(temp, n)
            return result

    def helper2(self, temp, n, cp=False):
        """See :meth:`StatFys.helper2`"""
        if temp == 0:
            raise NotImplementedError
        else:
            result = -1.5*temp**(n-2)
            if cp:
                result += self.gaslaw.helper2(temp, n)
            return result



class ExtRot(Info, StatFys):
    """The contribution from the external rotation"""
    def __init__(self, symmetry_number=0, im_threshold=1.0):
        self.symmetry_number = symmetry_number
        self.im_threshold = im_threshold
        Info.__init__(self, "rotational")

    def init_part_fun(self, nma, partf):
        """See :meth:`StatFys.init_part_fun`"""
        if nma.periodic:
            raise ValueError("There is no external rotation in periodic systems.")
        self.inertia_tensor = nma.inertia_tensor
        self.moments = numpy.linalg.eigvalsh(nma.inertia_tensor)
        if self.symmetry_number == 0:
            self.symmetry_number = nma.symmetry_number
            if self.symmetry_number == 0:
                from molmod import Molecule
                # compute the rotational symmetry number
                tmp_mol = Molecule(nma.numbers, nma.coordinates)
                self.symmetry_number = tmp_mol.compute_rotsym()
        self.factor = numpy.sqrt(numpy.product([
            2*numpy.pi*m*boltzmann for m in self.moments if m > self.im_threshold
        ]))/self.symmetry_number/numpy.pi
        self.count = (self.moments > self.im_threshold).sum()

    def dump(self, f):
        """See :meth:`Info.dump`"""
        Info.dump(self, f)
        print >> f, "    Rotational symmetry number: %i" % self.symmetry_number
        print >> f, "    Moments of inertia [amu*bohr**2]: %f  %f %f" % tuple(self.moments/amu)
        print >> f, "    Threshold for non-zero moments of inertia [amu*bohr**2]: %e" % (self.im_threshold/amu)
        print >> f, "    Non-zero moments of inertia: %i" % self.count

    def helper0(self, temp, n):
        """See :meth:`StatFys.helper0`"""
        if temp == 0:
            if n > 0:
                return 0.0
            else:
                raise NotImplementedError
        else:
            return temp**n*(numpy.log(temp)*0.5*self.count + numpy.log(self.factor))

    def helper1(self, temp, n, cp=False):
        """See :meth:`StatFys.helper1`"""
        return temp**(n-1)*0.5*self.count

    def helper2(self, temp, n, cp=False):
        """See :meth:`StatFys.helper2`"""
        return -temp**(n-2)*0.5*self.count


class PCMCorrection(Info, StatFys):
    """A correction to the free energy as function of the temperature

       The correction can be a constant shift of the free energy or a linear
       shift of the free energy as function of the temperature.
    """

    def __init__(self, point1, point2=None):
        """
        Argument:
         | point1  --  A 2-tuple with free energy and a temperature. A
                       correction for the free energy at the given temperature.
                       (If no second point is given, the same correction is
                       applied to all temperatures.)

        Optional argument:
         | point2  --  A 2-tuple with free energy and a temperature. In
                       combination with point1, a linear free energy correction
                       as function of the temperature is added.
        """
        if (not hasattr(point1, "__len__")) or len(point1) != 2:
            raise ValueError("The first argument must be a (delta_G, temp) pair.")
        if point2 is not None and ((not hasattr(point2, "__len__")) or len(point2)) != 2:
            raise ValueError("The second argument must be None or a (delta_G, temp) pair.")
        self.point1 = point1
        self.point2 = point2
        Info.__init__(self, "pcm_correction")

    def dump(self, f):
        """See :meth:`Info.dump`"""
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
        print >> f, "    Free energy contribution at T=0K [au]: %.7f" % self.free_energy(0.0)

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

    def helper0(self, temp, n):
        """See :meth:`StatFys.helper0`"""
        F, Fp, Fpp = self._eval_free(temp)
        return -F*temp**(n-1)/boltzmann

    def helper1(self, temp, n, cp=False):
        """See :meth:`StatFys.helper1`"""
        F, Fp, Fpp = self._eval_free(temp)
        return (F*temp**(n-2) - Fp*temp**(n-1))/boltzmann

    def helper2(self, temp, n, cp=False):
        """See :meth:`StatFys.helper2`"""
        F, Fp, Fpp = self._eval_free(temp)
        return (-Fpp*temp**(n-1) + 2*(Fp*temp**(n-2) - F*temp**(n-3)))/boltzmann


def helper0_vibrations(temp, n, freqs, classical=False, freq_scaling=1, zp_scaling=1):
    """Helper 0 function for a set of harmonic oscillators

       Returns T^n ln(Z), where Z is the partition function

       Arguments:
        | temp  --  the temperature
        | n  --  the power for the temperature factor
        | freqs  --  an array with frequencies

       Optional arguments:
        | classical  --  When True, the classical partition function is used
                         [default=False]
        | freq_scaling  --  Scale the frequencies with the given factor
                            [default=1]
        | freq_zp  --  Scale the zero-point energy correction with the given
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

def helper1_vibrations(temp, n, freqs, classical=False, freq_scaling=1, zp_scaling=1):
    """Helper 1 function for a set of harmonic oscillators

       Returns T^n (d ln(Z) / dT), where Z is the partition function

       Arguments:
        | temp  --  the temperature
        | n  --  the power for the temperature factor
        | energy_levels  --  an array with energy levels

       Optional arguments:
        | classical  --  When True, the classical partition function is used
                         [default=False]
        | freq_scaling  --  Scale the frequencies with the given factor
                            [default=1]
        | freq_zp  --  Scale the zero-point energy correction with the given
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

def helper2_vibrations(temp, n, freqs, classical=False, freq_scaling=1, zp_scaling=1):
    """Helper 2 function for a set of harmonic oscillators

       Returns T^n (d^2 ln(Z) / dT^2), where Z is the partition function

       Arguments:
        | temp  --  the temperature
        | n  --  the power for the temperature factor
        | energy_levels  --  an array with energy levels

       Optional arguments:
        | classical  --  When True, the classical partition function is used
                         [default=False]
        | freq_scaling  --  Scale the frequencies with the given factor
                            [default=1]
        | freq_zp  --  Scale the zero-point energy correction with the given
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
    """The vibrational contribution to the partition function"""
    def __init__(self, classical=False, freq_scaling=1, zp_scaling=1):
        """
        Optional arguments:
         | classical  --  When True, the vibrations are treated classically
                          [default=False]
         | freq_scaling  --  Scale factor for the frequencies [default=1]
         | zp_scaling  --  Scale factor for the zero-point energy correction
                           [default=1]
        """
        self.classical = classical
        self.freq_scaling = freq_scaling
        self.zp_scaling = zp_scaling
        Info.__init__(self, "vibrational")

    def init_part_fun(self, nma, partf):
        """See :meth:`StatFys.init_part_fun`"""
        zero_indexes = nma.zeros
        nonzero_mask = numpy.ones(len(nma.freqs), dtype=bool)
        nonzero_mask[zero_indexes] = False

        self.freqs = nma.freqs[nonzero_mask]
        self.zero_freqs = nma.freqs[zero_indexes]
        self.positive_freqs = nma.freqs[(nma.freqs > 0) & nonzero_mask]
        self.negative_freqs = nma.freqs[(nma.freqs < 0) & nonzero_mask]

        StatFysTerms.__init__(self, len(self.positive_freqs))

    def dump(self, f):
        """See :meth:`Info.dump`"""
        Info.dump(self, f)
        print >> f, "    Number of zero wavenumbers: %i " % (len(self.zero_freqs))
        print >> f, "    Number of real wavenumbers: %i " % (len(self.positive_freqs))
        print >> f, "    Number of imaginary wavenumbers: %i" % (len(self.negative_freqs))
        print >> f, "    Frequency scaling factor: %.4f" % self.freq_scaling
        print >> f, "    Zero-Point scaling factor: %.4f" % self.zp_scaling
        self.dump_values(f, "Zero Wavenumbers [1/cm]", self.zero_freqs/(lightspeed/centimeter), "% 8.1f", 8)
        self.dump_values(f, "Real Wavenumbers [1/cm]", self.positive_freqs/(lightspeed/centimeter), "% 8.1f", 8)
        self.dump_values(f, "Imaginary Wavenumbers [1/cm]", self.negative_freqs/(lightspeed/centimeter), "% 8.1f", 8)
        print >> f, "    Free energy contribution at T=0 [kJ/mol]: %.7f" % (self.free_energy(0.0)/kjmol)

    def helper0_terms(self, temp, n):
        """See :meth:`StatFysTerms.helper0_terms`"""
        return helper0_vibrations(
            temp, n, self.positive_freqs, self.classical, self.freq_scaling,
            self.zp_scaling
        )

    def helper1_terms(self, temp, n, cp=False):
        """See :meth:`StatFysTerms.helper1_terms`"""
        return helper1_vibrations(
            temp, n, self.positive_freqs, self.classical, self.freq_scaling,
            self.zp_scaling
        )

    def helper2_terms(self, temp, n, cp=False):
        """See :meth:`StatFysTerms.helper2_terms`"""
        return helper2_vibrations(
            temp, n, self.positive_freqs, self.classical, self.freq_scaling,
            self.zp_scaling
        )


class PartFun(Info, StatFys):
    """The partition function

       This object contains all contributions to the partition function in
       self.terms and makes sure they are properly initialized. It also
       implements all the methods defined in StatFys, e.g. it can compute
       the entropy, the free energy and so on.
    """
    __reserved_names__ = set(["terms"])

    def __init__(self, nma, terms=None, gaslaw=None):
        """
        Arguments:
          | nma  --  NMA object
        Optional arguments:
          | terms  --  list to select the contributions to the partition function
                       e.g. [Vibrations(classical=True), ExtRot(1)]
          | gaslaw  --  the gas law that the system under study obeys. This is
                        used to evaluation the PV term in the enthalpy and the
                        Gibbs free energy, and also to compute the derivative of
                        the volume towards the temperature under constant
                        pressure (required for the heat capacity at constant
                        pressure). By default, the ideal gas law is used.
        """
        if terms is None:
            terms = []
        self.terms = terms
        if gaslaw is None:
            self.gaslaw = IdealGasLaw()
        else:
            self.gaslaw = gaslaw
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
        Info.__init__(self, "total")

    def helper0(self, temp, n):
        """See :meth:`StatFys.helper0`"""
        return sum(term.helper0(temp, n) for term in self.terms)

    def helper1(self, temp, n, cp=False):
        """See :meth:`StatFys.helper1`"""
        return sum(term.helper1(temp, n, cp) for term in self.terms)

    def helper2(self, temp, n, cp=False):
        """See :meth:`StatFys.helper2`"""
        return sum(term.helper2(temp, n, cp) for term in self.terms)

    def internal_energy(self, temp):
        """Compute the internal energy

           Arguments:
            | temp  --  the temperature
        """
        return StatFys.internal_energy(self, temp) + self.energy

    def enthalpy(self, temp):
        """Compute the enthalpy

           Arguments:
            | temp  --  the temperature
        """
        return self.internal_energy(temp) + self.gaslaw.pv(temp)

    def entropy(self, temp):
        """Compute the total entropy

           Arguments:
            | temp  --  the temperature
        """
        # This is a bit tricky: There is a contribution to the total entropy
        # that can not be associated with one of the factors in the partition
        # function. It is purely related to the N factorial in the denominator
        # of the partition function. This controbution is one boltzmann constant
        # per partical (or one universal gas constant per mol)
        return StatFys.entropy(self, temp) + boltzmann

    def free_energy(self, temp):
        """This will raise an error. Use helmholtz_free_energy or gibbs_free_energy.

           At this level 'The Free Energy' is an ambiguous concept. The
           Helmholtz Free Energy and the Gibbs Free Energy are clearly defined.

           When one wants to compute a cobtribution to the Helmholtz or Gibbs
           Free Energy from a specific factor in the partition function, then
           one must use the method ``free_energy`` of that contribution::

           >>> pf = PartFun(...)
           >>> print pf.vibrations.free_energy(300)
           >>> print pf.gibbs_free_energy(300)
        """
        raise NotImplementedError("Please use helmholtz_free_energy or gibbs_free_energy.")

    def helmholtz_free_energy(self, temp):
        """Computes the Helmholtz free energy

           Arguments:
            | temp  --  the temperature
        """
        # Similar to the extra term in the entropy, there is also an extra term
        # here due to the N factorial in the denominator of the partition
        # function.
        # The molecular ground state energy is also added here. It is tempting
        # to include it in the electronic part of partition function.
        return StatFys.free_energy(self, temp) - boltzmann*temp + self.energy

    def gibbs_free_energy(self, temp):
        """Compute the Gibbs free energy

           Arguments:
            | temp  --  the temperature

           This quantity is also known as the Gibbs energy or the Free enthalpy.
        """
        return self.helmholtz_free_energy(temp) + self.gaslaw.pv(temp)

    def dump(self, f):
        """See :meth:`Info.dump`"""
        print >> f, "Energy at T=0K [au]: %.5f" % self.energy
        print >> f, "Energy at T=0K with zero-point if QM vibrations [au]: %.5f" % self.helmholtz_free_energy(0.0)
        print >> f, self.gaslaw.description
        print >> f, "Contributions to the partition function:"
        for term in self.terms:
            term.dump(f)

    def write_to_file(self, filename):
        """Write an extensive description of the parition function to a file

           Argument:
            | filename  --  The name of the file to write to.
        """
        f = file(filename, 'w')
        self.dump(f)
        f.close()


def compute_rate_coeff(pfs_react, pf_trans, temp, cp=True, do_log=False):
    """Computes a (forward) rate coefficient

       The implementation is based on transition state theory.

       Arguments:
         | pfs_react  --  a list of partition functions objects, one for each
                          reactant
         | pf_trans  --  the partition function of the transition state
         | temp  --  the temperature

       Optional argument:
         | cp  --  When True, the rate coefficient is compute at constant
                   pressure (default=True). When False, the rate coefficient
                   is computed at constant volume.
         | do_log  --  Return the logarithm of the rate coefficient instead of
                       just the rate coefficient itself.
    """
    if cp:
        delta_G = pf_trans.gibbs_free_energy(temp)
        delta_G -= sum(pf_react.gibbs_free_energy(temp) for pf_react in pfs_react)
        log_result = -delta_G/(boltzmann*temp)
    else:
        delta_A = pf_trans.helmholtz_free_energy(temp)
        delta_A -= sum(pf_react.helmholtz_free_energy(temp) for pf_react in pfs_react)
        log_result = -delta_A/(boltzmann*temp)
    log_result += sum(pf_react.gaslaw.helper0(temp,0) for pf_react in pfs_react)
    log_result -= pf_trans.gaslaw.helper0(temp,0)
    if do_log:
        return numpy.log(boltzmann*temp/(2*numpy.pi)) + log_result
    else:
        return boltzmann*temp/(2*numpy.pi)*numpy.exp(log_result)


def compute_equilibrium_constant(pfs_A, pfs_B, temp, cp=True, do_log=False):
    """Computes the logarithm of equilibrium constant between some reactants and
       some products

       Arguments:
         | pfs_A  --  a list of reactant partition functions
         | pfs_B  --  a list of product partition functions
         | temp  --  the temperature

       Optional argument:
         | cp  --  When True, the equilibrium constant is compute at constant
                   pressure (default=True). When False, the equilibrium constant
                   is computed at constant volume.
         | do_log  --  Return the logarithm of the equilibrium constant instead
                       of just the equilibrium constant itself.
    """
    if cp:
        delta_G = 0.0
        delta_G -= sum(pf_A.gibbs_free_energy(temp) for pf_A in pfs_A)
        delta_G += sum(pf_B.gibbs_free_energy(temp) for pf_B in pfs_B)
        log_K = -delta_G/(boltzmann*temp)
    else:
        delta_A = 0.0
        delta_A -= sum(pf_A.helmholtz_free_energy(temp) for pf_A in pfs_A)
        delta_A += sum(pf_B.helmholtz_free_energy(temp) for pf_B in pfs_B)
        log_K = -delta_A/(boltzmann*temp)

    log_K += sum(pf_A.gaslaw.helper0(temp,0) for pf_A in pfs_A)
    log_K -= sum(pf_B.gaslaw.helper0(temp,0) for pf_B in pfs_B)
    if do_log:
        return log_K
    else:
        return numpy.exp(log_K)
