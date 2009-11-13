# TAMkin is a post-processing toolkit for thermochemistry and kinetics analysis.
# Copyright (C) 2008-2009 Toon Verstraelen <Toon.Verstraelen@UGent.be>,
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

   The workhorse of this module is the PartFun class. PartFun objects represent
   a partition function with a interface that does not depend on the different
   types of terms that are present in the partition function. PartFun objects
   can be used to study chemical equilibrium, rate coefficients and various
   thermodynamic properties. All the applications work without explicitly relying
   on the individual contributions to the partition function.

   Abstract classes:
       Info, StatFys, StatFysTerms
   Contributions:
       Electronic, ExternalTranslation, ExternalRotation, Vibrations,
       Rotor (see rotor.py)
   Selection of the ensemble:
       IdealGasVolume, FixedVolume
       (these are used by ExternalTranslation)
   Helper functions:
       log_eval_vibrations, log_deriv_vibrations, log_deriv2_vibrations
   Applications:
       compute_rate_coeff, compute_equilibrium_constant
       (see tools.py for a friendly interface to these functions)
"""


from molmod.constants import boltzmann, lightspeed
from molmod.units import atm, bar, amu, cm

import numpy


__all__ = [
    "IdealGasVolume", "FixedVolume", "Info", "StatFys", "StatFysTerms",
    "Electronic", "ExternalTranslation", "ExternalRotation", "Vibrations",
    "log_eval_vibrations", "log_deriv_vibrations", "log_deriv2_vibrations",
    "PartFun", "compute_rate_coeff", "compute_equilibrium_constant"
]



class IdealGasVolume(object):
    """Computes the volume of a molecule in an ideal gas

       The volume is function of the temperature, at a fixed reference pressure.
       This law can be used to set up a constant pressure gas phase partition
       function.
    """

    def __init__(self, pressure=1*atm):
        self.pressure = pressure

    def __call__(self, temp):
        return boltzmann*temp/self.pressure

    def deriv(self, temp):
        return -self(temp)/temp

    def deriv2(self, temp):
        return -2*self.deriv(temp)/temp

    def get_description(self):
        return "Molecular volume: ideal gas law, reference pressure [bar] = %.5f" % (self.pressure/bar)

    description = property(get_description)


class FixedVolume(object):
    """Computes the volume of a molecule in system with fixed size

       This law can be used to set up a constant volume gas phase partition
       function.
    """

    def __init__(self, temp=298.15, pressure=1*atm):
        self.temp = temp
        self.pressure = pressure

    def __call__(self, temp):
        return boltzmann*self.temp/self.pressure

    def deriv(self, temp):
        return 0.0

    def deriv2(self, temp):
        return 0.0

    def get_description(self):
        return "Molecular volume: fixed, reference pressure [bar] = %.5f, reference temperature [K] = %.2f" % (self.pressure/bar, self.temp)

    description = property(get_description)


class Info(object):
    """An object that has a name and that can dump info to a file."""
    def __init__(self, name):
        self.name = name

    def dump(self, f):
        print >> f, "  %s" % self.name.upper()

    def dump_values(self, f, label, values, format, num_col=8):
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

       The constructor (__init__) and the first four methods (init_part_fun,
       log_eval, log_deriv and log_deriv2) must be implemented in derived
       classes.
    """
    def init_part_fun(self, nma):
        pass

    def log_eval(self, temp):
        """The logarithm of the partition function"""
        raise NotImplementedError

    def log_deriv(self, temp):
        """The derivative of the logarithm of the partition function"""
        raise NotImplementedError

    def log_deriv2(self, temp):
        """The second derivative of the logarithm of the partition function"""
        raise NotImplementedError

    def internal_energy(self, temp, log_deriv=None):
        """Computes the internal energy per molecule"""
        if log_deriv is None:
            log_deriv = self.log_deriv
        return boltzmann*temp**2*log_deriv(temp)

    def heat_capacity(self, temp, log_deriv=None, log_deriv2=None):
        """Computes the heat capacity per molecule"""
        if log_deriv is None:
            log_deriv = self.log_deriv
        if log_deriv2 is None:
            log_deriv2 = self.log_deriv2
        return boltzmann*temp*(2*log_deriv(temp) + temp*log_deriv2(temp))

    def entropy(self, temp, log_eval=None, log_deriv=None):
        """Computes the entropy contribution per molecule"""
        if log_eval is None:
            log_eval = self.log_eval
        if log_deriv is None:
            log_deriv = self.log_deriv
        return boltzmann*(log_eval(temp) + temp*log_deriv(temp))

    def free_energy(self, temp, log_eval=None):
        """Computes the free energy per molecule"""
        if log_eval is None:
            log_eval = self.log_eval
        return -boltzmann*temp*log_eval(temp)


class StatFysTerms(StatFys):
    """Abstract class for (contributions to) the parition function with multiple terms

       The constructor (__init__) and the four methods (init_part_fun,
       log_eval_terms, log_deriv_terms and log_deriv2_terms) must be implemented
       in derived classes.
    """
    def __init__(self, num_terms):
        self.num_terms = num_terms

    def log_eval(self, temp):
        return self.log_eval_terms(temp).sum()

    def log_deriv(self, temp):
        return self.log_deriv_terms(temp).sum()

    def log_deriv2(self, temp):
        return self.log_deriv2_terms(temp).sum()

    def log_eval_terms(self, temp):
        """The logarithm of the partition function (separate terms)"""
        raise NotImplementedError

    def log_deriv_terms(self, temp):
        """The derivative of the logarithm of the partition function (separate terms)"""
        raise NotImplementedError

    def log_deriv2_terms(self, temp):
        """The second derivative of the logarithm of the partition function (separate terms)"""
        raise NotImplementedError

    def internal_energy_terms(self, temp):
        """Computes the internal energy per molecule (separate terms)"""
        return self.internal_energy(temp, self.log_deriv_terms)

    def heat_capacity_terms(self, temp):
        """Computes the heat capacity per molecule (separate terms)"""
        return self.heat_capacity(temp, self.log_deriv_terms, self.log_deriv2_terms)

    def entropy_terms(self, temp):
        """Computes the entropy contribution per molecule (separate terms)"""
        return self.entropy(temp, self.log_eval_terms, self.log_deriv_terms)

    def free_energy_terms(self, temp, log_eval=None):
        """Computes the free energy per molecule (separate terms)"""
        return self.free_energy(temp, self.log_eval_terms)


class Electronic(Info, StatFys):
    """The electronic contribution to the partition function

       TODO: this should also include the potential energy from the ab initio
       computation
    """
    def __init__(self, multiplicity=None):
        self.multiplicity = multiplicity
        Info.__init__(self, "electronic")

    def init_part_fun(self, nma):
        if self.multiplicity is None:
            self.multiplicity = nma.multiplicity
            if self.multiplicity is None:
                raise ValueError("Spin multiplicity is not defined.")

    def dump(self, f):
        Info.dump(self, f)
        print >> f, "    Multiplicity: %i" % self.multiplicity

    def log_eval(self, temp):
        return numpy.log(self.multiplicity)

    def log_deriv(self, temp):
        return 0.0

    def log_deriv2(self, temp):
        return 0.0


class ExternalTranslation(Info, StatFys):
    """The contribution from the external translation"""
    def __init__(self, mol_volume=None):
        if mol_volume is None:
            self.mol_volume = FixedVolume()
        else:
            self.mol_volume = mol_volume
        Info.__init__(self, "translational")

    def init_part_fun(self, nma):
        self.mass = nma.mass

    def dump(self, f):
        Info.dump(self, f)
        print >> f, "    %s" % self.mol_volume.description
        print >> f, "    Mass [amu]: %f" % (self.mass/amu)

    def log_eval(self, temp):
        return (
            1.5*numpy.log(0.5*self.mass*boltzmann*temp/numpy.pi)+
            numpy.log(self.mol_volume(temp))
        )

    def log_deriv(self, temp):
        return 1.5/temp + self.mol_volume.deriv(temp)/self.mol_volume(temp)

    def log_deriv2(self, temp):
        V = self.mol_volume(temp)
        V_deriv = self.mol_volume.deriv(temp)
        V_deriv2 = self.mol_volume.deriv2(temp)
        return -1.5/temp**2 + (V_deriv2 - V_deriv**2/V)/V


class ExternalRotation(Info, StatFys):
    """The contribution from the external rotation"""
    def __init__(self, symmetry_number=None, im_threshold=1.0):
        self.symmetry_number = symmetry_number
        self.im_threshold = im_threshold
        Info.__init__(self, "rotational")

    def init_part_fun(self, nma):
        if nma.periodic:
            raise ValueError("There is no external rotation in periodic systems.")
        self.inertia_tensor = nma.inertia_tensor
        self.moments = numpy.linalg.eigvalsh(nma.inertia_tensor)
        if self.symmetry_number is None:
            self.symmetry_number = nma.symmetry_number
            if self.symmetry_number is None:
                raise ValueError("Symmetry number is not defined.")
        self.factor = numpy.sqrt(numpy.product([
            2*numpy.pi*m*boltzmann for m in self.moments if m > self.im_threshold
        ]))/self.symmetry_number/numpy.pi
        self.count = (self.moments > self.im_threshold).sum()

    def dump(self, f):
        Info.dump(self, f)
        print >> f, "    Rotational symmetry number: %i" % self.symmetry_number
        print >> f, "    Moments of inertia [amu*bohr**2]: %f  %f %f" % tuple(self.moments/amu)
        print >> f, "    Threshold for non-zero moments of inertia [amu*bohr**2]: %e" % (self.im_threshold/amu)
        print >> f, "    Non-zero moments of inertia: %i" % self.count

    def log_eval(self, temp):
        return numpy.log(temp)*0.5*self.count + numpy.log(self.factor)

    def log_deriv(self, temp):
        return 0.5*self.count/temp

    def log_deriv2(self, temp):
        return -0.5*self.count/temp**2


def log_eval_vibrations(temp, freqs, classical=False):
    """The logarithm of the partition function for a harmonic vibration"""
    # this is defined as a function because multiple classes need it
    if classical:
        return numpy.log(0.5*boltzmann*temp/numpy.pi/freqs)
    else:
        # The zero point correction is included in the partition function and
        # should not be taken into account when computing the reaction barrier.
        exp_arg = -2*numpy.pi*freqs/boltzmann/temp
        return (exp_arg/2 - numpy.log(1-numpy.exp(exp_arg)))
        # This would be the version when the zero point energy corrections are
        # included in the energy difference when computing the reaction rate:
        #return -numpy.log(1-numpy.exp(exp_arg))

def log_deriv_vibrations(temp, freqs, classical=False):
    """The derivative of the logarithm of the partition function for a harmonic vibration"""
    # this is defined as a function because multiple classes need it
    if classical:
        return numpy.ones(len(freqs))/temp
    else:
        exp_arg = -2*numpy.pi*freqs/boltzmann/temp
        exp_arg_deriv = -exp_arg/temp
        return exp_arg_deriv*(0.5-1/(1-numpy.exp(-exp_arg)))

def log_deriv2_vibrations(temp, freqs, classical=False):
    """The second derivative of the logarithm of the partition function for a harmonic vibration"""
    # this is defined as a function because multiple classes need it
    if classical:
        return -numpy.ones(len(freqs))/temp**2
    else:
        exp_arg = -2*numpy.pi*freqs/boltzmann/temp
        exp_arg_deriv = -exp_arg/temp
        exp_arg_deriv2 = -2*exp_arg_deriv/temp
        e = numpy.exp(-exp_arg)
        x = 1/(1-e)
        return exp_arg_deriv2*(0.5-x) + (exp_arg_deriv*x)**2*e


class Vibrations(Info, StatFysTerms):
    """The vibrational contribution to the partition function"""
    def __init__(self, classical=False):
        self.classical = classical
        Info.__init__(self, "vibrational")

    def init_part_fun(self, nma):
        zero_indexes = nma.zeros
        nonzero_mask = numpy.ones(len(nma.freqs), dtype=bool)
        nonzero_mask[zero_indexes] = False

        self.freqs = nma.freqs[nonzero_mask]
        self.zero_freqs = nma.freqs[zero_indexes]
        self.positive_freqs = nma.freqs[(nma.freqs > 0) & nonzero_mask]
        self.negative_freqs = nma.freqs[(nma.freqs < 0) & nonzero_mask]

        StatFysTerms.__init__(self, len(self.positive_freqs))

    def dump(self, f):
        Info.dump(self, f)
        print >> f, "    Number of zero wavenumbers: %i " % (len(self.zero_freqs))
        print >> f, "    Number of real wavenumbers: %i " % (len(self.positive_freqs))
        print >> f, "    Number of imaginary wavenumbers: %i" % (len(self.negative_freqs))
        self.dump_values(f, "Zero Wavenumbers [1/cm]", self.zero_freqs/(lightspeed/cm), "% 8.1f", 8)
        self.dump_values(f, "Real Wavenumbers [1/cm]", self.positive_freqs/(lightspeed/cm), "% 8.1f", 8)
        self.dump_values(f, "Imaginary Wavenumbers [1/cm]", self.negative_freqs/(lightspeed/cm), "% 8.1f", 8)

    def log_eval_terms(self, temp):
        return log_eval_vibrations(temp, self.positive_freqs, self.classical)

    def log_deriv_terms(self, temp):
        return log_deriv_vibrations(temp, self.positive_freqs, self.classical)

    def log_deriv2_terms(self, temp):
        return log_deriv2_vibrations(temp, self.positive_freqs, self.classical)


class PartFun(Info, StatFys):
    """The partition function

       This object contains all contributions to the partition function in
       self.terms and makes sure they are properly initialized. It also
       implements all the methods defined in StatFys, e.g. it can compute
       the entropy, the free energy and so on.
    """
    __reserved_names__ = set([
        "terms"
    ])

    def __init__(self, nma, terms=None):
        """Initialize the PartFun object

        Arguments:
          nma  --  NMA object
        Optional arguments:
          terms  --  list to select the contributions to the partition function
                     e.g. [Vibrations(classical=True), ExternalRotation(1)]
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

        for term in self.terms:
            term.init_part_fun(nma)

        self.energy = nma.energy
        Info.__init__(self, "total")

    def log_eval(self, temp):
        return sum(term.log_eval(temp) for term in self.terms)

    def log_deriv(self, temp):
        return sum(term.log_deriv(temp) for term in self.terms)

    def log_deriv2(self, temp):
        return sum(term.log_deriv2(temp) for term in self.terms)

    def entropy(self, temp):
        """Compute the total entropy"""
        return boltzmann + StatFys.entropy(self, temp)

    def free_energy(self, temp):
        """Computes the free energy"""
        return StatFys.free_energy(self, temp) + self.energy

    def dump(self, f):
        print >> f, "Contributions to the partition function:"
        for term in self.terms:
            term.dump(f)

    def write_to_file(self, filename):
        f = file(filename, 'w')
        self.dump(f)
        f.close()


def compute_rate_coeff(pfs_react, pf_trans, temp, mol_volume=None):
    """Computes a (forward) rate coefficient

       The implementation is based on transition state theory.

       Arguments:
         pfs_react  --  a list of partition functions objects, one for each
                        reactant
         pf_trans  --  the partition function of the transition state
         temp  --  the temperature

       Optional argument:
         mol_volume  --  function that computes the molecular volume as a
                         function of temperature. It should be the same as
                         the mol_volume function for the ExternalTranslation
                         object of all the partition functions.
    """
    delta_G = pf_trans.free_energy(temp) - sum(pf_react.free_energy(temp) for pf_react in pfs_react)
    log_result = -delta_G/(boltzmann*temp)
    if len(pfs_react) > 1:
        if mol_volume is None:
            mol_volume = FixedVolume()
        log_result += numpy.log(mol_volume(temp))*(len(pfs_react)-1)
    return boltzmann*temp/(2*numpy.pi)*numpy.exp(log_result)


def compute_equilibrium_constant(pfs_A, pfs_B, temp, mol_volume=None):
    """Computes the equilibrium constant between some reactants and some products

       Arguments:
         pfs_A  --  a list of reactant partition functions
         pfs_B  --  a list of product partition functions
         temp  --  the temperature

       Optional argument:
         mol_volume  --  function that computes the molecular volume as a
                         function of temperature. It should be the same as
                         the mol_volume function for the ExternalTranslation
                         object of all the partition functions.
    """
    delta_G = 0.0
    delta_G += sum(pf_A.free_energy(temp) for pf_A in pfs_A)
    delta_G -= sum(pf_B.free_energy(temp) for pf_B in pfs_B)
    log_K = -delta_G/(boltzmann*temp)
    if len(pfs_A) != len(pfs_B):
        if mol_volume is None:
            mol_volume = FixedVolume()
        log_K += (len(pfs_B)-len(pfs_A))*numpy.log(mol_volume(temp))
    return log_K


