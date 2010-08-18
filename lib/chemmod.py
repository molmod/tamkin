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
"""Convenient interfaces to define thermodynamic and kinetic models"""


import numpy

from molmod.units import kjmol, second, meter, mol

from tamkin.partf import compute_rate_coeff, compute_equilibrium_constant


__all__ = [
    "BaseModel", "ThermodynamicModel", "BaseKineticModel", "KineticModel",
    "ActivationKineticModel"
]


class BaseModel(object):
    """Base class for all physico-chemical models."""
    def __init__(self, pfs_all):
        """
           Argument:
            | pfs_all  --  All partition functions involved in the model.

           The methods in the base class are mainly used by the Monte Carlo
           routine in the ReactionAnalysis class.
        """
        self.pfs_all = pfs_all

    def backup_freqs(self):
        """Keep a backup copy of the frequencies and the energy of each partition function."""
        for pf in self.pfs_all:
            pf.vibrational.positive_freqs_orig = pf.vibrational.positive_freqs.copy()
            pf.vibrational.negative_freqs_orig = pf.vibrational.negative_freqs.copy()
            pf.energy_backup = pf.energy

    def alter_freqs(self, freq_error, scale_energy):
        """Randomly distort the frequencies and energies.

           Arguments:
            | freq_error  --  The absolute error to be introduced in the
                              frequencies.
            | scale_energy  --  The relative error to be introduced in the
                                (electronic) energies.
        """
        for pf in self.pfs_all:
            N = len(pf.vibrational.positive_freqs)
            freq_shift = numpy.random.normal(0, freq_error, N)
            pf.vibrational.positive_freqs = pf.vibrational.positive_freqs_orig + freq_shift
            pf.vibrational.positive_freqs[pf.vibrational.positive_freqs<=0] = 0.01
            N = len(pf.vibrational.negative_freqs)
            freq_shift = numpy.random.normal(0, freq_error, N)
            pf.vibrational.negative_freqs = pf.vibrational.negative_freqs_orig + freq_shift
            pf.vibrational.negative_freqs[pf.vibrational.negative_freqs>=0] = -0.01
            pf.energy = pf.energy_backup*scale_energy

    def restore_freqs(self):
        """Restore the backup of the frequencies and the energy of each partition function."""
        for pf in self.pfs_all:
            pf.vibrational.positive_freqs = pf.vibrational.positive_freqs_orig
            pf.vibrational.negative_freqs = pf.vibrational.negative_freqs_orig
            pf.energy = pf.energy_backup
            del pf.vibrational.positive_freqs_orig
            del pf.vibrational.negative_freqs_orig
            del pf.energy_backup

    def get_free_energy_symbol(self):
        """Return the symbol for the free energy."""
        raise NotImplementedError

    def compute_delta_free_energy(self, temp):
        """Compute the free energy difference for this chemical model.

           Arguments:
            | temp  -- The temperature.
        """
        raise NotImplementedError

    def dump(self, f):
        """Write all info about the model to a file."""
        raise NotImplementedError


class ThermodynamicModel(BaseModel):
    """A model for a thermodynamic equilibrium."""
    def __init__(self, pfs_react, pfs_prod, cp=True):
        """
           Arguments:
            | pfs_react  --  A list with reactant partition functions.
            | pfs_prod  --  A list with product partition functions.

           Optional argument:
            | cp  --  When True, the equilibrium is modeled at constant pressure,
                      otherwise at constant volume.
        """
        self.pfs_react = pfs_react
        self.pfs_prod = pfs_prod
        self.cp = cp
        BaseModel.__init__(self, pfs_react + pfs_prod)

    def get_free_energy_symbol(self):
        """Return the symbol for the free energy."""
        return {True: "G", False: "A"}[self.cp]

    def compute_equilibrium_constant(self, temp, do_log=False):
        """Compute the equilibrium constant at the given temperature.

           Argument:
            | temp  --  The temperature.

           Optional argument:
            | do_log  --  When True, the logarithm of the equilibrium constant
                          is returned instead of just the equilibrium constant
                          itself. [default=False]
        """
        return compute_equilibrium_constant(self.pfs_react, self.pfs_prod, temp, self.cp, do_log)

    def compute_delta_free_energy(self, temp):
        """Compute the free energy difference between (+) products and (-) reactants.

           Arguments:
            | temp  --  The temperature.
        """
        if self.cp:
            return self._compute_delta_G(temp)
        else:
            return self._compute_delta_A(temp)

    def _compute_delta_G(self, temp):
        """Compute the Gibbs free energy difference between (+) products and (-) reactants.

           Arguments:
            | temp  --  The temperature.
        """
        return sum(pf_prod.gibbs_free_energy(temp) for pf_prod in self.pfs_prod) - \
               sum(pf_react.gibbs_free_energy(temp) for pf_react in self.pfs_react)

    def _compute_delta_A(self, temp):
        """Compute the Helmholts free energy difference between (+) products and (-) reactants.

           Arguments:
            | temp  --  The temperature.
        """
        return sum(pf_prod.helmholtz_free_energy(temp) for pf_prod in self.pfs_prod) - \
               sum(pf_react.helmholtz_free_energy(temp) for pf_react in self.pfs_react)

    def _compute_delta_E(self):
        """Compute the classical (microscopic) energy difference between (+) products and (-) reactants."""
        return sum(pf_prod.energy for pf_prod in self.pfs_prod) - \
               sum(pf_react.energy for pf_react in self.pfs_react)

    def dump(self, f):
        """Write all info about the thermodynamic model to a file."""
        delta_E = self._compute_delta_E()
        print >> f, "Delta E at T=0K [kJ/mol] = %.1f" % (delta_E/kjmol)
        delta_A0K = self._compute_delta_A(0.0)
        print >> f, "Delta E0 at T=0K (with zero-point if QM vibrations) [kJ/mol] = %.1f" % (delta_A0K/kjmol)
        for counter, pf_react in enumerate(self.pfs_react):
            print >> f, "Reactant %i partition function" % counter
            pf_react.dump(f)
            print >> f
        for counter, pf_prod in enumerate(self.pfs_prod):
            print >> f, "Prodcut %i partition function" % counter
            pf_react.dump(f)
            print >> f


class BaseKineticModel(BaseModel):
    """A generic model for the rate of a chemical reaction."""
    def compute_rate_coeff(self, temp, do_log=False):
        """Compute the rate coefficient of the reaction in this analysis

           Arguments:
            | temp  -- The temperature.

           Optional argument:
            | do_log  --  When True, the logarithm of the rate coefficient is
                          returned instead of just the rate coefficient itself.
                          [default=False]
        """
        raise NotImplementedError


class KineticModel(BaseKineticModel):
    """A model for the rate of a single-step chemical reaction."""
    def __init__(self, pfs_react, pf_trans, cp=True, tunneling=None):
        """
           Arguments:
            | pfs_react  --  a list of partition functions for the reactants
            | pf_trans  --  the partition function of the transition state

           Optional arguments:
            | cp  --  When True, the rate coefficients are compute at constant
                      pressure [default=True]. When False, the rate coefficients
                      are computed at constant volume.
            | tunneling  --  A tunneling correction object. If not given, no
                             tunneling correction is applied.


           Useful attributes:
            | unit_name  --  A string describing the conventional unit of
                             the rate coefficient
            | unit  --  The conversion factor to transform self.A into
                        conventional units (rate/self.unit)
        """
        if len(pfs_react) == 0:
            raise ValueError("At least one reactant must be given.")
        self.pfs_react = pfs_react
        self.pf_trans = pf_trans
        self.cp = cp
        self.tunneling = tunneling

        self.unit = (meter**3/mol)**(len(self.pfs_react)-1)/second
        if len(self.pfs_react)==1:
            self.unit_name = "1/s"
        elif len(self.pfs_react)==2:
            self.unit_name = "(m**3/mol)/s"
        else:
            self.unit_name = "(m**3/mol)**%i/s" % (len(self.pfs_react)-1)
        BaseKineticModel.__init__(self, pfs_react + [pf_trans])

    def get_free_energy_symbol(self):
        """Return the symbol for the free energy."""
        return {True: "G", False: "A"}[self.cp]

    def compute_rate_coeff(self, temp, do_log=False):
        """See :meth:`BaseKineticModel.compute_rate_coeff`"""
        result = compute_rate_coeff(self.pfs_react, self.pf_trans, temp, self.cp, do_log)
        if self.tunneling is not None:
            if do_log:
                result += numpy.log(self.tunneling(temp))
            else:
                result *= self.tunneling(temp)
        return result

    def compute_delta_free_energy(self, temp):
        """Compute the free energy barrier of the reaction.

           Arguments:
            | temp  -- the temperature
        """
        if self.cp:
            return self._compute_delta_G(temp)
        else:
            return self._compute_delta_A(temp)

    def _compute_delta_G(self, temp):
        """Compute the Gibbs free energy barrier of the reaction.

           Arguments:
            | temp  -- the temperature
        """
        return self.pf_trans.gibbs_free_energy(temp) - \
               sum(pf_react.gibbs_free_energy(temp) for pf_react in self.pfs_react)

    def _compute_delta_A(self, temp):
        """Compute the Helmholtz free energy barrier of the reaction.

           Arguments:
            | temp  -- the temperature
        """
        return self.pf_trans.helmholtz_free_energy(temp) - \
               sum(pf_react.helmholtz_free_energy(temp) for pf_react in self.pfs_react)

    def _compute_delta_E(self):
        """Compute the classical (microscopic) energy barrier of the reaction."""
        return self.pf_trans.energy - \
               sum(pf_react.energy for pf_react in self.pfs_react)

    def dump(self, f):
        """Write all info about the kinetic model to a file."""
        delta_E = self._compute_delta_E()
        print >> f, "Delta E at T=0K [kJ/mol] = %.1f" % (delta_E/kjmol)
        delta_A0K = self._compute_delta_A(0.0)
        print >> f, "Delta E0 at T=0K (with zero-point if QM vibrations) [kJ/mol] = %.1f" % (delta_A0K/kjmol)
        for counter, pf_react in enumerate(self.pfs_react):
            print >> f, "Reactant %i partition function" % counter
            pf_react.dump(f)
            print >> f
        print >> f, "Transition state partition function"
        self.pf_trans.dump(f)


class ActivationKineticModel(BaseKineticModel):
    """A model for the rate of a single-step chemical reaction with a pre-reactive complex."""
    def __init__(self, tm, km):
        """
           Arguments:
            | tm  --  The thermodynamic model for the pre-reactive complex.
            | km  --  The kinetic model for the single-step reaction.
        """
        self.tm = tm
        self.km = km
        assert(km.cp==tm.cp)
        BaseKineticModel.__init__(self, tm.pfs_all + km.pfs_all)

    def get_free_energy_symbol(self):
        """Return the symbol for the free energy."""
        return self.km.get_free_energy_symbol()

    def compute_delta_free_energy(self, temp):
        """Compute the free energy barrier of the entire reaction.

           Arguments:
            | temp  -- the temperature
        """
        return self.tm.compute_delta_free_energy(temp) + \
               self.km.compute_delta_free_energy(temp)

    def dump(self, f):
        """Write all info about the kinetic model to a file."""
        print >> f, "Thermodynamic model"
        self.tm.dump(f)
        print >> f, "Kinetic model"
        self.km.dump(f)

    def compute_rate_coeff(self, temp, do_log=False):
        """See :meth:`BaseKineticModel.compute_rate_coeff`"""
        if do_log:
            return self.tm.compute_rate_coeff(temp, True) + \
                   self.km.compute_rate_coeff(temp, True)
        else:
            return self.tm.compute_rate_coeff(temp, False)* \
                   self.km.compute_rate_coeff(temp, False)
