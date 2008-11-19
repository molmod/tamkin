# TAMkin is a post-processing toolkit for thermochemistry and kinetics analysis.
# Copyright (C) 2008 Toon Verstraelen <Toon.Verstraelen@UGent.be>,
# Matthias Vandichel <Matthias.Vandichel@UGent.be> and
# An Ghysels <An.Ghysels@UGent.be>
#
# This file is part of TAMkin.
#
# TAMkin is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
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


from tamkin.partf import compute_rate_coeff

from molmod.units import kjmol, second, meter, mol
from molmod.constants import boltzmann

import numpy, pylab


__all__ = ["ReactionAnalysis"]


class ReactionAnalysis(object):
    def __init__(self, pfs_react, pf_trans, temp_low, temp_high, mol_volume=None, num_temp=100):
        if len(pfs_react) == 0:
            raise ValueError("At least one reactant must be given.")
        self.pfs_react = pfs_react
        self.pf_trans = pf_trans
        self.temp_low = temp_low
        self.temp_high = temp_high
        self.mol_volume = mol_volume
        self.num_temp = num_temp

        self.temps = numpy.arange(num_temp,dtype=float)/(num_temp-1)*(temp_high-temp_low) + temp_low
        self.rate_coeffs = numpy.array([
            compute_rate_coeff(pfs_react, pf_trans, temp, mol_volume)
            for temp in self.temps
        ])

        self.temps_inv = 1/numpy.array(self.temps,float)
        self.ln_rate_coeffs = numpy.log(self.rate_coeffs)

        design_matrix = numpy.zeros((len(self.temps),2), float)
        design_matrix[:,0] = 1
        design_matrix[:,1] = -self.temps_inv
        self.parameters, residuals, rank, s = numpy.linalg.lstsq(design_matrix, self.ln_rate_coeffs)

        SSE = (residuals**2).sum()
        SST = ((self.ln_rate_coeffs - self.ln_rate_coeffs.mean())**2).sum()
        self.R2 = 1-SSE/SST

        self.A = numpy.exp(self.parameters[0])
        self.Ea = self.parameters[1]*(boltzmann)

        self.unit = (meter**3/mol)**(len(self.pfs_react)-1)/second
        if len(self.pfs_react)==1:
            self.unit_name = "1/s"
        elif len(self.pfs_react)==2:
            self.unit_name = "(m**3/mol)/s"
        else:
            self.unit_name = "(m**3/mol)**%i/s" % (len(self.pfs_react)-1)

    def dump(self, f):
        print >> f, "number of reactants: %i" % len(self.pfs_react)
        print >> f, "A [%s] = %10.5e" % (self.unit_name, self.A/self.unit)
        print >> f, "Ea [kJ/mol] = %10.2f" % (self.Ea/kjmol)
        print >> f
        print >> f, "T_low [K] = %.1f" % self.temp_low
        print >> f, "T_high [K] = %.1f" % self.temp_high
        print >> f, "number of T steps = %i" % self.num_temp
        print >> f, "R2 (Pearson) = %.2f%%" % (self.R2*100)
        print >> f
        delta_E = self.pf_trans.molecule.energy - sum(pf_react.molecule.energy for pf_react in self.pfs_react)
        print >> f, "Delta E [kJ/mol] = %.1f" % (delta_E/kjmol)
        print >> f
        print >> f, "Reaction rate coefficients"
        print >> f, "    T [K]     Delta G [kJ/mol]       k(T) [%s]" % self.unit_name
        for i in xrange(self.num_temp):
            temp = self.temps[i]
            delta_G = self.pf_trans.free_energy(temp) - sum(pf_react.free_energy(temp) for pf_react in self.pfs_react)
            print >> f, "% 10.2f      %8.1f             % 10.5e" % (
                temp, delta_G/kjmol, self.rate_coeffs[i]/self.unit
            )
        print >> f
        for counter, pf_react in enumerate(self.pfs_react):
            print >> f, "Reactant %i partition function" % counter
            pf_react.dump(f)
            print >> f
        print >> f, "Transition state partition function"
        self.pf_trans.dump(f)
        print >> f

    def write_to_file(self, filename):
        f = file(filename, "w")
        self.dump(f)
        f.close()

    def plot(self, filename):
        temps_inv_line = numpy.linspace(self.temps_inv.min(),self.temps_inv.max(),100)
        ln_rate_coeffs_line = self.parameters[0] - self.parameters[1]*temps_inv_line

        pylab.clf()
        pylab.title("Arrhenius plot: A [%s] = %.5e    Ea [kJ/mol] = %.2f" % (
            self.unit_name, self.A/self.unit, self.Ea/kjmol
        ))
        pylab.xlabel("1/T [1/K]")
        pylab.ylabel("Rate coefficient [%s]" % self.unit_name)
        pylab.plot(temps_inv_line,numpy.exp(ln_rate_coeffs_line)/self.unit,"r-",label="Fitted curve")
        pylab.plot(self.temps_inv,numpy.exp(self.ln_rate_coeffs)/self.unit,"ro",label="Computed values")
        pylab.semilogy()
        pylab.legend(loc=0)
        pylab.savefig(filename)


