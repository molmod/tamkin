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

from molmod.units import kjmol, second, meter, mol, K
from molmod.constants import boltzmann

import sys, numpy, pylab, types


__all__ = ["ThermoTable", "ThermoAnalysis", "ReactionAnalysis"]


class ThermoTable(object):
    def __init__(self, label, format, unit, unit_name, method_name, pf, temps):
        self.label = label
        self.format = format
        self.unit = unit
        self.unit_name = unit_name
        self.method_name = method_name
        self.pf = pf
        self.temps = temps

        self.keys = []
        data =  []
        for term in (pf,) + pf.contributions:
            self.keys.append(term.name)
            method = getattr(term, method_name, None)
            if isinstance(method, types.MethodType):
                row = []
                for temp in temps:
                    row.append(method(temp))
                data.append(numpy.array(row,ndmin=2))
            method = getattr(term, "%s_terms" % method_name, None)
            if isinstance(method, types.MethodType):
                columns = []
                for i in xrange(term.num_terms):
                    self.keys.append("%s (%i)" % (term.name, i))
                for temp in temps:
                    columns.append(method(temp))
                data.append(numpy.array(columns).transpose())
        self.data = numpy.concatenate(data)

    def dump(self, f):
        """Dumps the table in csv format"""
        print >> f, '"%s","[%s]"' % (self.label, self.unit_name)
        print >> f, '"Temperatures",', ",".join("%.1f" % temp for temp in self.temps)
        for key, row in zip(self.keys, self.data):
            print >> f, '"%s",' % key, ",".join(self.format % (value/self.unit) for value in row)


class ThermoAnalysis(object):
    def __init__(self, pf, temps):
        from molmod.units import kjmol, J, mol, K
        self.pf = pf
        self.temps = temps
        self.tables = [
            ThermoTable("Energy", "%.1f", kjmol, "kJ/mol", "internal_energy", pf, temps),
            ThermoTable("Free energy", "%.1f", kjmol, "kJ/mol", "free_energy", pf, temps),
            ThermoTable("Heat capacity", "%.3f", J/mol/K, "J/(mol*K)", "heat_capacity", pf, temps),
            ThermoTable("Entropy", "%.1f",  J/mol/K, "J/(mol*K)", "entropy", pf, temps),
            ThermoTable("log(q)", "%.1f", 1, "1", "log_eval", pf, temps),
            ThermoTable("d log(q) / dT", "%.3e", 1/K, "1/K", "log_deriv", pf, temps),
            ThermoTable("d^2 log(q) / dT^2", "%.1e", 1/K, "1//K", "log_deriv2", pf, temps)
        ]

    def write_to_file(self, filename):
        f = file(filename, "w")
        self.dump(f)
        f.close()

    def dump(self, f=sys.stdout):
        for table in self.tables:
            table.dump(f)
            print >> f


class ReactionAnalysis(object):
    def __init__(self, pfs_react, pf_trans, temp_low, temp_high, mol_volume=None, temp_step=10*K):
        if len(pfs_react) == 0:
            raise ValueError("At least one reactant must be given.")
        self.pfs_react = pfs_react
        self.pf_trans = pf_trans
        self.temp_low = float(temp_low)
        self.temp_high = float(temp_high)
        self.temp_step = float(temp_step)
        self.temp_high = numpy.ceil((self.temp_high-self.temp_low)/self.temp_step)*self.temp_step+self.temp_low
        self.mol_volume = mol_volume

        # make sure that the final temperature is included
        self.temps = numpy.arange(self.temp_low,self.temp_high+0.5*self.temp_step,self.temp_step)
        self.rate_coeffs = numpy.array([
            compute_rate_coeff(pfs_react, pf_trans, temp, mol_volume)
            for temp in self.temps
        ])

        self.temps_inv = 1/numpy.array(self.temps,float)
        self.ln_rate_coeffs = numpy.log(self.rate_coeffs)

        design_matrix = numpy.zeros((len(self.temps),2), float)
        design_matrix[:,0] = 1
        design_matrix[:,1] = -self.temps_inv
        self.parameters, SSE, rank, s = numpy.linalg.lstsq(design_matrix, self.ln_rate_coeffs)

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
        print >> f, "T_step [K] = %.1f" % self.temp_step
        print >> f, "Number of steps = %i" % len(self.temps)
        print >> f, "R2 (Pearson) = %.2f%%" % (self.R2*100)
        print >> f
        delta_E = self.pf_trans.molecule.energy - sum(pf_react.molecule.energy for pf_react in self.pfs_react)
        print >> f, "Delta E [kJ/mol] = %.1f" % (delta_E/kjmol)
        print >> f
        print >> f, "Reaction rate coefficients"
        print >> f, "    T [K]     Delta G [kJ/mol]       k(T) [%s]" % self.unit_name
        for i in xrange(len(self.temps)):
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

    def plot(self, filename=sys.stdout):
        temps_inv_line = numpy.linspace(self.temps_inv.min(),self.temps_inv.max(),100)
        ln_rate_coeffs_line = self.parameters[0] - self.parameters[1]*temps_inv_line

        pylab.clf()
        pylab.title("Arrhenius plot: A [%s] = %.3e    Ea [kJ/mol] = %.2f" % (
            self.unit_name, self.A/self.unit, self.Ea/kjmol
        ))
        pylab.xlabel("1/T [1/K]")
        pylab.ylabel("Rate coefficient [%s]" % self.unit_name)
        pylab.plot(temps_inv_line,numpy.exp(ln_rate_coeffs_line)/self.unit,"r-",label="Fitted curve")
        pylab.plot(self.temps_inv,numpy.exp(self.ln_rate_coeffs)/self.unit,"ro",label="Computed values")
        pylab.semilogy()
        pylab.legend(loc=0)
        pylab.savefig(filename)


