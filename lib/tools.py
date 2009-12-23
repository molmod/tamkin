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


from tamkin.partf import compute_rate_coeff

from molmod.units import kjmol, second, meter, mol, K, J
from molmod.constants import boltzmann

import sys, numpy, pylab, types


__all__ = ["ThermoAnalysis", "ThermoTable", "ReactionAnalysis"]


class ThermoAnalysis(object):
    def __init__(self, pf, temps):
        """Perform a regular thermochemistry analysis.

           Arguments:
             pf  --  A partition function
             temps  --  An array with temperatures to consider.

           The tables with energy, free energy, heat capacity, entropy,
           logarithm of the partition function and the first and second order
           derivative of the logarithm of the partition functions are computed
           and stored in self.tables. The latter attribute is a list of
           ThermoTable objects.

           The results can be written to a csv file with the method
           write_to_file.
        """
        self.pf = pf
        self.temps = temps
        self.tables = [
            ThermoTable("Energy", "%.5f", kjmol, "kJ/mol", "internal_energy", pf, temps),
            ThermoTable("Free energy", "%.5f", kjmol, "kJ/mol", "free_energy", pf, temps),
            ThermoTable("Heat capacity", "%.5f", J/mol/K, "J/(mol*K)", "heat_capacity", pf, temps),
            ThermoTable("Entropy", "%.5f",  J/mol/K, "J/(mol*K)", "entropy", pf, temps),
            ThermoTable("log(q)", "%.1f", 1, "1", "log_eval", pf, temps),
            ThermoTable("d log(q) / dT", "%.3e", 1/K, "1/K", "log_deriv", pf, temps),
            ThermoTable("d^2 log(q) / dT^2", "%.1e", 1/K, "1//K", "log_deriv2", pf, temps)
        ]

    def write_to_file(self, filename):
        """Write the entire thermochemistry analysis to a csv file.

           Argument:
             filename  --  the file to write the output.
        """
        f = file(filename, "w")
        self.dump(f)
        f.close()

    def dump(self, f=sys.stdout):
        """Write the entire thermochemistry analysis to screen or to a stream in csv format.

           Optional argument:
             f  --  the stream to write to. [default=sys.stdout]
        """
        for table in self.tables:
            table.dump(f)
            print >> f


class ThermoTable(object):
    def __init__(self, label, format, unit, unit_name, method_name, pf, temps):
        """Initialize a thermo table, i.e. the thermochemistry analysis for one
           specific thermodynamic quantity.

           This object is used by the ThermoAnalysis class and should probably
           never be used directly.

           Arguments:
             label  --  a string to identify the thermodynamic quantity.
             format  --  the floating point format, e.g. "%.3f"
             unit  --  the conversion factor from the conventional unit to atomic units
             unit_name  --  a human readable string that describes the conventional unit
             method_name  --  the method of the partition function that computes the quantity of interest
             temps  --  the temperatures at which the quantity has to be computed.

           The results are stored in an array self.data of which the columns
           correspond to the given temperatures and the rows correspond to the
           different terms in the partition function.

           The attribute self.keys is a list describing the rows, i.e. the each
           contribution from the partition function.
        """
        self.label = label
        self.format = format
        self.unit = unit
        self.unit_name = unit_name
        self.method_name = method_name
        self.pf = pf
        self.temps = temps

        self.keys = []
        data =  []
        for term in [pf] + pf.terms:
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
        print >> f, '"Temperatures",%s' % ",".join("%.1f" % temp for temp in self.temps)
        for key, row in zip(self.keys, self.data):
            print >> f, '"%s",%s' % (key, ",".join(self.format % (value/self.unit) for value in row))


class ReactionAnalysis(object):
    def __init__(self, pfs_react, pf_trans, temp_low, temp_high, temp_step=10*K, mol_volume=None, tunneling=None):
        """Initialize a reaction analysis object

           Arguments:
             pfs_react  --  a list of partition functions for the reactants
             pf_trans  --  the partition function of the transition state
             temp_low  --  the lower bound of the temperature interval in Kelvin
             temp_high  --  the upper bound of the temperature interval in
                            Kelvin

           Optional arguments:
             temp_step  --  The resolution of the temperature grid.
                            [default=10K]
             mol_volume  --  An object that computes the molecular volume as
                             function of the temperature. When not given, a
                             constant volume corresponding to normal conditions
                             is assumed.
             tunneling  --  A tunneling correction object. If not given, no
                            tunneling correction is applied.

           The rate coefficients are computed on the specified temperature grid
           and afterwards the kinetic parameters are fitted to these data. All
           the results are stored as attributes of the reaction analysis object
           and can be written to text files (method write_to_file) or plotted
           (methods plot and plot_parameters). The results from multiple
           reactions can be gathered in a single plot when this is desirable.

           The following attributes might be useful:
             self.A and self.Ea  --  The kinetic parameters in atomic units.
             self.unit_name  --  A string describing the conventional unit of
                                 self.A
             self.unit  --  The conversion factor to transform self.A into
                            conventional units (self.A/self.unit)
             self.R2  --  The Pearson R^2 of the fit.
             self.temps  --  An array with the temperature grid in Kelvin
             self.temps_inv  --  An array with the inverse temperatures
             self.ln_rate_coeffs  -- the logarithm of 'the rate coefficients in
                                     atomic units'
        """
        if len(pfs_react) == 0:
            raise ValueError("At least one reactant must be given.")
        self.pfs_react = pfs_react
        self.pf_trans = pf_trans
        self.temp_low = float(temp_low)
        self.temp_high = float(temp_high)
        self.temp_step = float(temp_step)
        self.temp_high = numpy.ceil((self.temp_high-self.temp_low)/self.temp_step)*self.temp_step+self.temp_low
        self.mol_volume = mol_volume
        self.tunneling = tunneling

        # make sure that the final temperature is included
        self.temps = numpy.arange(self.temp_low,self.temp_high+0.5*self.temp_step,self.temp_step)
        self.rate_coeffs = numpy.array([
            compute_rate_coeff(pfs_react, pf_trans, temp, mol_volume)
            for temp in self.temps
        ])
        if self.tunneling is None:
            self.corrections = None
        else:
            self.corrections = self.tunneling(self.temps)
            self.rate_coeffs *= self.corrections

        self.temps_inv = 1/numpy.array(self.temps,float)
        self.ln_rate_coeffs = numpy.log(self.rate_coeffs)

        design_matrix = numpy.zeros((len(self.temps),2), float)
        design_matrix[:,0] = 1
        design_matrix[:,1] = -self.temps_inv/boltzmann
        expected_values = self.ln_rate_coeffs
        self.hessian = numpy.dot(design_matrix.transpose(), design_matrix)
        self.parameters, SSE, rank, s = numpy.linalg.lstsq(design_matrix, self.ln_rate_coeffs)

        SST = ((self.ln_rate_coeffs - self.ln_rate_coeffs.mean())**2).sum()
        self.R2 = 1-SSE/SST

        self.A = numpy.exp(self.parameters[0])
        self.Ea = self.parameters[1]

        self.unit = (meter**3/mol)**(len(self.pfs_react)-1)/second
        if len(self.pfs_react)==1:
            self.unit_name = "1/s"
        elif len(self.pfs_react)==2:
            self.unit_name = "(m**3/mol)/s"
        else:
            self.unit_name = "(m**3/mol)**%i/s" % (len(self.pfs_react)-1)

        self.covariance = None # see monte_carlo method

    def compute_rate_coeff(self, temp):
        return compute_rate_coeff(self.pfs_react, self.pf_trans, temp, self.mol_volume)

    def compute_delta_G(self, temp):
        return self.pf_trans.free_energy(temp) - \
               sum(pf_react.free_energy(temp) for pf_react in self.pfs_react)

    def compute_delta_E(self):
        return self.pf_trans.energy - \
               sum(pf_react.energy for pf_react in self.pfs_react)

    def dump(self, f=sys.stdout):
        """Write the results in text format on screen or to another stream.

           Optional argument:
             f  --  the file object to write to. This is by default the standard
                    output.
        """
        print >> f, "Summary"
        print >> f, "number of reactants: %i" % len(self.pfs_react)
        print >> f, "A [%s] = %.5e" % (self.unit_name, self.A/self.unit)
        print >> f, "ln(A [a.u.]) = %.2f" % (self.parameters[0])
        print >> f, "Ea [kJ/mol] = %.2f" % (self.Ea/kjmol)
        delta_E = self.compute_delta_E()
        print >> f, "Delta E at T=0K (classical nuclei) [kJ/mol] = %.1f" % (delta_E/kjmol)
        delta_G0K = self.compute_delta_G(0.0)
        print >> f, "Delta E at T=0K (partition function) [kJ/mol] = %.1f" % (delta_G0K/kjmol)
        print >> f, "R2 (Pearson) = %.2f%%" % (self.R2*100)
        print >> f
        if self.covariance is not None:
            print >> f, "Error analysis"
            print >> f, "Number of Monte Carlo iterations = %i" % self.monte_carlo_iter
            print >> f, "Relative systematic error on the frequencies = %.2f" % self.freq_error
            print >> f, "Relative systematic error on the energy = %.2f" % self.energy_error
            print >> f, "Error on A [%s] = %10.5e" % (self.unit_name, numpy.sqrt(self.covariance[0,0])*self.A/self.unit)
            print >> f, "Error on ln(A [a.u.]) = %.2f" % numpy.sqrt(self.covariance[0,0])
            print >> f, "Error on Ea [kJ/mol] = %.2f" % (numpy.sqrt(self.covariance[1,1])/kjmol)
            print >> f, "Parameter correlation = %.2f" % (self.covariance[0,1]/numpy.sqrt(self.covariance[0,0]*self.covariance[1,1]))
            print >> f
        print >> f, "Temperature grid"
        print >> f, "T_low [K] = %.1f" % self.temp_low
        print >> f, "T_high [K] = %.1f" % self.temp_high
        print >> f, "T_step [K] = %.1f" % self.temp_step
        print >> f, "Number of temperatures = %i" % len(self.temps)
        print >> f
        print >> f, "Reaction rate coefficients"
        if self.tunneling is None:
            print >> f, "    T [K]     Delta G [kJ/mol]       k(T) [%s]" % self.unit_name
        else:
            print >> f, "    T [K]     Delta G [kJ/mol]       k(T) [%s]         cor(T) [1]         k_tun(T) [%s]" % (self.unit_name, self.unit_name)
        for i in xrange(len(self.temps)):
            temp = self.temps[i]
            delta_G = self.pf_trans.free_energy(temp) - sum(pf_react.free_energy(temp) for pf_react in self.pfs_react)
            if self.tunneling is None:
                print >> f, "% 10.2f      %8.1f             % 10.5e" % (
                    temp, delta_G/kjmol, self.rate_coeffs[i]/self.unit
                )
            else:
                print >> f, "% 10.2f      %8.1f             % 10.5e       % 10.5e       % 10.5e" % (
                    temp, delta_G/kjmol,
                    self.rate_coeffs[i]/self.unit/self.corrections[i], # division undoes the correction
                    self.corrections[i], # the correction factor
                    self.rate_coeffs[i]/self.unit, # with correction
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
        """Write the entire analysis to a text file.

           One argument:
             filename  --  the file to write the output.
        """
        f = file(filename, "w")
        self.dump(f)
        f.close()

    def plot(self, filename=None, label=None, color="red"):
        """Plot the rate coefficients and the fitted line.

           Optional arguments:
             filename  --  When given, the plot is written to that file, other-
                           wise this plot method can be called multiple times
                           with different reaction analysis objects to put all
                           the results in one plot.
             label -- When multiple fits are put in one figure, this label is
                      used distinguish between the various results with a
                      legend.
             color -- Determines the color of the plotted data points and line.
                      [default="red"]. Common color names, html codes and RGB
                      tuples are accepted. (See matplotlib docs for more info.)
        """
        temps_inv_line = numpy.linspace(self.temps_inv.min(),self.temps_inv.max(),100)
        ln_rate_coeffs_line = self.parameters[0] - self.parameters[1]/boltzmann*temps_inv_line

        if filename is not None:
            pylab.clf()
            pylab.title("Arrhenius plot: A [%s] = %.3e    Ea [kJ/mol] = %.2f" % (
                self.unit_name, self.A/self.unit, self.Ea/kjmol
            ))
        pylab.xlabel("1/T [1/K]")
        pylab.ylabel("Rate coefficient [%s]" % self.unit_name)
        if label is None:
            label_fit = "Fitted line"
            label_data = "Computed values"
        else:
            label_fit = label
            label_data = "_nolegend_"
        pylab.plot(
            temps_inv_line,numpy.exp(ln_rate_coeffs_line)/self.unit,
            color=color, linestyle="-", marker="None",label=label_fit
        )
        pylab.plot(
            self.temps_inv,numpy.exp(self.ln_rate_coeffs)/self.unit,
            color=color, linestyle="None", marker="o",label=label_data)
        pylab.semilogy()
        if label is None:
            pylab.legend(loc=0)
        if filename is not None:
            pylab.savefig(filename)


    def monte_carlo(self, freq_error=0.05, energy_error=0.00, num_iter=100):
        """Estimate the uncertainty on the parameters due to a systematic error
           in the frequencies.

           Optional argument:
             freq_error  --  The systematic relative error on the frequencies
                             [default=0.05]
             energy_error  --  The systematic relative error on the energies
                               [default=0.00]
             num_iter  --  The number of Monte Carlo iterations [default=1000]

        """
        if freq_error < 0.0 or freq_error >= 1.0:
            raise ValueError("The argument freq_error must be in the range [0,1[.")
        if energy_error < 0.0 or freq_error >= 1.0:
            raise ValueError("The argument energy_error must be in the range [0,1[.")
        if num_iter <= 0:
            raise ValueError("The number of iterations must be strictly positive.")
        self.freq_error = freq_error
        self.energy_error = energy_error
        self.monte_carlo_iter = num_iter

        def backup_freqs(pf):
            pf.vibrational.positive_freqs_orig = pf.vibrational.positive_freqs.copy()
            pf.vibrational.negative_freqs_orig = pf.vibrational.negative_freqs.copy()
            pf.energy_backup = pf.energy

        def alter_freqs(pf, scale_freq, scale_energy):
            pf.vibrational.positive_freqs = pf.vibrational.positive_freqs_orig*scale_freq
            pf.vibrational.negative_freqs = pf.vibrational.negative_freqs_orig*scale_freq
            pf.energy = pf.energy_backup*scale_energy

        def restore_freqs(pf):
            pf.vibrational.positive_freqs = pf.vibrational.positive_freqs_orig
            pf.vibrational.negative_freqs = pf.vibrational.negative_freqs_orig
            pf.energy = pf.energy_backup
            del pf.vibrational.positive_freqs_orig
            del pf.vibrational.negative_freqs_orig
            del pf.energy_backup

        for pf_react in self.pfs_react:
            backup_freqs(pf_react)
        backup_freqs(self.pf_trans)

        solutions = numpy.zeros((num_iter, 2), float)
        for i in xrange(num_iter):
            scale_freq = 1.0 + numpy.random.normal(0.0, 1.0)*freq_error
            scale_energy = 1.0 + numpy.random.normal(0.0, 1.0)*energy_error
            for pf_react in self.pfs_react:
                alter_freqs(pf_react, scale_freq, scale_energy)
            alter_freqs(self.pf_trans, scale_freq, scale_energy)
            altered_ra = ReactionAnalysis(
                self.pfs_react, self.pf_trans, self.temp_low, self.temp_high,
                self.temp_step, self.mol_volume, self.tunneling
            )
            solutions[i] = altered_ra.parameters

        self.monte_carlo_samples = solutions.copy()
        solutions -= self.parameters
        self.covariance = numpy.dot(solutions.transpose(), solutions)/num_iter

        for pf_react in self.pfs_react:
            restore_freqs(pf_react)
        restore_freqs(self.pf_trans)

    def plot_parameters(self, filename=None, label=None, color="red", error=True):
        """Plot the kinetic parameters.

           Optional arguments:
             filename  --  When given, the plot is written to that file, other-
                           wise this plot method can be called multiple times
                           with different reaction analysis objects to put all
                           the results in one plot.
             label -- When multiple fits are put in one figure, this label is
                      used distinguish between the various results with a
                      legend.
             color -- Determines the color of the plotted data points and line.
                      [default="red"]. Common color names, html codes and RGB
                      tuples are accepted. (See matplotlib docs for more info.)
             error -- A boolean that determines whether the monte carlo results
                      are plotted when they are available. [default=True]
        """
        if filename is not None:
            pylab.clf()
            pylab.title("Parameter plot: A [%s] = %.3e    Ea [kJ/mol] = %.2f" % (
                self.unit_name, self.A/self.unit, self.Ea/kjmol
            ))
        pylab.xlabel("E_a [kJ/mol]")
        pylab.ylabel("A [%s]" % self.unit_name)
        if label is None:
            label_point = "Optimal parameters"
            label_error = "Relative uncertainty"
            label_scatter = "Monte Carlo samples"
        else:
            label_point = label
            label_error = "_nolegend_"
            label_scatter = "_nolegend_"

        if self.covariance is not None and error:
            ## Let us assume a relative Gaussian error of 400% on the k values.
            ## This is ab absolute Gaussian error of ln(4) on ln(k).
            ## We divide the Hessian with ln(4)**2 and then invert it to
            ## obtain the covariance matrix of the parameters.
            #error = numpy.log(4.0)
            #covar = numpy.linalg.inv(self.hessian/error**2)
            evals, evecs = numpy.linalg.eigh(self.covariance)
            angles = numpy.arange(0.0,360.5,1.0)/180*numpy.pi
            data = numpy.outer(evecs[:,0],numpy.cos(angles))*numpy.sqrt(evals[0]) + \
                   numpy.outer(evecs[:,1],numpy.sin(angles))*numpy.sqrt(evals[1])
            pylab.plot(
                (self.Ea + data[1])/kjmol,
                (numpy.exp(self.parameters[0] + data[0]))/self.unit,
                color=color, linestyle="-", marker="None",label=label_error
            )
            pylab.plot(
                self.monte_carlo_samples[:,1]/kjmol,
                numpy.exp(self.monte_carlo_samples[:,0])/self.unit,
                color=color, marker=".", label=label_scatter, linestyle="None",
                markersize=2.0
            )
        pylab.plot([self.Ea/kjmol],[self.A/self.unit],color=color, marker="o",label=label_point)
        pylab.semilogy()
        if label is None:
            pylab.legend(loc=0)
        if filename is not None:
            pylab.savefig(filename)


