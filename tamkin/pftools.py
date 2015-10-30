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
"""High level utilities for partition functions"""


import sys, numpy, types, csv

from molmod.units import kjmol, mol, kelvin, joule, centimeter
from molmod.constants import boltzmann, lightspeed

from tamkin.partf import PartFun


__all__ = ["ThermoAnalysis", "ThermoTable", "ReactionAnalysis"]


class ThermoAnalysis(object):
    """Perform a regular thermochemistry analysis."""

    def __init__(self, pf, temps):
        """
           Arguments:
            | ``pf`` -- A partition function
            | ``temps`` -- An array with temperatures to consider.

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
            ThermoTable("Internal heat", kjmol, "kJ/mol", "internal_heat", pf, temps),
            ThermoTable("Heat capacity", joule/mol/kelvin, "J/(mol*K)", "heat_capacity", pf, temps),
            ThermoTable("Free energy", kjmol, "kJ/mol", "free_energy", pf, temps),
            ThermoTable("Chemical potential", kjmol, "kJ/mol", "chemical_potential", pf, temps),
            ThermoTable("Entropy", joule/mol/kelvin, "J/(mol*K)", "entropy", pf, temps),
            ThermoTable("ln(Z_N)/N", 1.0, "1", "log", pf, temps),
            ThermoTable("1/N d ln(Z_N) / dT", 1.0/kelvin, "1/K", "logt", pf, temps),
            ThermoTable("1/N d^2 ln(Z_N) / dT^2", 1.0/kelvin**2, "1/K^2", "logtt", pf, temps),
            ThermoTable("d ln(Z_N) / dN", 1.0, "1", "logn", pf, temps),
            ThermoTable("d ln(Z_N) / dN - log(V/N)", 1.0, "mixed: 1 or ln(bohr^-dim)", "logv", pf, temps),
        ]

    def write_to_file(self, filename):
        """Write the entire thermochemistry analysis to a csv file.

           Argument:
            | ``filename`` -- the file to write the output.
        """
        with open(filename, 'w') as f:
            self.dump(f)

    def dump(self, f):
        """Write the entire thermochemistry analysis in csv format.

           Argument:
            | ``f`` -- the stream to write to.
        """
        for table in self.tables:
            table.dump(f)
            print >> f


class ThermoTable(object):
    """A thermo table, i.e. the thermochemistry analysis for one
       specific thermodynamic quantity.
    """

    def __init__(self, label, unit, unit_name, method_name, pf, temps, pf_method_name=None):
        """This object is used by the ThermoAnalysis class and should probably
           never be used directly.

           Arguments:
            | ``label`` -- a string to identify the thermodynamic quantity.
            | ``unit`` -- the conversion factor from the conventional unit to
                          atomic units
            | ``unit_name`` -- a human readable string that describes the
                               conventional unit
            | ``method_name`` -- the method of the partition function that
                                 computes the quantity of interest
            | ``temps`` -- the temperatures at which the quantity has to be
                           computed.

           Optional argument:
            | ``pf_method_name`` -- In case of the actual partition function
                                    object, this alternative method can be used
                                    compute to quantity of interest. This
                                    workaround is required due to poor naming
                                    conventions in statistical physics.

           The results are stored in an array self.data of which the columns
           correspond to the given temperatures and the rows correspond to the
           different terms in the partition function.

           The attribute self.keys is a list describing the rows, i.e. the each
           contribution from the partition function.
        """
        if pf_method_name is None:
            pf_method_name = method_name

        self.label = label
        self.unit = unit
        self.unit_name = unit_name
        self.method_name = method_name
        self.pf = pf
        self.temps = temps
        self.pf_method_name = pf_method_name

        self.keys = []
        data =  []
        for term in [pf] + pf.terms:
            self.keys.append(term.name)
            if isinstance(term, PartFun):
                method = getattr(term, pf_method_name, None)
            else:
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
        """Dumps the table in csv format

           Arguments:
            | ``f`` -- the file object to write to
        """
        c = csv.writer(f)
        c.writerow([self.label, "[%s]" % self.unit_name])
        c.writerow(["Temperatures"] + [temp for temp in self.temps])
        for key, row in zip(self.keys, self.data):
            c.writerow([key] + [value/self.unit for value in row])


class ReactionAnalysis(object):
    """A Reaction analysis object."""

    def __init__(self, kinetic_model, temp_low, temp_high, temp_step=10*kelvin):
        """
           Arguments:
            | ``kinetic_model`` -- A kinetic model object. See
                                   mod:`tamkin.chemmod`.
            | ``temp_low`` -- The lower bound of the temperature interval in
                              Kelvin.
            | ``temp_high`` -- The upper bound of the temperature interval in
                               Kelvin.

           Optional arguments:
            | ``temp_step`` -- The resolution of the temperature grid.
                               [default=10K]

           The rate constants are computed on the specified temperature grid
           and afterwards the kinetic parameters are fitted to these data. All
           the results are stored as attributes of the reaction analysis object
           and can be written to text files (method write_to_file) or plotted
           (methods plot and plot_parameters). The results from multiple
           reactions can be gathered in a single plot when this is desirable.

           The following attributes may be useful:
            | ``A`` and ``Ea`` -- The kinetic parameters in atomic units.
            | ``R2`` -- The Pearson R^2 of the fit.
            | ``temps`` -- An array with the temperature grid in Kelvin
            | ``temps_inv`` -- An array with the inverse temperatures
            | ``rate_consts`` -- the rate constants at the grid points in atomic
                                 units
            | ``ln_rate_consts`` -- the logarithm of `the rate constants in
                                    atomic units`
        """
        self.kinetic_model = kinetic_model
        self.temp_low = float(temp_low)
        self.temp_high = float(temp_high)
        self.temp_step = float(temp_step)
        self.temp_high = numpy.ceil((self.temp_high-self.temp_low)/self.temp_step)*self.temp_step+self.temp_low

        # make sure that the final temperature is included
        self.temps = numpy.arange(self.temp_low,self.temp_high+0.5*self.temp_step,self.temp_step,dtype=float)
        self.temps_inv = 1/self.temps
        self.ln_rate_consts = numpy.array([
            self.kinetic_model.rate_constant(temp, do_log=True)
            for temp in self.temps
        ])
        self.rate_consts = numpy.exp(self.ln_rate_consts)

        design_matrix = numpy.zeros((len(self.temps),2), float)
        design_matrix[:,0] = 1
        design_matrix[:,1] = -self.temps_inv/boltzmann
        expected_values = self.ln_rate_consts
        if not numpy.isfinite(expected_values).all():
            raise ValueError("non-finite rate constants. check your partition functions for errors.")
        self.hessian = numpy.dot(design_matrix.transpose(), design_matrix)
        self.parameters, SSE, rank, s = numpy.linalg.lstsq(design_matrix, self.ln_rate_consts)

        SST = ((self.ln_rate_consts - self.ln_rate_consts.mean())**2).sum()
        self.R2 = 1-SSE/SST

        self.A = numpy.exp(self.parameters[0])
        self.Ea = self.parameters[1]

        self.covariance = None # see monte_carlo method

    def dump(self, f):
        """Write the results in text format on screen or to another stream.

           Argument:
            | ``f`` -- the file object to write to.
        """
        print >> f, "Summary"
        print >> f, "A [%s] = %.5e" % (self.kinetic_model.unit_name, self.A/self.kinetic_model.unit)
        print >> f, "ln(A [a.u.]) = %.2f" % (self.parameters[0])
        print >> f, "Ea [kJ/mol] = %.2f" % (self.Ea/kjmol)
        print >> f, "R2 (Pearson) = %.2f%%" % (self.R2*100)
        print >> f
        if self.covariance is not None:
            print >> f, "Error analysis"
            print >> f, "Number of Monte Carlo iterations = %i" % self.monte_carlo_iter
            print >> f, "Relative systematic error on the frequencies = %.2f" % self.freq_error
            print >> f, "Relative systematic error on the energy = %.2f" % self.energy_error
            print >> f, "Error on A [%s] = %10.5e" % (self.kinetic_model.unit_name, numpy.sqrt(self.covariance[0,0])*self.A/self.kinetic_model.unit)
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
        if self.kinetic_model.tunneling is not None:
            self.kinetic_model.tunneling.dump(f)
        print >> f, "Reaction rate constants"
        print >> f, "    T [K]    Delta_r F [kJ/mol]      k(T) [%s]" % (self.kinetic_model.unit_name)
        for i in xrange(len(self.temps)):
            temp = self.temps[i]
            delta_free = self.kinetic_model.free_energy_change(temp)
            print >> f, "% 10.2f      %8.1f             % 10.5e" % (
                temp, delta_free/kjmol, self.rate_consts[i]/self.kinetic_model.unit
            )
        print >> f
        self.kinetic_model.dump(f)
        print >> f

    def write_to_file(self, filename):
        """Write the entire analysis to a text file.

           One argument:
            | ``filename`` -- the file to write the output.
        """
        f = file(filename, "w")
        self.dump(f)
        f.close()

    def plot_arrhenius(self, filename=None, label=None, color="red"):
        """Plot the rate constant and the fitted line.

           Optional arguments:
            | ``filename`` -- When given, the plot is written to that file,
                              otherwise this plot method can be called multiple
                              times with different reaction analysis objects to
                              put all the results in one plot.
            | ``label`` -- When multiple fits are put in one figure, this label
                           is used distinguish between the various results with
                           a legend.
            | ``color`` -- Determines the color of the plotted data points and
                           line. [default="red"]. Common color names, html codes
                           and RGB tuples are accepted. (See matplotlib docs for
                           more info.)
        """
        import matplotlib.pyplot as pt

        temps_inv_line = numpy.linspace(self.temps_inv.min(),self.temps_inv.max(),100)
        ln_rate_consts_line = self.parameters[0] - self.parameters[1]/boltzmann*temps_inv_line

        if filename is not None:
            pt.clf()
            pt.title('Arrhenius plot')
        pt.text(
            0.05, 0.05,
            "A [%s] = %.3e\nEa [kJ/mol] = %.1f\nR^2 [%%] = %.1f" % (
                self.kinetic_model.unit_name, self.A/self.kinetic_model.unit,
                self.Ea/kjmol, self.R2*100
                ),
            transform=pt.gca().transAxes
        )
        pt.xlabel("1/T [1/K]")
        pt.ylabel("Rate coefficient [%s]" % self.kinetic_model.unit_name)
        if label is None:
            label_fit = "Fitted line"
            label_data = "Computed values"
        else:
            label_fit = label
            label_data = "_nolegend_"
        pt.plot(
            temps_inv_line,numpy.exp(ln_rate_consts_line)/self.kinetic_model.unit,
            color=color, linestyle="-", marker="None",label=label_fit
        )
        pt.plot(
            self.temps_inv,numpy.exp(self.ln_rate_consts)/self.kinetic_model.unit,
            color=color, linestyle="None", marker="o",label=label_data
        )
        pt.semilogy()
        if label is None:
            pt.legend(loc=0)
        if filename is not None:
            pt.savefig(filename)


    def monte_carlo(self, freq_error=1*(lightspeed/centimeter), energy_error=0.00, num_iter=100):
        """Estimate the uncertainty on the parameters

           The uncertainties are modeled by stochastic errors in the frequencies
           and the barrier. These may be caused by limitiated in convergence of
           the wafunction, limited convergence of the geometry optimization,
           errors in the numerical integration in DFT methods, or numerical
           errors due to the computation of the Hessian with a finite difference
           method.

           Optional argument:
            | ``freq_error`` -- The with of the absolute gaussian distortion on
                                the frequencies [default=1*invcm]
            | ``energy_error`` -- The width of the relative gaussian error on
                                  the energy barrier [default=0.00]
            | ``num_iter`` -- The number of Monte Carlo iterations
                              [default=1000]
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


        self.kinetic_model.backup_freqs()

        solutions = numpy.zeros((num_iter, 2), float)
        for i in xrange(num_iter):
            scale_energy = 1.0 + numpy.random.normal(0.0, 1.0)*energy_error
            self.kinetic_model.alter_freqs(freq_error, scale_energy)
            altered_ra = ReactionAnalysis(
                self.kinetic_model, self.temp_low, self.temp_high,
                self.temp_step
            )
            solutions[i] = altered_ra.parameters

        self.monte_carlo_samples = solutions.copy()
        solutions -= self.parameters
        self.covariance = numpy.dot(solutions.transpose(), solutions)/num_iter

        self.kinetic_model.restore_freqs()

    def plot_parameters(self, filename=None, label=None, color="red", marker="o", error=True):
        """Plot the kinetic parameters.

           Optional arguments:
            | ``filename`` -- When given, the plot is written to that file,
                              otherwise this plot method can be called multiple
                              times with different reaction analysis objects to
                              put all the results in one plot.
            | ``label`` -- When multiple fits are put in one figure, this label
                           is used distinguish between the various results with
                           a legend.
            | ``color`` -- Determines the color of the plotted data points and
                           line. [default="red"]. Common color names, html codes
                           and RGB tuples are accepted. (See matplotlib docs for
                           more info.)
            | ``marker`` -- The marker used for the (original) fitted parameters
                            [default="o"] (See matplotlib docs for more info.)
            | ``error`` -- A boolean that determines whether the monte carlo
                           results are plotted when they are available.
                           [default=True]
        """
        import matplotlib.pyplot as pt

        if filename is not None:
            pt.clf()
            pt.title("Parameter plot")
        pt.text(
            0.05, 0.05,
            "A [%s] = %.3e\nEa [kJ/mol] = %.1f\nR^2 [%%] = %.1f" % (
                self.kinetic_model.unit_name, self.A/self.kinetic_model.unit,
                self.Ea/kjmol, self.R2*100
                ),
            transform=pt.gca().transAxes
        )
        pt.xlabel("E_a [kJ/mol]")
        pt.ylabel("ln(A) [ln(%s)]" % self.kinetic_model.unit_name)
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
            pt.plot(
                (self.Ea + data[1])/kjmol,
                self.parameters[0] + data[0] - numpy.log(self.kinetic_model.unit),
                color=color, linestyle="-", marker="None",label=label_error
            )
            pt.plot(
                self.monte_carlo_samples[:,1]/kjmol,
                self.monte_carlo_samples[:,0] - numpy.log(self.kinetic_model.unit),
                color=color, marker=".", label=label_scatter, linestyle="None",
                markersize=1.2
            )
        pt.plot([self.Ea/kjmol],[numpy.log(self.A/self.kinetic_model.unit)], color=color,
                   marker=marker, label=label_point, mew=2, mec="white", ms=10)
        if label is None:
            pt.legend(loc=0, numpoints=1)
        if filename is not None:
            pt.savefig(filename)
