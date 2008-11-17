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


from molmod.units import kjmol, second
from molmod.constants import boltzmann

import numpy, pylab

#from pylab import *
#from scipy import *


# this function fits kin_coeffs to the Arrhenius Equation
def fit_kin(temp,kin_coeffs,prefix="arrhenius"): # args=list of temperatures and list of kinetic coeffs
    temp_inv = 1/numpy.array(temp,float)
    ln_kin_coeffs = numpy.log(kin_coeffs/(1/second))

    design_matrix = numpy.zeros((len(temp),2), float)
    design_matrix[:,0] = 1
    design_matrix[:,1] = temp_inv
    parameters, residus, rank, s = numpy.linalg.lstsq(design_matrix, ln_kin_coeffs)

    A = parameters[0]
    Ea = parameters[1]*(-boltzmann)/kjmol
    #print "A [1/s] = %10.5e" % A
    #print "Ea [kJ/mol] = %10.2f" % Ea

    temp_inv_line = numpy.linspace(temp_inv.min(),temp_inv.max(),100)
    ln_kin_coeffs_line = parameters[0] + parameters[1]*temp_inv_line

    pylab.clf()
    pylab.title("Arrhenius plot: A [1/s] = %.5e    Ea [kJ/mol] = %.2f" % (A, Ea))
    pylab.xlabel("1/T [1/K]")
    pylab.ylabel("Unimolecular kinetic coefficient [1/s]")
    pylab.plot(temp_inv_line,ln_kin_coeffs_line,"r-",label="Fitted curve")
    pylab.plot(temp_inv,ln_kin_coeffs,"ro",label="Calculated values")
    pylab.legend(loc=0)
    pylab.savefig("%s.png" % prefix)


