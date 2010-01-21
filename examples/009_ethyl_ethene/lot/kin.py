
import numpy, sys

from molmod.units import bar, liter, meter, mol, second, kjmol
from molmod.constants import boltzmann

__all__ = [
    "load_summary", "get_error_color",
    "temps", "invtemps", "experimental_k", "experimental_A", "experimental_Ea",
    "covariance_parameters", "experimental_parameters",
]


def load_summary(fn):
    f = file(fn)
    result = tuple(float(word) for word in f.next().split())
    f.close()
    return result[:-4], result[-4:-2], result[-2:]


def get_error_color(ratio):
    x = abs(ratio)
    if x < 1:
        rgb = (x,1.0,0.0)
    elif x < 2:
        rgb = (1.0,2.0-x,0.0)
    elif x < 3:
        rgb = ((4-x)/2,0.0,(x-2.0)/2)
    else:
        rgb = (0.5,0.0,0.5)
    result = "#%s" % "".join(hex(int(c*255))[2:].rjust(2,"0") for c in rgb)
    return result


temps = numpy.array([300, 400, 500, 600], float)
invtemps = 1/temps

centimeter = meter*0.01
experimental_k = 7.19e-15 * (temps/298)**2.44*numpy.exp(-22.45*kjmol/(boltzmann*temps))
experimental_k *= (centimeter**3/second)
experimental_k /= (meter**3/mol/second)
experimental_lnk = numpy.log(experimental_k)
print "%e" % (7.19e-15*(centimeter**3/second)/(meter**3/mol/second))
print experimental_k
print experimental_lnk

# fit experimental A and Ea and also sensitivity to errors

design_matrix = numpy.array([numpy.ones(4), -invtemps/boltzmann*kjmol]).transpose()
expected_values = experimental_lnk

A = numpy.dot(design_matrix.transpose(), design_matrix)
B = numpy.dot(design_matrix.transpose(), expected_values)

experimental_parameters = numpy.linalg.solve(A, B)
experimental_A = numpy.exp(experimental_parameters[0])
experimental_Ea = experimental_parameters[1]
covariance_lnk = numpy.ones((4,4),float)*numpy.log(10)**2 + numpy.identity(4,float)*numpy.log(2)**2
#covariance_lnk = numpy.identity(4,float)*numpy.log(10)**2
sensitivity = numpy.linalg.solve(A, design_matrix.transpose())
covariance_parameters = numpy.dot(numpy.dot(sensitivity, covariance_lnk), sensitivity.transpose())

print "%e" % experimental_A
print experimental_Ea
#sys.exit()
# done fitting

