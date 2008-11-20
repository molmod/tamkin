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


from molmod.constants import boltzmann, lightspeed
from molmod.units import atm, bar, amu, cm

import numpy


__all__ = [
    "IdealGasVolume", "FixedVolume", "StatFys", "StatFysTerms",
    "Contribution", "Electronic", "PHVA", "ExternalTranslation", "ExternalRotation",
    "PartFun", "compute_rate_coeff"
]



class IdealGasVolume(object):
    """Computes the volume of a molecule in an ideal gas as a function of the
    temperature, at a fixed refrence pressure.

    This law can be used to set up a constant pressure gas phase partition function
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
    """Computes the volume of a molecule in an ideal gas at a fixed reference
    temperature and at a fixed reference pressure.

    This law can be used to set up a constant volume gas phase partition function
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


class StatFys(object):
    def __init__(self, name):
        self.name = name

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

    def free_energy(self, temp, log_eval=None):
        """Computes the free energy per molecule (separate terms)"""
        return self.free_energy(temp, self.log_eval_terms)


class Contribution(object):
    def init_part_fun(self, molecule):
        pass

    def dump(self, f):
        print >> f, " Name: %s" % self.name

    def get_fixed_basis(self, coordinates):
        raise NotImplementedError


class Electronic(Contribution, StatFys):
    # TODO: include this one automatically?
    def __init__(self, multiplicity=None):
        self.multiplicity = multiplicity
        StatFys.__init__(self, "electronic")

    def init_part_fun(self, molecule):
        if self.multiplicity is None:
            self.multiplicity = molecule.multiplicity

    def dump(self, f):
        Contribution.dump(self, f)
        print >> f, "  Multiplicity: %i" % self.multiplicity

    def get_fixed_basis(self, coordinates):
        return []

    def log_eval(self, temp):
        return numpy.log(self.multiplicity)

    def log_deriv(self, temp):
        return 0.0

    def log_deriv2(self, temp):
        return 0.0


class PHVA(Contribution, StatFys):
    def __init__(self, fixed_indices):
        self.fixed_indices = fixed_indices
        StatFys.__init__(self, "phva")

    def dump(self, f):
        Contribution.dump(self, f)
        print >> f, "  Fixed atoms: %s" % (" ".join(str(i) for i in self.fixed_indices))

    def get_fixed_basis(self, coordinates):
        result = numpy.zeros((len(self.fixed_indices)*3, coordinates.size), float)
        for i, fixed_index in enumerate(self.fixed_indices):
            result[3*i  ,fixed_index*3  ] = 1
            result[3*i+1,fixed_index*3+1] = 1
            result[3*i+2,fixed_index*3+2] = 1
        return result

    def log_eval(self, temp):
        return 0.0

    def log_deriv(self, temp):
        return 0.0

    def log_deriv2(self, temp):
        return 0.0


class ExternalTranslation(Contribution, StatFys):
    def __init__(self, mol_volume=None):
        if mol_volume is None:
            self.mol_volume = FixedVolume()
        else:
            self.mol_volume = mol_volume
        StatFys.__init__(self, "translational")

    def init_part_fun(self, molecule):
        self.mass = molecule.mass

    def dump(self, f):
        Contribution.dump(self, f)
        print >> f, "  %s" % self.mol_volume.description

    def get_fixed_basis(self, coordinates):
        result = numpy.zeros((3, coordinates.size), float)
        result[0, 0::3] = 1
        result[1, 1::3] = 1
        result[2, 2::3] = 1
        return result

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



class ExternalRotation(Contribution, StatFys):
    def __init__(self, symmetry_number, im_threshold=1.0):
        self.symmetry_number = symmetry_number
        self.im_threshold = im_threshold
        StatFys.__init__(self, "rotational")

    def init_part_fun(self, molecule):
        self.com = molecule.com
        evals = numpy.linalg.eigvalsh(molecule.inertia_tensor)
        self.factor = numpy.sqrt(numpy.product([
            2*v*boltzmann for v in evals if v > self.im_threshold
        ])*numpy.pi)/self.symmetry_number
        self.count = (evals > self.im_threshold).sum()

    def dump(self, f):
        Contribution.dump(self, f)
        print >> f, "  Rotational symmetry number: %i" % self.symmetry_number
        print >> f, "  Inertia moment threshold [amu*bohr**2]: %e" % (self.im_threshold/amu)
        print >> f, "  Non-zero moments of inertia: %i" % self.count

    def get_fixed_basis(self, coordinates):
        result = numpy.zeros((3, coordinates.size), float)
        result[0, 1::3] =  coordinates[:,2] - self.com[2]
        result[0, 2::3] = -coordinates[:,1] + self.com[1]
        result[1, 2::3] =  coordinates[:,0] - self.com[0]
        result[1, 0::3] = -coordinates[:,2] + self.com[2]
        result[2, 0::3] =  coordinates[:,1] - self.com[1]
        result[2, 1::3] = -coordinates[:,0] + self.com[0]
        return result

    def log_eval(self, temp):
        return numpy.log(temp)*0.5*self.count + numpy.log(self.factor)

    def log_deriv(self, temp):
        return 0.5*self.count/temp

    def log_deriv2(self, temp):
        return -0.5*self.count/temp**2


class Vibrations(Contribution, StatFysTerms):
    def __init__(self, other_contributions, classical=False):
        self.other_contributions = other_contributions
        self.classical = classical
        StatFysTerms.__init__(self, "vibrational")

    def init_part_fun(self, molecule):
        self.hessian = molecule.hessian
        self.m3sqrt = numpy.sqrt(numpy.array([molecule.masses, molecule.masses, molecule.masses]).transpose().ravel())

        fixed_basis = []
        for other in self.other_contributions:
            for row in other.get_fixed_basis(molecule.coordinates):
                fixed_basis.append(self.to_weighted(row))

        if len(fixed_basis) > 0:
            fixed_basis = numpy.array(fixed_basis)
            U, V, Wt = numpy.linalg.svd(fixed_basis, full_matrices=True)

            # from now on, the basis vectors are orthogonal...
            num_fixed = sum(V > V.max()*1e-10)
            self.fixed_basis = Wt[:num_fixed]
            self.free_basis = Wt[num_fixed:]
        else:
            self.fixed_basis = []
            self.free_basis = None
        # solve the vibrational problem in the free basis
        free_hessian = self.to_free(self.to_weighted(molecule.hessian))
        evals, evecs = numpy.linalg.eigh(free_hessian)

        self.freqs = numpy.sqrt(abs(evals))/(2*numpy.pi)
        self.freqs *= (evals > 0)*2-1
        self.positive_freqs = self.freqs[self.freqs > 0]
        self.negative_freqs = self.freqs[self.freqs < 0]

        # convert the eigenmodes back to normal coordinates
        self.eigen_modes = []
        for evec in evecs.transpose():
            eigen_mode = self.from_weighted(self.from_free(evec))
            eigen_mode /= numpy.linalg.norm(eigen_mode)
            self.eigen_modes.append(eigen_mode)

    num_fixed = property(lambda self: len(self.fixed_basis))
    num_free = property(lambda self: len(self.hessian)-len(self.fixed_basis))
    num_dim = property(lambda self: len(self.hessian))

    def to_free(self, mat):
        if self.num_fixed == 0:
            return mat
        elif len(mat.shape) == 2 and mat.shape[0] == self.num_dim and mat.shape[1] == self.num_dim:
            return numpy.dot(self.free_basis, numpy.dot(mat, self.free_basis.transpose()))
        elif len(mat.shape) == 2 and mat.shape[0] == self.num_dim and mat.shape[1] == 1:
            return numpy.dot(self.free_basis, mat.transpose()).transpose()
        elif len(mat.shape) == 1 and mat.shape[0] == self.num_dim:
            return numpy.dot(self.free_basis, mat)
        else:
            raise NotImplementedError

    def from_free(self, mat):
        if self.num_fixed == 0:
            return mat
        elif len(mat.shape) == 2 and mat.shape[0] == self.num_free and mat.shape[1] == self.num_free:
            return numpy.dot(self.free_basis.transpose(), numpy.dot(mat, self.free_basis))
        elif len(mat.shape) == 2 and mat.shape[0] == self.num_free and mat.shape[1] == 1:
            return numpy.dot(self.free_basis.transpose(), mat.transpose()).transpose()
        elif len(mat.shape) == 1 and mat.shape[0] == self.num_free:
            return numpy.dot(self.free_basis.transpose(), mat)
        else:
            raise NotImplementedError

    def to_weighted(self, mat):
        if len(mat.shape) == 2 and mat.shape[0] == self.num_dim and mat.shape[1] == self.num_dim:
            return (mat/self.m3sqrt).transpose()/self.m3sqrt         # mass weighted hessian
        elif len(mat.shape) == 2 and mat.shape[0] == self.num_dim and mat.shape[1] == 1:
            return (mat.transpose()*self.m3sqrt).transpose()
        elif len(mat.shape) == 1 and mat.shape[0] == self.num_dim:
            return mat*self.m3sqrt
        else:
            raise NotImplementedError

    def from_weighted(self, mat):
        if len(mat.shape) == 2 and mat.shape[0] == self.num_dim and mat.shape[1] == self.num_dim:
            return (mat*self.m3sqrt).transpose()*self.m3sqrt         # mass weighted hessian
        elif len(mat.shape) == 2 and mat.shape[0] == self.num_dim and mat.shape[1] == 1:
            return (mat.transpose()/self.m3sqrt).transpose()
        elif len(mat.shape) == 1 and mat.shape[0] == self.num_dim:
            return mat/self.m3sqrt
        else:
            raise NotImplementedError

    def dump(self, f):
        Contribution.dump(self, f)
        def print_freqs(label, freqs):
            parts = ["  "]
            for counter, freq in enumerate(freqs):
                parts.append("% 8.1f" % (freq/lightspeed/(1/cm)))
                if counter % 6 == 5:
                    parts.append("\n  ")
            if len(parts) > 1:
                if parts[-1] == "\n  ":
                    parts.pop()
                parts.insert(0, "  %s Wavenumbers [1/cm]:\n" % label)
                print >> f, "".join(parts)

        print >> f, "  Total DOF: %i" % self.num_dim
        print >> f, "  Vibrational DOF: %i" % self.num_free
        print >> f, "  Other DOF: %i" % self.num_fixed
        print_freqs("Real", self.positive_freqs)
        print_freqs("Imaginary", self.negative_freqs)

    def log_eval_terms(self, temp):
        if self.classical:
            return numpy.log(0.5*boltzmann*temp/numpy.pi/self.positive_freqs)
        else:
            # The zero point correction is included in the partition function and
            # should not be taken into account when computing the reaction barrier.
            exp_arg = -2*numpy.pi*self.positive_freqs/boltzmann/temp
            return (exp_arg/2 - numpy.log(1-numpy.exp(exp_arg)))
            # This would be the version when the zero point energy corrections are
            # included in the energy difference when computing the reaction rate:
            #return -numpy.log(1-numpy.exp(exp_arg))

    def log_deriv_terms(self, temp):
        if self.classical:
            return numpy.ones(len(self.positive_freqs))/temp
        else:
            exp_arg = -2*numpy.pi*self.positive_freqs/boltzmann/temp
            exp_arg_deriv = -exp_arg/temp
            return exp_arg_deriv*(0.5-1/(1-numpy.exp(-exp_arg)))

    def log_deriv2_terms(self, temp):
        if self.classical:
            return -numpy.ones(len(self.positive_freqs))/temp**2
        else:
            exp_arg = -2*numpy.pi*self.positive_freqs/boltzmann/temp
            exp_arg_deriv = -exp_arg/temp
            exp_arg_deriv2 = -2*exp_arg_deriv/temp
            e = numpy.exp(-exp_arg)
            x = 1/(1-e)
            return exp_arg_deriv2*(0.5-x) + (exp_arg_deriv*x)**2*e


class PartFun(StatFys):
    __reserved_names__ = set(["other_contributions", "vibrational"])

    def __init__(self, molecule, other_contributions, classical_vib=False):
        self.molecule = molecule
        # perform a sanity check on the names of the contributions:
        for other in other_contributions:
            if other.name in self.__reserved_names__:
                raise ValueError("A partition function contribution can not have the name '%s'" % other.name)
        # done testing, start initialization
        self.other_contributions = tuple(other_contributions)
        for other in self.other_contributions:
            other.init_part_fun(self.molecule)
        self.vibrational = Vibrations(self.other_contributions, classical_vib)
        self.vibrational.init_part_fun(self.molecule)
        # use convenient attribute names:
        for other in self.other_contributions:
            self.__dict__[other.name] = other
        StatFys.__init__(self, "total")

    def log_eval(self, temp):
        result = self.vibrational.log_eval(temp)
        for other in self.other_contributions:
            result += other.log_eval(temp)
        return result

    def log_deriv(self, temp):
        result = self.vibrational.log_deriv(temp)
        for other in self.other_contributions:
            result += other.log_deriv(temp)
        return result

    def log_deriv2(self, temp):
        result = self.vibrational.log_deriv2(temp)
        for other in self.other_contributions:
            result += other.log_deriv2(temp)
        return result

    def entropy(self, temp):
        """Compute the total entropy"""
        return boltzmann + StatFys.entropy(self, temp)

    def free_energy(self, temp):
        """Computes the free energy"""
        return StatFys.free_energy(self, temp) + self.molecule.energy

    def dump(self, f):
        print >> f, "Chemical formula: %s" % self.molecule.chemical_formula
        print >> f, "Number of atoms: %s" % self.molecule.size
        print >> f, "Mass [amu]: %f" % (self.molecule.mass/amu)
        print >> f, "Moments of inertia [amu*bohr**2]: %f  %f %f" % tuple(
            numpy.linalg.eigvalsh(self.molecule.inertia_tensor)/amu
        )
        print >> f, "Contributions:"
        self.vibrational.dump(f)
        for other in self.other_contributions:
            other.dump(f)

    def write_to_file(self, filename):
        f = file(filename, 'w')
        self.dump(f)
        f.close()


def compute_rate_coeff(pfs_react, pf_trans, temp, mol_volume=None):
    delta_G = pf_trans.free_energy(temp) - sum(pf_react.free_energy(temp) for pf_react in pfs_react)
    log_result = -delta_G/(boltzmann*temp)
    if len(pfs_react) > 1:
        if mol_volume is None:
            mol_volume = FixedVolume()
        log_result += numpy.log(mol_volume(temp))*(len(pfs_react)-1)
    return boltzmann*temp/(2*numpy.pi)*numpy.exp(log_result)


