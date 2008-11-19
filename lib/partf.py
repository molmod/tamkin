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
    "Contribution", "Electronic", "PHVA", "ExternalTranslation", "ExternalRotation",
    "PartFun", "compute_rate_coeff"
]



class IdealGasVolume(object):
    """Computes the volume of a molecule in an ideal gas as a function of the
    temperature, at a fixed refrence pressure.
    """

    def __init__(self, pressure=1*atm):
        self.pressure = pressure

    def __call__(self, temp):
        return boltzmann*temp/self.pressure

    def get_description(self):
        return "Molecular volume: ideal gas law, reference pressure [bar] = %.5f" % (self.pressure/bar)

    description = property(get_description)


class Contribution():
    # TODO: implement and test partial derivatives of the partition function
    # with respect to the temperature. Also implement functions like entropy
    # and heat capacity in a generic fashion based on the partial derivatives.
    def __init__(self, name):
        self.name = name

    def init_part_fun(self, molecule):
        pass

    def dump(self, f):
        print >> f, " Name: %s" % self.name

    def get_fixed_basis(self, coordinates):
        raise NotImplementedError

    def log_eval(self, temp):
        return 0.0


class Electronic(Contribution):
    # TODO: include this one automatically?
    def __init__(self, multiplicity=None):
        self.multiplicity = multiplicity
        Contribution.__init__(self, "electronic")

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


class PHVA(Contribution):
    def __init__(self, fixed_indices):
        self.fixed_indices = fixed_indices
        Contribution.__init__(self, "phva")

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


class ExternalTranslation(Contribution):
    def __init__(self, mol_volume=None):
        if mol_volume is None:
            self.mol_volume = IdealGasVolume()
        else:
            self.mol_volume = mol_volume
        Contribution.__init__(self, "translational")

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


class ExternalRotation(Contribution):
    def __init__(self, symmetry_number, im_threshold=1.0):
        self.symmetry_number = symmetry_number
        self.im_threshold = im_threshold
        Contribution.__init__(self, "rotational")

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


class Vibrations(Contribution):
    def __init__(self, other_contributions):
        self.other_contributions = other_contributions
        Contribution.__init__(self, "vibrational")

    def init_part_fun(self, molecule):
        self.hessian = molecule.hessian
        self.m3sqrt = numpy.sqrt(numpy.array([molecule.masses, molecule.masses, molecule.masses]).transpose().ravel())

        fixed_basis = []
        for other in self.other_contributions:
            for row in other.get_fixed_basis(molecule.coordinates):
                fixed_basis.append(self.to_weighted(row))

        fixed_basis = numpy.array(fixed_basis)
        U, V, Wt = numpy.linalg.svd(fixed_basis, full_matrices=True)

        # from now on, the basis vectors are orthogonal...
        num_fixed = sum(V > V.max()*1e-10)
        self.fixed_basis = Wt[:num_fixed]
        self.free_basis = Wt[num_fixed:]

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

    num_fixed = property(lambda self: self.fixed_basis.shape[0])
    num_free = property(lambda self: self.free_basis.shape[0])
    num_dim = property(lambda self: len(self.hessian))

    def to_free(self, mat):
        if len(mat.shape) == 2 and mat.shape[0] == self.num_dim and mat.shape[1] == self.num_dim:
            return numpy.dot(self.free_basis, numpy.dot(mat, self.free_basis.transpose()))
        elif len(mat.shape) == 2 and mat.shape[0] == self.num_dim and mat.shape[1] == 1:
            return numpy.dot(self.free_basis, mat.transpose()).transpose()
        elif len(mat.shape) == 1 and mat.shape[0] == self.num_dim:
            return numpy.dot(self.free_basis, mat)
        else:
            raise NotImplementedError

    def from_free(self, mat):
        if len(mat.shape) == 2 and mat.shape[0] == self.num_free and mat.shape[1] == self.num_free:
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

    def log_eval(self, temp):
        return self.log_eval_terms(temp).sum()

    def log_eval_terms(self, temp):
        # always with zero point energy correction...
        exp_arg = -2*numpy.pi*self.positive_freqs/boltzmann/temp
        return (exp_arg/2 - numpy.log(1-numpy.exp(exp_arg)))
        # This would be the version when the zero point energy corrections are
        # included in the energy difference when computing the reaction rate:
        #return -numpy.log(1-numpy.exp(exp_arg))


class PartFun(object):
    __reserved_names__ = set(["other_contributions", "vibrational"])

    def __init__(self, molecule, other_contributions):
        self.molecule = molecule
        # perform a sanity check on the names of the contributions:
        for other in other_contributions:
            if other.name in self.__reserved_names__:
                raise ValueError("A partition function contribution can not have the name '%s'" % other.name)
        # done testing, start initialization
        self.other_contributions = tuple(other_contributions)
        for other in self.other_contributions:
            other.init_part_fun(self.molecule)
        self.vibrational = Vibrations(self.other_contributions)
        self.vibrational.init_part_fun(self.molecule)
        # use convenient attribute names:
        for other in self.other_contributions:
            self.__dict__[other.name] = other

    def log_eval(self, temp):
        log_result = self.vibrational.log_eval(temp)
        for other in self.other_contributions:
            log_result += other.log_eval(temp)
        return log_result

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
    log_result = 0
    log_result += pf_trans.log_eval(temp)
    for pf_react in pfs_react:
        log_result -= pf_react.log_eval(temp)
    delta_E = pf_trans.molecule.energy - sum(pf_react.molecule.energy for pf_react in pfs_react)
    log_result += -delta_E/(boltzmann*temp)
    if len(pfs_react) > 1:
        if mol_volume is None:
            mol_volume = IdealGasVolume()
        log_result += numpy.log(mol_volume(temp))*(len(pfs_react)-1)
    return boltzmann*temp/(2*numpy.pi)*numpy.exp(log_result)


