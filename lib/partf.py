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


from molmod.constants import boltzmann

import numpy


__all__ = ["Contribution", "PHVA", "PartFun", "calc_kin_coeff"]



class Contribution():
    def get_fixed_basis(self, coordinates):
        raise NotImplementedError

    def init_part_fun(self, data):
        pass

    def log_part_fun(self, temp):
        return 0.0


class PHVA(Contribution):
    def __init__(self, fixed_indices):
        self.fixed_indices = fixed_indices

    def get_fixed_basis(self, coordinates):
        result = numpy.zeros((len(self.fixed_indices)*3, coordinates.size), float)
        for i, fixed_index in enumerate(self.fixed_indices):
            result[3*i  ,fixed_index*3  ] = 1
            result[3*i+1,fixed_index*3+1] = 1
            result[3*i+2,fixed_index*3+2] = 1
        return result


class PartFun(object):
    def __init__(self, data, contributions):
        self.__dict__.update(data)
        self.contributions = contributions # The non-vibrational contributions

        self._m3sqrt = numpy.sqrt(numpy.array([self.masses, self.masses, self.masses]).transpose().ravel())
        self._init_vibrational()

    def _init_vibrational(self):
        fixed_basis = []
        for contribution in self.contributions:
            for row in contribution.get_fixed_basis(self.coordinates):
                fixed_basis.append(self.to_weighted(row))

        fixed_basis = numpy.array(fixed_basis)
        U, V, Wt = numpy.linalg.svd(fixed_basis, full_matrices=True)

        # from now on, the basis vectors are orthogonal...
        num_fixed = sum(V > V.max()*1e-10)
        self.fixed_basis = Wt[:num_fixed]
        self.free_basis = Wt[num_fixed:]

        # solve the vibrational problem in the free basis
        free_hessian = self.to_free(self.to_weighted(self.hessian))
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
            return (mat/self._m3sqrt).transpose()/self._m3sqrt         # mass weighted hessian
        elif len(mat.shape) == 2 and mat.shape[0] == self.num_dim and mat.shape[1] == 1:
            return (mat.transpose()*self._m3sqrt).transpose()
        elif len(mat.shape) == 1 and mat.shape[0] == self.num_dim:
            return mat*self._m3sqrt
        else:
            raise NotImplementedError

    def from_weighted(self, mat):
        if len(mat.shape) == 2 and mat.shape[0] == self.num_dim and mat.shape[1] == self.num_dim:
            return (mat*self._m3sqrt).transpose()*self._m3sqrt         # mass weighted hessian
        elif len(mat.shape) == 2 and mat.shape[0] == self.num_dim and mat.shape[1] == 1:
            return (mat.transpose()/self._m3sqrt).transpose()
        elif len(mat.shape) == 1 and mat.shape[0] == self.num_dim:
            return mat/self._m3sqrt
        else:
            raise NotImplementedError

    def log_eval(self, temp):
        # vibrational part
        exp_arg = -2*numpy.pi*self.positive_freqs/(boltzmann*temp)
        log_result = (exp_arg/2 - numpy.log(1-numpy.exp(exp_arg))).sum()
        # other parts:
        for contribution in self.contributions:
            log_result += contribution.log_part_fun(temp)
        return log_result


def calc_kin_coeff(pfs_react, pf_trans, temp):
    delta_E = pf_trans.energy - sum(pf_react.energy for pf_react in pfs_react)
    log_result = -delta_E/(boltzmann*temp)
    log_result += pf_trans.log_eval(temp)
    for pf_react in pfs_react:
        log_result -= pf_react.log_eval(temp)
    return boltzmann*temp/(2*numpy.pi)*numpy.exp(log_result)


def construct_q_rot(detI,T):
    prefactor = 2*numpy.sqrt(2*math.pi) / 1   # default symmetrynb is 1
    q_rot       = prefactor *float(T)*(1.5) * float(detI)**(0.5) # TODO controle op correctheid en eenheden
    return q_rot

def construct_q_vib(freqs,T):
    f = freqs
    T = T
    for i in freqs:
        q_vib *= numpy.exp(-float(f[i])/(2.0*float(T))) / ( 1.0 - numpy.exp(-float(f[i])/(T)) )
    return q_vib

def construct_q_trans(T,totalmass):
    totm = totalmass
    pressure  = 101.325 * 1000 * (length_au_to_m)**3 / energy_au_to_J #TODO Ans eenheden vervangen door die van Toon
    volume    = T / pressure
    prefactor = (1/(2.0*math.pi))**(1.5)  * volume
    q_trans       = prefactor * (totalmass *float(T) )**(1.5)
    return q_trans

def q_total(detI,freqs,T,totalmass,freqmin):
    q_elec = 1
    q_vib = construct_q_vib(freqs,T)
    q_trans = construct_q_trans(p,T,totalmass)
    q_rot = construct_q_rot(detI,T)
    return q_elec*q_vib*q_rot*q_trans, freqmin



