# TAMkin is a post-processing toolkit for thermochemistry and kinetics analysis.
# Copyright (C) 2008-2009 Toon Verstraelen <Toon.Verstraelen@UGent.be>,
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


from tamkin.data import Molecule

from molmod.io.gaussian03.fchk import FCHKFile
from molmod.io.xyz import XYZFile
from molmod.units import amu, calorie, avogadro, angstrom
from molmod.data.periodic import periodic

import numpy


__all__ = [
    "load_fixed_g03com", "load_molecule_g03fchk", "load_molecule_cp2k",
    "load_molecule_cpmd", "load_molecule_charmm",
    "load_chk", "dump_chk",
    "load_fixed_txt", "load_subs_txt", "load_envi_txt", "load_blocks_txt",
]


def load_fixed_g03com(filename):
    f = file(filename)
    for line in f:
        # iterate until we reach the title line
        line = line.strip()
        if not (len(line)==0 or line[0] == "#" or line[0] == "%"):
            break
    for line in f:
        # skip lines until charge/multiplicty is reached
        words = line.split()
        if len(words) == 2:
            try:
                int(words[0])
                int(words[1])
                break
            except ValueError:
                pass
    counter = 0
    fixed_atoms = []
    for line in f:
        # go trough the list of atoms and store the fixed ones.
        words = line.split()
        if len(words) == 0:
            break
        if words[1] == "-1":
            fixed_atoms.append(counter)
        counter += 1
    return fixed_atoms


def load_molecule_g03fchk(filename_freq,filename_ener=None): # if one file contains every information, give the name twice
    fchk_freq = FCHKFile(filename_freq, ignore_errors=True, field_labels=[
        "Cartesian Force Constants", "Real atomic weights", "Total Energy",
        "Multiplicity", "Cartesian Gradient"
    ])
    if filename_ener is None:
        fchk_ener = fchk_freq
    else:
        fchk_ener = FCHKFile(filename_ener, ignore_errors=True, field_labels=[
            "Total Energy"
        ])
    return Molecule(
        fchk_freq.molecule.numbers,
        fchk_freq.molecule.coordinates,
        fchk_freq.fields["Real atomic weights"]*amu,
        fchk_ener.fields["Total Energy"],
        fchk_freq.fields["Cartesian Gradient"],
        fchk_freq.get_hessian(),
        fchk_freq.fields["Multiplicity"],
        None, # gaussian is very poor at computing the rotational symmetry number
        False,
    )


def load_molecule_cp2k(fn_xyz, fn_sp, fn_freq, multiplicity=1, is_periodic=True):
    molecule = XYZFile(fn_xyz).get_molecule()
    masses = numpy.array([periodic[number].mass for number in molecule.numbers])

    # go trhough the single point file: energy and gradient
    energy = None
    gradient = None
    f = file(fn_sp)
    while True:
        line = f.readline()
        if line == "":
            break
        if line.startswith(" ENERGY|"):
            energy = float(line[60:])
        elif line.startswith(" FORCES|"):
            tmp = []
            while True:
                line = f.readline()
                if line == "\n":
                    break
                if line == "":
                    raise IOError("End of file while reading gradient (forces).")
                words = line.split()
                tmp.append([float(words[1]), float(words[2]), float(words[3])])
            gradient = -numpy.array(tmp) # force to gradient
            break
    if energy is None or gradient is None:
        raise IOError("Could not read energy and/or gradient (forces) from single point file.")
    f.close()

    # go trhough the freq file: hessian
    hessian = None
    f = file(fn_freq)
    while True:
        line = f.readline()
        if line.startswith(" VIB| Hessian in cartesian coordinates"):
            block_len = molecule.size*3
            tmp = numpy.zeros((block_len,block_len), float)
            i2 = 0
            while i2 < block_len:
                num_cols = min(5, block_len-i2)
                f.readline() # skip two lines
                f.readline()
                for j in xrange(block_len):
                    line = f.readline()
                    if line == "":
                        raise IOError("End of file while reading hessian.")
                    words = line.split()
                    for i1 in xrange(num_cols):
                        tmp[i2+i1,j] = float(words[i1+2])
                i2 += num_cols
            hessian = tmp
            break
    f.close()
    if hessian is None:
        raise IOError("Could not read hessian from freq file.")

    hessian = 0.5*(hessian+hessian.transpose())
    # cp2k prints a transformed hessian, here we convert it back to the normal
    # hessian in atomic units.
    conv = 1e-3*numpy.array([masses, masses, masses]).transpose().ravel()**0.5
    hessian *= conv
    hessian *= conv.reshape((-1,1))

    return Molecule(
        molecule.numbers, molecule.coordinates, masses, energy, gradient,
        hessian, multiplicity, None, is_periodic
    )


def load_molecule_cpmd(fn_out, fn_geometry, fn_hessian, multiplicity=1, is_periodic=True):
    # go through the output file: grep the total energy
    energy = None
    f = file(fn_out)
    while True:
        line = f.readline()
        if line == "":
            raise IOError("Could not find final results in %s. Is the output file truncated?" % fn_out)
        if line == " *                        FINAL RESULTS                         *\n":
            break
    while True:
        line = f.readline()
        if line == "":
            raise IOError("Could not find total energy in %s. Is the output file truncated?" % fn_out)
        if line.startswith(" (K+E1+L+N+X)           TOTAL ENERGY ="):
            words= (line.strip()).split()
            energy = float(words[4])
            break
    f.close()

    f = file(fn_geometry)
    num_atoms = int(f.readline())
    numbers = numpy.zeros(num_atoms, int)
    coordinates = numpy.zeros((num_atoms,3), float)
    gradient = numpy.zeros((num_atoms,3), float)

    f.readline()
    i = 0
    while True:
        line = f.readline()
        if line == "":
            break # end of file
        words = (line.strip()).split()
        if len(words) == 7:
            numbers[i] = periodic[words[0]].number
            coordinates[i][0] = float(words[1])*angstrom
            coordinates[i][1] = float(words[2])*angstrom
            coordinates[i][2] = float(words[3])*angstrom
            gradient[i][1] = float(words[4])
            gradient[i][1] = float(words[5])
            gradient[i][2] = float(words[6])
            i += 1
        else:
            raise IOError("Expecting seven words at each atom line in %s." % fn_geometry)
    if i != num_atoms:
        raise IOError("The number of atoms is incorrect in %s." % fn_geometry)
    f.close()

    # go trhough the freq file: hessian
    f = file(fn_hessian)

    line = f.readline()
    if not line.startswith(" &CART"):
        raise IOError("File %s does not start with &CART." % fn_hessian)
    masses = numpy.zeros(num_atoms, float)
    for i in xrange(num_atoms):
        line = f.readline()
        words = line.split()
        masses[i] = float(words[4])*amu
    f.readline() # &END

    line = f.readline()
    if not line.startswith(" &FCON"):
        raise IOError("File %s does not contain section &FCON." % fn_hessian)
    num_cart = num_atoms*3
    hessian = numpy.zeros((num_cart, num_cart), float)
    for i in xrange(num_cart):
        line = f.readline()
        words = line.split()
        for j in xrange(num_cart):
            hessian[i,j] = float(words[j])

    f.close()

    return Molecule(
        numbers, coordinates, masses, energy, gradient, hessian, multiplicity,
        None, is_periodic
    )


def load_molecule_charmm(charmmfile_cor, charmmfile_hess,
                  is_periodic = False):
    # Units: CHARMM gradient in kcal/mol/angstrom, TAMkin internally all in atomic units

    # Read from first CHARMM file
    # format:  nb of atoms = N
    #          energy
    #          gradient (N lines, 3 elements on each line)
    #          Hessian  (upper triangular form, 1 element on each line)
    #          coordinates (N lines, 3 elements on each line)
    f = file(charmmfile_hess)

    N = int(f.readline().split()[-1])   # nb of atoms
    assert N > 0   # nb of atoms should be > 0

    energy = float(f.readline().split()[-1]) * 1000*calorie/avogadro

    gradient = numpy.zeros((N,3),float)
    for i,line in enumerate(f):
        words = line.split()
        gradient[i,:] = [float(word) for word in words]
        if i == (N-1):
            break
    gradient *= 1000*calorie/avogadro/angstrom
    hessian = numpy.zeros((3*N,3*N),float)
    row = 0
    col = 0
    for line in f:
        element = float(line.split()[-1])
        hessian[row,col]=element
        hessian[col,row]=element
        col += 1
        if col>=3*N:        #if this new col doesn't exist
            row += 1        #go to next row
            col = row       #to diagonal element
            if row >= 3*N:  #if this new row doesn't exist
               break
    hessian = hessian * 1000*calorie/avogadro /angstrom**2

    positions = numpy.zeros((N,3),float)
    for i,line in enumerate(f):
        words = line.split()
        positions[i,:] = [float(word)*angstrom for word in words]
        if i == (N-1):
            break
    f.close()

    # Read from second CHARMM file
    # format:  header lines, which start with *
    #          N lines with   - mass in last column
    #                         - atomic type in 4th column
    f = file(charmmfile_cor)
    masses = numpy.zeros(N,float)
    symbols  = []
    for line in f:
        if not line.startswith("*"): # skip header lines
            break
    for i,line in enumerate(f):
        words = line.split()
        masses[i] = float(words[-1])*amu   # mass
        symbols.append( words[3] )         # symbol
        if i == (N-1):
            break
    f.close()

    # get corresponding atomic numbers
    mass_table = numpy.zeros(len(periodic))
    for i in xrange(1, len(mass_table)):
        m1 = periodic[i].mass
        if m1 is None:
            m1 = 200000.0
        m2 = periodic[i+1].mass
        if m2 is None:
            m2 = 200000.0
        mass_table[i] = 0.5*(m1+m2)
    atomicnumbers = numpy.zeros(N, int)
    for i,mass in enumerate(masses):
        atomicnumbers[i] = mass_table.searchsorted(mass)

    return Molecule(
        atomicnumbers, positions, masses, energy, gradient,
        hessian, 1, None, is_periodic
    )


def load_chk(filename):
    f = file(filename)
    result = {}
    while True:
        line = f.readline()
        if line == "":
            break
        if len(line) < 54:
            raise IOError("Header lines must be at least 54 characters long.")
        key = line[:40].strip()
        kind = line[47:52].strip()
        value = line[53:-1] # discard newline
        if kind == 'str':
            result[key] = value
        elif kind == 'int':
            result[key] = int(value)
        elif kind == 'flt':
            result[key] = float(value)
        elif kind[3:5] == 'ar':
            if kind[:3] == 'int':
                dtype = int
            elif kind[:3] == 'flt':
                dtype = float
            else:
                raise IOError("Unsupported kind: %s" % kind)
            shape = tuple(int(i) for i in value.split(","))
            array = numpy.zeros(shape, dtype)
            work = array.ravel()
            counter = 0
            while True:
                short = f.readline().split()
                for s in short:
                    work[counter] = dtype(s)
                    counter += 1
                    if counter == array.size:
                        break
                if counter == array.size:
                    break
            result[key] = array
        elif kind == 'none':
            result[key] = None
        else:
            raise IOError("Unsupported kind: %s" % kind)
    f.close()
    return result


def dump_chk(filename, data):
    f = file(filename, "w")
    for key, value in sorted(data.iteritems()):
        if not isinstance(key, str):
            raise TypeError("The keys must be strings.")
        if len(key) > 40:
            raise ValueError("Key strings can not be longer than 40 characters.")
        if isinstance(value, str):
            if len(value) > 256:
                raise TypeError("Only small strings are supported (256 chars).")
            if "\n" in value:
                raise ValueError("The string can not contain new lines.")
            print >> f, "%40s  kind=str   %s" % (key.ljust(40), value)
        elif isinstance(value, int):
            print >> f, "%40s  kind=int   %i" % (key.ljust(40), value)
        elif isinstance(value, float):
            print >> f, "%40s  kind=flt   %22.15e" % (key.ljust(40), value)
        elif isinstance(value, numpy.ndarray):
            if value.dtype.fields is not None:
                raise TypeError("Arrays with fields are not supported.")
            shape_str = ",".join(str(i) for i in value.shape)
            if issubclass(value.dtype.type, int):
                print >> f, "%40s  kind=intar %s" % (key.ljust(40), shape_str)
                format_str = "%22i"
            elif issubclass(value.dtype.type, float):
                print >> f, "%40s  kind=fltar %s" % (key.ljust(40), shape_str)
                format_str = "%22.15e"
            else:
                raise TypeError("Numpy array dtype %s not supported." % value.dtype)
            short_len = 4
            short = []
            for x in value.ravel():
                short.append(x)
                if len(short) == 4:
                    print >> f, " ".join(format_str  % s for s in short)
                    short = []
            if len(short) > 0:
                print >> f, " ".join(format_str  % s for s in short)
        elif value is None:
            print >> f, "%40s  kind=none   None" % key.ljust(40)
        else:
            raise TypeError("Type %s not supported." % type(value))
    f.close()


def load_fixed_txt(filename,shift=-1):
    """Read the fixed atoms into the list fixed_atoms.
    Empty lines are skipped.

    Arguments:
    filename  --  file that contains the fixed atoms:
                  1 atom on every line, empty lines are skipped
    shift  --  default on -1, because numbering in Python starts with -1
    """
    fixed_atoms = []
    f = file(filename)
    for line in f:
        if line != "\n":    # skip empty lines
            fixed_atoms.append(int(line)+shift)
    f.close()
    # TODO
    # check that every atom appears once
    return fixed_atoms


def load_subs_txt(filename,shift=-1):
    """Read the subsystem atoms into a list (for VSA treatment).
    Empty lines are skipped.

    Arguments:
    filename  --  file that contains the subsystem atoms:
                  1 atom on every line, empty lines are skipped
    shift  --  default on -1, because numbering in Python starts with -1
    """
    return load_fixed_txt(filename,shift=shift)


def load_envi_txt(filename,shift=-1):
    """Read the environment atoms into a list (for VSA treatment).
    Empty lines are skipped.

    Arguments:
    filename  --  file that contains the environment atoms:
                  1 atom on every line, empty lines are skipped
    shift  --  default on -1, because numbering in Python starts with -1
    """
    return load_fixed_txt(filename,shift=shift)


def load_blocks_txt(filename,shift=-1):
    pass






