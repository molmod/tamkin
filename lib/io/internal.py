# TAMkin is a post-processing toolkit for thermochemistry and kinetics analysis.
# Copyright (C) 2008-2010 Toon Verstraelen <Toon.Verstraelen@UGent.be>,
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


import numpy


__all__ = [
    "load_chk", "dump_chk",
    "load_indices", "dump_indices",
]



def load_chk(filename):
    """Load a TAMKin checkpoint file

       Argument:
        | filename  --  the file to load from

       The return value is a dictionary whose keys are field labels and the
       values can be None, string, integer, float, boolean or an array of
       strings, integers, booleans or floats.

       The file format is similar to the Gaussian fchk format, but has the extra
       feature that the shapes of the arrays are also stored.
    """
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
        elif kind == 'bln':
            result[key] = bool(value)
        elif kind == 'flt':
            result[key] = float(value)
        elif kind[3:5] == 'ar':
            if kind[:3] == 'str':
                dtype = str
            elif kind[:3] == 'int':
                dtype = int
            elif kind[:3] == 'bln':
                dtype = bool
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
                    if dtype is bool and s.lower() in ["True", "1", "Yes"]:
                        work[counter] = True
                    else:
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
    """Dump a TAMKin checkpoint file

       Argument:
        | filename  --  the file to write to
        | data  -- a dictionary whose keys are field labels and the values can
                   be None, string, integer, float, boolean, an array/list of
                   strings, integers, floats or booleans.

       The file format is similar to the Gaussian fchk format, but has the extra
       feature that the shapes of the arrays are also stored.
    """
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
        elif isinstance(value, bool):
            print >> f, "%40s  kind=bln   %s" % (key.ljust(40), value)
        elif isinstance(value, float):
            print >> f, "%40s  kind=flt   %22.15e" % (key.ljust(40), value)
        elif isinstance(value, numpy.ndarray) or isinstance(value, list):
            if isinstance(value, list):
                value = numpy.array(value)
            if value.dtype.fields is not None:
                raise TypeError("Arrays with fields are not supported.")
            shape_str = ",".join(str(i) for i in value.shape)
            if issubclass(value.dtype.type, str):
                for cell in value.flat:
                    if len(cell) >= 22:
                        raise ValueError("In case of string arrays, a string may contain at most 21 characters.")
                    if " " in cell or "\n" in cell:
                        raise ValueError("In case of string arrays, a string may not contain spaces or new lines.")
                print >> f, "%40s  kind=strar %s" % (key.ljust(40), shape_str)
                format_str = "%22s"
            elif issubclass(value.dtype.type, int):
                print >> f, "%40s  kind=intar %s" % (key.ljust(40), shape_str)
                format_str = "%22i"
            elif issubclass(value.dtype.type, numpy.bool_):
                print >> f, "%40s  kind=blnar %s" % (key.ljust(40), shape_str)
                format_str = "%22s"
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


def load_indices(filename, shift=-1, groups=False):
    """Load atom indexes from file

       Individual indexes or separated by white space or by one new-line
       character. The atom indexes can be separated by one or more empty lines
       to define groups of atoms that belong toghether in one block.

       Arguments:
         | filename  --  The file to load the indexes from

       Optional arguments:
         | shift  --  A constant shift applied to all atom indexes to convert
                      between numbers starting from zero and numbers starting from
                      one.
         | groups  --  When True, the function always returns a list of lists,
                       even when only one group of indexes is found in the file.
                       Otherwise only a single list of indexes is returned, even
                       when multiple groups of indexes are encountered.

    """
    blocks = []
    block  = []
    f = file(filename)
    for line in f:
        words = line.split()
        if len(words) == 0 and len(block) > 0:
            blocks.append(block)
            # start new block
            block = []
        else:
            for word in words:
                block.append(int(word)+shift)
    if len(block) > 0:
        # add last block to blocks list
        blocks.append(block)
    f.close()
    if groups:
        return blocks
    else:
        return sum(blocks, [])


def dump_indices(filename, indices, shift=1, compact=True):
    """Dump atom indexes to file

       Arguments:
         | filename  --  the file to dump the atom indexes to

       Optional arguments:
         | indices  --  a list of atom indices or a list of lists of atom indices
                        (the latter is used to define blocks of atoms)
         | shift  --  A constant shift applied to all atom indexes to convert
                      between numbers starting from zero and numbers starting
                      from one.
         | compact  --  When True, the numbers belonging to the same block are
                        put on one line, otherwise they are separated by a
                        newline character. Different blocks are separated by an
                        empty line. [default=True]
    """

    if len(indices) > 0 and not hasattr(indices[0], "__len__"):
        indices = [indices]

    separator = {True: " ", False: "\n"}[compact]

    f = file(filename, "w")
    for l in indices:
        group_str = separator.join(str(index+shift) for index in l)
        print >> f, group_str
        print >> f
    f.close()
