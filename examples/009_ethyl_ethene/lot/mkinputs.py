#!/usr/bin/python


import sys, glob, string, os, numpy

from molmod.data.periodic import periodic
from molmod.io.xyz import XYZFile
from molmod.io.gaussian03.fchk import FCHKFile
from molmod.units import angstrom


def make_input(fn_com_template, dirname, lot_basis, mol, random=False):
    t = string.Template("".join(file(fn_com_template)))

    atom_lines = []
    for i in xrange(mol.size):
        if random:
            x, y, z = mol.coordinates[i]/angstrom + numpy.random.uniform(-0.1,0.1,3)
        else:
            x, y, z = mol.coordinates[i]/angstrom
        atom_lines.append(" %2s   % 10.5f   % 10.5f   % 10.5f" % (
            periodic[mol.numbers[i]].symbol, x, y, z,
        ))
    atom_str = "\n".join(atom_lines)

    fn_com = fn_com_template.replace("templates", dirname)
    dir_com = os.path.dirname(fn_com)
    if not os.path.isdir(dir_com):
        os.makedirs(dir_com)
    if not os.path.isfile(fn_com):
        print "Creating", fn_com
        f = file(fn_com, "w")
        f.write(t.substitute(lot_basis=lot_basis, atom_str=atom_str))
        f.close()
        return dir_com


def process_group(fn_com_opt, dirname, lot_basis, random=False):
    prefix = fn_com_opt.split("/")[1][:-4]
    fn_fchk_opt = "%s/%s_opt/gaussian.fchk" % (dirname, prefix)
    if os.path.isfile(fn_fchk_opt):
        # load the optimized geometry
        opt_mol = FCHKFile(fn_fchk_opt, field_labels=[]).molecule
        # write the other inputs
        for fn_com in glob.glob("templates/%s_*/gaussian.com" % prefix):
            if fn_com != fn_com_opt:
                dir_com = make_input(fn_com, dirname, lot_basis, opt_mol)
                if not (dir_com is None or os.path.isfile("%s/gaussian.in.fchk" % dir_com)):
                    os.symlink(
                        "../%s_opt/gaussian.fchk" % prefix,
                        "%s/gaussian.in.fchk" % dir_com,
                    )
    else:
        # load the initial geometry
        init_mol = XYZFile("templates/%s_init.xyz" % prefix).get_molecule()
        # write the optimization input
        make_input(fn_com_opt, dirname, lot_basis, init_mol, random)


lot_basis_to_dirname = lambda lot_basis: lot_basis.replace("/", "__").replace("(","").replace(")","")


usage = """USAGE: ./mkinputs.py [-s suffix] [-r] lot/basis

Prepares input files for a combination of level of theory (LOT) and basis set.
When executed for the first time, the optimization inputs are prepared. When
executed with completed optimizations, the corresponding inputs for frequency
and rotational scan computations are prepared. As soon as all these are
completed, one can use TAMkin to compute the kinetic parameters.

Each combination of LOT and basis is stored in its proper directory. The slash
is replaced by a double underscore, brackets are removed and the directory is
put in lower case, e.g. HF/STO-3G becomes hf__sto-3g and B3LYP/6-31G(D) becomes
b3lyp__6-31gd.
"""


def main():
    from optparse import OptionParser
    parser = OptionParser(usage)
    parser.add_option(
        "-r", "--random", default=False, action='store_true',
        help="Randomly distort the initial geometries with uniform "
             "displacements of [-0.1*angstrom,0.1*angstrom] in x, y, and z "
             "directions."
    )
    parser.add_option(
        "-s", "--suffix", default=None,
        help="Append __SUFFIX to the directory name"
    )
    options, args = parser.parse_args()

    if len(args) != 1:
        parser.error("Expecting exactly on argument")

    lot_basis = args[0].lower()
    dirname = lot_basis_to_dirname(lot_basis)
    if options.suffix is not None:
        dirname += '__%s' % options.suffix

    for fn_com_opt in glob.glob("templates/*_opt/gaussian.com"):
        process_group(fn_com_opt, dirname, lot_basis, options.random)


if __name__ == "__main__":
    main()

