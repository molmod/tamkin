#!/usr/bin/python


import glob, os
from optparse import OptionParser

from states import *

usage = "TODO"

def main():
    parser = OptionParser(usage)
    parser.add_option(
        "-r", "--random", default=False, action='store_true',
        help="Randomly distort the initial geometries with uniform "
             "displacements of [-0.1*angstrom,0.1*angstrom] in x, y, and z "
             "directions."
    )
    parser.add_option(
        "-s", "--suffix", default="",
        help="Append __SUFFIX to the directory name"
    )
    options, args = parser.parse_args()

    if len(args) == 0:
        parser.error("Expecting at least one argument")

    if args[0] == "create":
        if len(args) != 3:
            parser.error("create needs two extra arguments")
        lot_label, basis_label = args[1:]
        root = get_root(lot_label, basis_label, options.suffix)
        for state in states:
            for job in state.jobs:
                if job.name == "sp":
                    continue
                job.write_input(state, root, lot_label, basis_label, options.suffix, options.random)
    elif args[0] == "clone":
        if len(args) != 5:
            parser.error("create needs four extra arguments")
        lot0_label, basis0_label, lot1_label, basis1_label = args[1:]
        root0 = get_root(lot0_label, basis0_label, options.suffix)
        root1 = "GEO__%s__ENERGY__%s" % (
            get_root(lot0_label, basis0_label, ""),
            get_root(lot1_label, basis1_label, options.suffix)
        )
        if not os.path.isdir(root1):
            os.mkdir(root1)
        for dirname in glob.glob("%s/*" % root0):
            if os.path.isdir(dirname) and not (dirname.endswith("sp") or dirname.endswith("bsse")):
                jobdir = os.path.basename(dirname)
                source = "../%s/%s" % (root0, jobdir)
                destination = "%s/%s" % (root1, jobdir)
                if not os.path.exists(destination):
                    os.symlink(source, destination)
        for state in states:
            for job in state.jobs:
                if job.name == "bsse" or job.name == "sp":
                    job.write_input(state, root1, lot1_label, basis1_label, options.suffix, options.random)
    else:
        parser.error("Unknown command: %s" % args[0])



if __name__ == "__main__":
    main()

