

from lot_basis import *

from molmod.units import angstrom
from molmod.data.periodic import periodic
from molmod.io.xyz import XYZFile
from molmod.io.gaussian03.fchk import FCHKFile

import os, numpy


__all__ = ["get_root", "states"]


def get_root(lot_label, basis_label, suffix):
    result = "%s__%s" % (lot_label, basis_label)
    if len(suffix) > 0:
        result = "%s__%s" % (result, suffix)
    return result


class State(object):
    def __init__(self, name, charge, mult, jobs):
        self.name = name
        self.charge = charge
        self.mult = mult
        self.jobs = jobs


class G03Job(object):
    def __init__(self, name, cmd, post=""):
        self.name = name
        self.cmd = cmd
        if isinstance(post, list):
            self.post = "\n".join(post)
        else:
            self.post = post

    def get_mol(self, state, root):
        if self.name == "opt":
            # load the initial geometry
            return XYZFile("init/%s.xyz" % state.name).get_molecule()
        else:
            # load the optimized geometry
            fn_fchk = "%s/%s__opt/gaussian.fchk" % (root, state.name)
            return FCHKFile(fn_fchk, field_labels=[]).molecule

    def write_input(self, state, root, lot_label, basis_label, suffix="", random=False):
        if lot_label.startswith("ro"):
            lot = get_lot(lot_label[2:], state.mult, restricted=True)
        else:
            lot = get_lot(lot_label, state.mult)
        basis = get_basis(basis_label)

        try:
            mol = self.get_mol(state, root)
        except IOError:
            return

        dirname = "%s/%s__%s" % (root, state.name, self.name)
        if not os.path.isdir(dirname):
            os.makedirs(dirname)

        fn_com = "%s/gaussian.com" % dirname
        if os.path.isfile(fn_com):
            return

        # copy an initial guess of the wavefunction
        destination = "%s/gaussian.in.fchk" % dirname
        if not self.name.startswith("opt"):
            if os.path.isfile(destination):
                os.remove(destination)
            os.symlink(
                "../%s__opt/gaussian.fchk" % state.name,
                destination,
            )

        # write the input file
        f = file(fn_com, "w")
        print >> f, "%chk=gaussian.chk"
        print >> f, "%nproc=1"
        print >> f, "%mem=1000MB"
        print >> f, "#p %s/%s %s %s maxdisk=5GB" % (lot.name, basis.name, self.cmd, lot.iop),
        if basis.diffuse:
            print >> f, "scf=tight",
        if not self.name.startswith("opt"):
            print >> f, "guess=read",
        if len(lot.extra_overlay) > 0:
            print >> f, "extraoverlay",
        print >> f
        print >> f
        if len(lot.extra_overlay) > 0:
            print >> f, lot.extra_overlay
            print >> f
        print >> f, dirname
        print >> f
        print >> f, "%i %i" % (state.charge, state.mult)
        for i in xrange(mol.size):
            if random and self.name == "opt":
                x, y, z = mol.coordinates[i]/angstrom + numpy.random.uniform(-0.1,0.1,3)
            else:
                x, y, z = mol.coordinates[i]/angstrom
            print >> f, " %2s   % 10.5f   % 10.5f   % 10.5f" % (
                periodic[mol.numbers[i]].symbol, x, y, z,
            )
        print >> f
        print >> f, self.post
        print >> f
        print >> f
        f.close()

        print dirname




states = [
    State("ethene", 0, 1, [
        G03Job("opt", "opt"),
        G03Job("freq", "freq(noraman)"),
        G03Job("sp", "sp"),
    ]),

    State("ethyl", 0, 2, [
        G03Job("opt", "opt"),
        G03Job("freq", "freq(noraman)"),
        G03Job("sp", "sp"),
        G03Job("scan_methyl", "opt(modredundant)", "4 1 2 6 S 72 5.0"),
    ]),

    State("ts_ad1_trans", 0, 2, [
        G03Job("opt", "opt(ts,calcfc,noeigentest)"),
        G03Job("freq", "freq(noraman)"),
        G03Job("sp", "sp"),
        G03Job("scan_forming_bond", "opt(ts,calcfc,modredundant,noeigentest)", ["1 2 3 4 S 72 5.0"]),
        G03Job("scan_methyl", "opt(ts,calcfc,modredundant,noeigentest)", ["5 1 2 8 S 72 5.0", "1 2 3 5 B"]),
    ]),

    State("ts_ad1_gauche", 0, 2, [
        G03Job("opt", "opt(ts,calcfc,noeigentest)"),
        G03Job("freq", "freq(noraman)"),
        G03Job("sp", "sp"),
        G03Job("scan_methyl", "opt(ts,calcfc,modredundant,noeigentest)", ["5 1 2 8 S 72 5.0", "1 2 3 5 B"]),
    ]),
]


