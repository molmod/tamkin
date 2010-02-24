

from lot_basis import *

from molmod import angstrom, Molecule
from molmod.periodic import periodic
from molmod.io import FCHKFile

import os, numpy


__all__ = ["get_root", "states"]


def get_root(lot_label, basis_label, suffix):
    result = "%s__%s" % (lot_label, basis_label)
    if len(suffix) > 0:
        result = "%s__%s" % (result, suffix)
    return result


class State(object):
    def __init__(self, name, jobs):
        self.name = name
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
            mol = Molecule.from_file("init/%s.xyz" % state.name)
        else:
            # load the optimized geometry
            fn_fchk = "%s/%s__opt/gaussian.fchk" % (root, state.name)
            mol = FCHKFile(fn_fchk, field_labels=[]).molecule
        f = file("init/%s.fragments" % state.name)
        mol.charge_mult = f.readline().split()
        mol.tags = f.readline().split()
        f.close()
        return mol

    def write_input(self, state, root, lot_label, basis_label, suffix="", random=False):
        try:
            mol = self.get_mol(state, root)
        except IOError:
            return

        if self.name.startswith("cps"):
            i = int(self.name[-1])
            mult = int(mol.charge_mult[2*i+1])
        else:
            mult = int(mol.charge_mult[1])
        if lot_label.startswith("ro"):
            lot = get_lot(lot_label[2:], mult, restricted=True)
        else:
            lot = get_lot(lot_label, mult)
        basis = get_basis(basis_label)

        dirname = "%s/%s__%s" % (root, state.name, self.name)
        if not os.path.isdir(dirname):
            os.makedirs(dirname)

        fn_com = "%s/gaussian.com" % dirname
        if os.path.isfile(fn_com):
            return

        # copy an initial guess of the wavefunction
        destination = "%s/gaussian.in.fchk" % dirname
        if not (self.name == "opt" or self.name == "bsse" or self.name.startswith("cps")):
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
        if not (self.name == "opt" or self.name == "bsse" or self.name.startswith("cps")):
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
        if self.name == "bsse":
            print >> f, " ".join(mol.charge_mult)
        elif self.name.startswith("cps"):
            i = int(self.name[-1])
            print >> f, " ".join(mol.charge_mult[2*i:2*i+2])
        else:
            print >> f, " ".join(mol.charge_mult[:2])
        for i in xrange(mol.size):
            if random and self.name == "opt":
                x, y, z = mol.coordinates[i]/angstrom + numpy.random.uniform(-0.1,0.1,3)
            else:
                x, y, z = mol.coordinates[i]/angstrom
            symbol = periodic[mol.numbers[i]].symbol
            if self.name == "bsse":
                tag = mol.tags[i]
            else:
                tag = ""
            if self.name.startswith("cps"):
                if self.name.startswith("cps_full"):
                    if self.name[-1] != mol.tags[i] and self.name[-1] != "0":
                        symbol += "-Bq"
                else:
                    if self.name[-1] != mol.tags[i]:
                        continue
            print >> f, " %2s   % 10.5f   % 10.5f   % 10.5f  %s" % (
                symbol, x, y, z, tag
            )
        print >> f
        print >> f, self.post
        print >> f
        print >> f
        f.close()

        print dirname




states = [
    State("ethene", [
        G03Job("opt", "opt"),
        G03Job("freq", "freq(noraman)"),
        G03Job("sp", "sp"),
    ]),

    State("ethyl", [
        G03Job("opt", "opt"),
        G03Job("freq", "freq(noraman)"),
        G03Job("sp", "sp"),
        G03Job("scan_methyl", "opt(modredundant)", "4 1 2 6 S 72 5.0"),
    ]),

    State("ts_ad1_trans", [
        G03Job("opt", "opt(ts,calcfc,noeigentest)"),
        G03Job("freq", "freq(noraman)"),
        G03Job("bsse", "counterpoise=2"),
        G03Job("cps_full_0", "sp"),
        G03Job("cps_full_1", "sp"),
        G03Job("cps_full_2", "sp"),
        G03Job("cps_sole_1", "sp"),
        G03Job("cps_sole_2", "sp"),
        G03Job("sp", "sp"),
        G03Job("scan_forming_bond", "opt(ts,calcfc,modredundant,noeigentest)", ["1 2 3 4 S 72 5.0"]),
        G03Job("scan_methyl", "opt(ts,calcfc,modredundant,noeigentest)", ["5 1 2 8 S 72 5.0", "1 2 3 5 B"]),
    ]),

    State("ts_ad1_gauche", [
        G03Job("opt", "opt(ts,calcfc,noeigentest)"),
        G03Job("freq", "freq(noraman)"),
        G03Job("bsse", "counterpoise=2"),
        G03Job("cps_full_0", "sp"),
        G03Job("cps_full_1", "sp"),
        G03Job("cps_full_2", "sp"),
        G03Job("cps_sole_1", "sp"),
        G03Job("cps_sole_2", "sp"),
        G03Job("sp", "sp"),
        G03Job("scan_methyl", "opt(ts,calcfc,modredundant,noeigentest)", ["5 1 2 8 S 72 5.0", "1 2 3 5 B"]),
    ]),
]


