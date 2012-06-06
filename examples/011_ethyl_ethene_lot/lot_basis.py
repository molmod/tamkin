# -*- coding: utf-8 -*-
# TAMkin is a post-processing toolkit for normal mode analysis, thermochemistry
# and reaction kinetics.
# Copyright (C) 2008-2012 Toon Verstraelen <Toon.Verstraelen@UGent.be>, An Ghysels
# <An.Ghysels@UGent.be> and Matthias Vandichel <Matthias.Vandichel@UGent.be>
# Center for Molecular Modeling (CMM), Ghent University, Ghent, Belgium; all
# rights reserved unless otherwise stated.
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
# parts of this program are required to cite the following article:
#
# "TAMkin: A Versatile Package for Vibrational Analysis and Chemical Kinetics",
# An Ghysels, Toon Verstraelen, Karen Hemelsoet, Michel Waroquier and Veronique
# Van Speybroeck, Journal of Chemical Information and Modeling, 2010, 50,
# 1736-1750W
# http://dx.doi.org/10.1021/ci100099g
#
# TAMkin is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
#--


__all__ = ["get_lot", "get_basis"]


class G03LevelOfTheory(object):
    def __init__(self, spin, label, color, linestyle, name, iop="", extra_overlay=""):
        self.spin = spin
        self.label = label
        self.color = color
        self.linestyle = linestyle
        self.name = name
        self.iop = iop
        self.extra_overlay = extra_overlay
        if spin =="ROS":
            self.key = "ro%s" % self.label
        else:
            self.key = self.label


colors = [
    "#4040FF",
    "#000080",
    "#80A0FF",

    "#EE00FF",

    "#77FF77",
    "#008000",

    "#0090A0",

    "deeppink",
    #"magenta",
    "red",
    "darkorange",
    "gold",

    "pink",
]

lots_list = [
    # hartree fock
    G03LevelOfTheory("DEF", "hf",        colors[0],  "-",  "HF"),
    G03LevelOfTheory("ROS", "hf",        colors[0],  "--", "ROHF"),
    # moller plesset
    G03LevelOfTheory("DEF", "mp2",       colors[1],  "-",  "MP2"),
    G03LevelOfTheory("ROS", "mp2",       colors[1],  "--", "ROMP2"),
    G03LevelOfTheory("DEF", "mp3",       colors[1],  ":",  "MP3"),
    G03LevelOfTheory("DEF", "mp4",       colors[1],  "-.", "MP4"),
    # coupled cluster
    G03LevelOfTheory("DEF", "ccd",       colors[2],  "-",  "CCD"),
    G03LevelOfTheory("DEF", "ccsd",      colors[2],  "--", "CCSD"),
    #G03LevelOfTheory("DEF", "ccsdt",     colors[2],  ":",  "CCSD-T"),

    # lda
    G03LevelOfTheory("DEF", "hfs",       colors[3],  "-",  "HFS"),
    # gga (generalizef gradient approximation)
    G03LevelOfTheory("DEF", "hcth93",    colors[4],  "-",  "HCTH93"),
    G03LevelOfTheory("DEF", "hcth147",   colors[4],  "--", "HCTH147"),
    G03LevelOfTheory("DEF", "hcth407",   colors[4],  ":",  "HCTH407"),
    G03LevelOfTheory("DEF", "blyp",      colors[4],  "-.", "BLYP"),
    G03LevelOfTheory("DEF", "olyp",      colors[5],  "-",  "OLYP"),
    G03LevelOfTheory("DEF", "pw91pw91",  colors[5],  "--", "PW91PW91"),
    G03LevelOfTheory("DEF", "hfb",       colors[5],  ":",  "HFB"),
    G03LevelOfTheory("DEF", "pbepbe",    colors[5],  "-.", "PBEPBE"),
    G03LevelOfTheory("DEF", "tpsstpss",  colors[6],  "-",  "TPSSTPSS"),
    G03LevelOfTheory("DEF", "vsxc",      colors[6],  "--", "VSXC"),
    # hybrid-gga (with fock exchange)
    G03LevelOfTheory("DEF", "b1lyp",     colors[7],  "--", "B1LYP"),
    G03LevelOfTheory("DEF", "b3pw91",    colors[7],  ":",  "B3PW91"),
    G03LevelOfTheory("DEF", "b3p86",     colors[7],  "-.", "B3P86"),
    G03LevelOfTheory("DEF", "b3lyp",     colors[8],  "-",  "B3LYP"),
    G03LevelOfTheory("DEF", "o3lyp",     colors[8],  "--", "O3LYP"),
    G03LevelOfTheory("DEF", "bhandh",    colors[8],  ":",  "BHANDH"),
    G03LevelOfTheory("DEF", "bhandhlyp", colors[8],  "-.", "BHANDHLYP"),
    G03LevelOfTheory("DEF", "b971",      colors[9],  ":",  "B971"),
    G03LevelOfTheory("DEF", "b972",      colors[9],  "-.", "B972"),
    G03LevelOfTheory("DEF", "pbe1pbe",   colors[10], "-",  "PBE1PBE"),
    G03LevelOfTheory("DEF", "mpwpw91",   colors[10], "--", "MPWPW91"),
    G03LevelOfTheory("DEF", "mpw1pw91",  colors[10], ":",  "MPW1PW91"),
    G03LevelOfTheory("DEF", "mpw1k",     colors[11], "-",  "MPWPW91",  "iop(3/76=0572004280)"),
    # meta-hybrid ()
    G03LevelOfTheory("DEF", "b98",       colors[9],  "--", "B98"),
    G03LevelOfTheory("DEF", "b1b95",     colors[7],  "-",  "B1B95"),
    G03LevelOfTheory("DEF", "mpwb95",    colors[10], "-.", "MPWB95"),
    G03LevelOfTheory("DEF", "bmk",       colors[9],  "-",  "BMK"),
    G03LevelOfTheory("DEF", "mpwb1k",    colors[11], "--", "MPWB95",   "iop(3/76=0560004400)"),
    # double-hybrid (also with mp2 based corrections)
    G03LevelOfTheory("DEF", "b2plyp",    colors[11], ":",  "BLYP",     "iop(3/76=0470005300,3/78=0730007300)", "8/10=90/1;9/16=-3/6;"),
    G03LevelOfTheory("DEF", "mpw2plyp",  colors[11], "-.", "MPWLYP",   "iop(3/76=0450005500,3/78=0750007500)", "8/10=90/1;9/16=-3/6;"),
    G03LevelOfTheory("ROS", "b2plyp",    colors[11], ":",  "ROBLYP",   "iop(3/76=0470005300,3/78=0730007300)", "8/5=-1,10=90/1;9/16=-3/6;"),
    G03LevelOfTheory("ROS", "mpw2plyp",  colors[11], "-.", "ROMPWLYP", "iop(3/76=0450005500,3/78=0750007500)", "8/5=-1,10=90/1;9/16=-3/6;"),

    #G03LevelOfTheory("", "", "", ""),
    #G03LevelOfTheory("", "", "", "", ""),
]

lots = dict(((lot.spin, lot.label), lot) for lot in lots_list)

def get_lot(label, mult, restricted=None):
    if (mult > 1):
        if restricted is not True:
            result = lots.get(("DEF", label))
        else:
            result = lots.get(("ROS", label))
    else:
        if restricted is not False:
            result = lots.get(("DEF", label))
        else:
            result = lots.get(("UCS", label))
    if result is None:
        raise KeyError("Could not find a level of theory matching label=%s, mult=%s, restricted=%s" % (mult, label, restricted))
    return result



class G03Basis(object):
    def __init__(self, label, name, diffuse=False):
        self.label = label
        self.name = name
        self.diffuse = diffuse

basisses = [
    G03Basis("sto-3g", "STO-3G"),
    G03Basis("6-31gd", "6-31G(d)"),
    G03Basis("6-311+g3df2p", "6-311+G(3df,2p)", True),
]

basisses = dict((basis.label, basis) for basis in basisses)

def get_basis(label):
    result = basisses.get(label)
    if result is None:
        raise KeyError("Could not find a basis set matching label=%s" % label)
    return result
