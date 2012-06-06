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
"""
This script finds the overlap between the 7th mode and the 7th mode
and writes the result to a file delta.overlaps.7-7.csv.
First line/column contains the frequencies.
Next column contains the overlap (still to be squared).

Script does this for every delta.*overlap.csv file present in the directory.

usage:  python %progr
"""


def get_freqs_and_overlap_7_7(filename):
    f = file(filename)
    print "-"*20
    print "Reading file...", filename

    # read 7th freq on first line (first element is 0, so is 8th freq)
    for line in f:
        words = line.split(";")
        #freq1 = float(words[7])  # !!!!!!!!!!!!!
        freq1 = float(words[-1])
        break

    # get 7th freq in first column (we already skipped first line)
    # and overlap 7th-7th
    count = 0
    for line in f:
        count +=1
        if count==7:
            words = line.split(";")
            freq2 = float(words[0])
            #overlap = float(words[7])  # !!!!!!!!!
            overlap = float(words[-1])
            break

    f.close()
    print "freq1: ", freq1
    print "freq2: ", freq2
    print "overlap: ", overlap

    return freq1, freq2, overlap


import sys, glob

freqs1 = []
freqs2 = []
overlaps = []
filenames = glob.glob("delta.*overlap.csv")

filenames.sort()
print "Considered files: ", filenames
for filename in filenames:
    #y = str(x)
    #if x < 10:
    #    y = "0"+str(x)
    freq1, freq2, overlap = get_freqs_and_overlap_7_7(filename)
    freqs1.append(freq1)
    freqs2.append(freq2)
    overlaps.append(overlap)

conv = 7.25163277859109E-07  # conversion factor such that frequencies are printed in 1/cm

filename_out = "delta.overlaps.7-7.csv"
f = file(filename_out,"w+")
for i in range(len(freqs1)):
    print >> f, filenames[i]+";"+str(freqs1[i]/conv)+";"+str(freqs2[i]/conv)+";"+str(overlaps[i])
f.close()
print "file written:", filename_out
