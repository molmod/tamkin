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


from tamkin import *

corfile = "open/adk.open.new.cor"
fixedfile = "fixed_files/fixed."


for i in range(1,6):
    blocks = create_blocks_peptide_charmm(corfile, BlocksPeptideMBH("RTB", blocksize=i))
    blocks_write_to_file(blocks, fixedfile+str(i)+".txt")

for i in range(1,6):
    subs = create_subs_peptide_charmm(corfile, SubsPeptideVSA(frequency=i))
    subs_write_to_file(subs, fixedfile+str(i+5)+".txt")


blocks = create_blocks_peptide_charmm(corfile, BlocksPeptideMBH("RHbending"))
blocks_write_to_file(blocks, fixedfile+"11.txt")

blocks = create_blocks_peptide_charmm(corfile, BlocksPeptideMBH("dihedral"))
blocks_write_to_file(blocks, fixedfile+"12.txt")

blocks = create_blocks_peptide_charmm(corfile, BlocksPeptideMBH())
blocks_write_to_file(blocks, fixedfile+"13.txt")


#for x in 6 7 8 9 10 ; do
#python analyse-adk.py --job vsa --filecor open/adk.open.new.cor  --filehess open/adk.open.hess.full \
#                      --filechk chkfiles/adk.open.vsa.$x.chk  \
#                      --filefixed fixed_files/fixed.$x.txt

#done
