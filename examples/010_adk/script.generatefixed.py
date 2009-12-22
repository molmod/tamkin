#  #!/bin/usr/python


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
