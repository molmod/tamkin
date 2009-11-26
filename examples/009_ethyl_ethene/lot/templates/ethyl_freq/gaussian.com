%chk=gaussian.chk
%nproc=1
%mem=1000MB
#p ${lot_basis} freq(noraman) gfinput iop(6/7=3) guess=read maxdisk=5GB

Frequency computation

0 2
${atom_str} 


