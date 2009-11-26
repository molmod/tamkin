%chk=gaussian.chk 
%nproc=1
%mem=1000MB
#p ${lot_basis} opt gfinput iop(6/7=3) maxdisk=5GB

Geometry optimization

0 1
${atom_str}

 
