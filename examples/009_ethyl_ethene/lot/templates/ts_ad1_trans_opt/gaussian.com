%chk=gaussian.chk
%nproc=1
%mem=1000MB
#p ${lot_basis} opt(ts,modredundant,noeigentest) gfinput iop(6/7=3) maxdisk=5GB

Geometry optimization a transition state

0 2
${atom_str}



