%chk=gaussian.chk
%nproc=1
%mem=1000MB
#p ${lot_basis} opt(ts,modredundant,noeigentest) gfinput iop(6/7=3) guess=read maxdisk=5GB

Rotational scan of the methyl group

0 2
${atom_str}

5 1 2 9 S 72 5.0


