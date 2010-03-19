#!/usr/bin/env python


# Import the tamkin library.
from tamkin import *
from molmod import * # for the units
from molmod.isotopes import ame2003

def compute_k(mol_react, mol_trans):
    # Perform normal mode analysis on the three molecules
    nma_react = NMA(mol_react, ConstrainExt())
    nma_trans = NMA(mol_trans, ConstrainExt())
    # Construct the two partition functions.
    pf_react = PartFun(nma_react, [ExtTrans(), ExtRot(1), Vibrations(classical=True, freq_scaling=0.9085)])
    pf_trans = PartFun(nma_trans, [ExtTrans(), ExtRot(1), Vibrations(classical=True, freq_scaling=0.9085)])

    return compute_rate_coeff([pf_react], pf_trans, 303)



old_mol_react = load_molecule_g03fchk("reactant.fchk")
old_mol_trans = load_molecule_g03fchk("trans.fchk")
k_orig = compute_k(old_mol_react, old_mol_trans)
print "Original rate coefficient =", k_orig/((meter**3/mol)/second), "(m**3/mol)/s"


# does not work:
# mol_react1.masses[0] = 13*amu

new_masses = old_mol_react.masses.copy()
new_masses[20] = ame2003.masses[7][15]
new_mol_react = old_mol_react.copy_with(masses=new_masses)
new_mol_trans = old_mol_trans.copy_with(masses=new_masses)

k_new = compute_k(new_mol_react, new_mol_trans)
print "New rate coefficient =", k_new/((meter**3/mol)/second), "(m**3/mol)/s"

print "Ratio =", k_orig/k_new
