
## Usage

import PMPformat as PMP

## Open PMP molecule file
pmpf = PMP.PMPFormator(pmpfile="Data/Delaney_1.pmp")

## Get the molecule. You can replace "MolFromSmiles" to this function to load molecule.
mol = pmpf.GetMol()

## Get Property/Feature on the Molecular Level
mol_prop = mol.GetPropsAsDict()
print mol_prop
##### Get a certain feature by name. 
## Recommend use the Molecular Area/Volume/ElecSolvEnergy/Dipole and so on
mol_prop['MolArea']

## Get Property/Feature on the Atom Level
atom0 = mol.GetAtomWithIdx(0)
atom_prop = atom0.GetPropsAsDict()
print atom_prop

atom_prop['AtomArea']





