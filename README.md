# Data for Delaney logS in PMP format

The PMP format contains following information:

- `HF/6-31G*` optimized structure
- `HF/6-31G*` based RESP charge
- Molecular ELectrostatic Solvation Energy and Atomic ELectrostatic Solvation Energy (RESP charge, mbondi radius)
- Dipole and Quadrapole based on GDMA analysis (H radius at default 0.325)
- Molecular Area/Volume (0.2 grid size, 1.4 probe radius)
- Atomic Area

## Usage

```python
import PMPformat as PMP
```

### Read file as molecule

- Open PMP molecule file
```python
pmpf = PMP.PMPFormator(pmpfile="Data/Delaney_1.pmp")
```

- Get the molecule.  
> You can replace "MolFromSmiles" to this function to load molecule.

```python
mol = pmpf.GetMol()
```
### Get molecular or atomic features
- Get Property/Feature on the Molecular Level

```python
mol_prop = mol.GetPropsAsDict()
print mol_prop
```

-  Get a certain feature by name.   
> Recommend use the `MoleArea/Volume/ElecSolvEng/Dipole` and so on.  
> I don't recommend to use element based area, such as `Area_F`, `Area_S` et al. But you can also try.

```python
mol_prop['MolArea']
```

An alternative method is to use `GetProp`, `GetIntProp` or `GetDoubleProp` methods to obtain properties when you exactly know the data type for the proporty.

```python
mol_area = mol.GetDoubleProp('MolArea')
```

- Get Property/Feature on the Atom Level
> The atom type `AT_gaff`, `AT_sybyl` somehow can be represented by NGF and it's string based data. If you want to use it, you can use int value to represent it.

```python
atom0 = mol.GetAtomWithIdx(0)
atom_prop = atom0.GetPropsAsDict()
print atom_prop

atom_prop['AtomArea']

## or using exactly method

features = [
    atom0.GetProps('AT_gaff'), atom0.GetDoubleProps('AtomArea')
]

```




