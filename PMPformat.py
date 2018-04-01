#! /usr/bin/env python

'''
Module to generate and read the PMP Format file. 

Notice: If using PDB/PQR format as input, you should use a 
      reference molecule to define bond type.

TODO:
1. The bond order now is saving in CONECT, when double bond with two connected atoms:
   > CONECT   1   2   2
   This format can be read by Rdkit as PDB format normally. However, other software?  
   Is it essential to save bond order in 81~84? 

Update: 2018.4.1 by Zhixiong Zhao
XtalPi Inc.
'''

#%%
__DEBUG=False

import os,sys,math,StringIO
from collections import OrderedDict
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.RDLogger import logger

logger = logger()

class PMPFormator(object):
    '''PMP Formator'''
    _mol = None
    _pmp = None

    def __init__(self, pmpfile=None, pdbfile=None, mol2file=None, smiles=None):
        '''Init by mol2file or pdb file with smiles/mol2file'''
        if pmpfile:
            self.MolFromPMPFile(pmpfile)
        elif pdbfile and mol2file:
            self.MolFromPDBFile(pdbfile, refmol2file=mol2file)
        elif pdbfile and smiles:
            self.MolFromPDBFile(pdbfile, refsmiles=smiles)
        elif mol2file:
            self.MolFromMol2File(mol2file)
        elif smiles:
            self.MolFromSmiles(smiles)


        #mol=Chem.MolFromPDBFile(pdbfile)
        #molc=mol.GetConformer(0)
        #x, y, z=molc.GetAtomPosition(1)
        #print Chem.MolToSmiles(mol)
        pass

    def GetMol(self):
        ''' Return saved molecule '''
        return self._mol
    def GetPMP(self):
        ''' Return saved PMP format '''
        return self._pmp

    def GetValueByType(self,value, ptype='s'):
        '''Get Value based on single ptype'''
        if ptype =='s':
            return str(value)
        elif ptype =='i' or ptype =='d':
            return int(value)
        elif ptype =='f' :
            return float(value)

    def SetAtomValue(self, atomidx, value, pname, ptype='s'):
        '''
        Atom idx, value, property name, property type'''
        atom = self._mol.GetAtomWithIdx(atomidx)
        if ptype =='s':
            atom.SetProp(pname,str(value).strip())
        elif ptype =='i' or ptype =='d':
            atom.SetIntProp(pname,int(value))
        elif ptype =='f' :
            atom.SetDoubleProp(pname,float(value))

    def SetMolValue(self, value, pname, ptype='s'):
        '''value, property name, property type'''
        if ptype =='s':
            self._mol.SetProp(pname,str(value).strip())
        elif ptype =='i' or ptype =='d':
            self._mol.SetIntProp(pname,int(value))
        elif ptype =='f' :
            self._mol.SetDoubleProp(pname,float(value))

    def SetAtomsProp(self, values, pname, ptype='s', plen=6,floatPoint=3):
        '''Set Atom property : values, property name, property type, string length'''
        if len(values) != self._mol.GetNumAtoms():
            logger.warning('SetAtomsProp: '+self._mol.GetProp('Filename')+'The length of values not equal to  atom number!')
            return
        if not self._pmp:
            logger.warning('SetAtomsProp: No PMP information!')
            return            
        lines=self._pmp.splitlines()
        newpmp=""
        findLabel = False
        plen = str(plen)
        if ptype == 'str' or ptype == 'string' or ptype == 's': 
            ptype = 's'
        elif ptype == 'float' or ptype == 'double' or ptype =='f': 
            ptype = 'f'
        elif ptype == 'int' or ptype == 'long' or ptype=='i' or ptype=='d': 
            ptype = 'd'
        else :
            ptype = 's'
        ## Formatting string format
        if ptype!='f':
            sformat = "{:>"+plen+ptype+"}"
        else:
            sformat = "{:>"+plen+"."+str(floatPoint)+ptype+"}"
        atomcount = 0
        for line in lines:
            if not findLabel and line[:17]=="REMARK 666 LABELS":
                newpmp += line+"|"+pname+":"+plen+ptype+'\n'
            elif line[:6]=="ATOM  " or line[:6]=="HETATM":
                newpmp += line+sformat.format(values[atomcount])+'\n'
                self.SetAtomValue(atomcount, values[atomcount],pname, ptype)
                atomcount += 1
            else:
                newpmp += line+'\n'
        self._pmp = newpmp
        #print "Successfully add atomic property: "+pname+":"+plen+ptype

    def SetMolProp(self, value, pname, ptype='s', plen=6,floatPoint=3):
        if not self._pmp:
            logger.warning('SetAtomsProp: No PMP information!')
            return            
        lines=self._pmp.splitlines()
        plen = str(plen)
        newpmp=""
        if ptype == 'str' or ptype == 'string' or ptype == 's': 
            ptype = 's'
        elif ptype == 'float' or ptype == 'double' or ptype =='f': 
            ptype = 'f'
        elif ptype == 'int' or ptype == 'long' or ptype=='i' or ptype=='d': 
            ptype = 'd'
        else :
            ptype = 's'
        ## Formatting string format
        if ptype!='f':
            sformat = "{:>"+plen+ptype+"}"
        else:
            sformat = "{:>"+plen+"."+str(floatPoint)+ptype+"}"
        for line in lines:
            if line[:17] == "REMARK 666 LABELS":
                newpmp += "REMARK 111 MOLP: "+pname+": "+ptype+": "+sformat.format(value)+"\n"
                newpmp += line+'\n'
                self.SetMolValue(value, pname, ptype)
            else:
                newpmp += line+'\n'
        self._pmp = newpmp
        #print "Successfully add molecular property: "+pname+":"+plen+ptype

    def PMPfromPDBblock(self, PDBblock):
        '''Convert PDB Block to PMP format.
        PDB block can be generated from Chem.MolToPDBBlock(mol)
        '''
        lines = PDBblock.splitlines()
        connectDict={}
        connectAtoms=[]
        for line in lines:
            if line[:6] == 'CONECT':
                tmp=line[6:].rstrip()
                # Maybe some bug
                if len(tmp)%5 != 0:
                    logger.warning("PMPfromPDBblock(): Connect length maybe something wrong!")
                if len(tmp)/5 >= 2:
                    atomid=tmp[:5].strip()
                    connectAtoms.append(atomid)
                    connectDict[atomid]=[tmp[5:10].strip()]
                    for i in range(len(tmp)/5-2):
                         connectDict[atomid].append(tmp[10+5*i:15+5*i].strip())
                    connectDict[atomid]=list(OrderedDict.fromkeys(connectDict[atomid]))
                else:
                    logger.warning("PMPfromPDBblock(): Connect length maybe too short!")
        newblock = ""
        firstline = False
        for line in lines:
            if line[:6]=="ATOM  " or line[:6]=="HETATM":
                if not firstline:
                    newblock += "REMARK 666 LABELS|PMP_INFO:80s\n"
                    firstline = True
                atomid = line[6:11].strip()
                if connectDict.has_key(atomid):
                    connects = connectDict[atomid]
                    bondnum = len(connects)
                    newblock += line[:55]+str(bondnum)
                    if bondnum > 4: 
                        newblock += " "*20
                    else:
                        for i in range(4):
                            if i <= bondnum -1:
                                newblock += "%5s" % connects[i]
                            else:
                                newblock += " "*5
                    newblock += line[76:]+"\n"
                else:
                    newblock += line[:55]+"0"+" "*20+line[76:]+"\n"
            else:
                newblock += line+"\n"
        return newblock

    def MolMatchBondBySmiles(self,smiles, mol=None):
        '''Using ref smiles to modify the bond type of given molecule "mol".
        If "mol" not given, using the inner molecule and update it!
        ''' 
        # TODO: Some bug may be wrong even using 
        #       "AllChem.RemoveHs(,updateExplicitCount=True)"
        #       to remove H. If smiles contain explicit H, can't match!!
        #
        refmol = Chem.MolFromSmiles(smiles, sanitize = True)
        refmol = Chem.MolFromSmiles(Chem.MolToSmiles(refmol), sanitize = True)
        if mol:
            mdone  = AllChem.AssignBondOrdersFromTemplate(refmol,mol)
            return mdone
        elif self._mol:
            mdone  = AllChem.AssignBondOrdersFromTemplate(refmol,self._mol)
            pdbblock = Chem.MolToPDBBlock(mdone, flavor = 4)
            self._pmp = self.PMPfromPDBblock(pdbblock)
            self._mol = mdone
            return mdone
        else:
            logger.warning("MolMatchBondBySmiles: No input molecule was given!")
            return None

    def MolMatchBondByMol2File(self,mol2file, mol=None):
        '''Using ref mol2 file to modify the bond type of given molecule "mol".
        If "mol" not given, using the inner molecule and update it!
        ''' 
        refmol = Chem.MolFromMol2File(mol2file, sanitize = True, removeHs=False)
        if mol:
            mdone  = AllChem.AssignBondOrdersFromTemplate(refmol,mol)
            return mdone
        elif self._mol:
            mdone  = AllChem.AssignBondOrdersFromTemplate(refmol,self._mol)
            pdbblock = Chem.MolToPDBBlock(mdone, flavor = 4)
            self._pmp = self.PMPfromPDBblock(pdbblock)
            self._mol = mdone
            return mdone
        else:
            logger.warning("MolMatchBondByMol2File: No input molecule was given!")
            return None

    def MolFromMol2File(self, fname):
        '''Read the molecule by mol2 file. 
        Save molecule in object. '''
        mol = Chem.MolFromMol2File(fname, sanitize=True, removeHs=False)
        pdbblock = Chem.MolToPDBBlock(mol, flavor = 4)
        self._pmp = self.PMPfromPDBblock(pdbblock)
        self._mol = mol
        self.SetMolProp(fname, 'Filename', ptype='s', plen=10)
        return mol

    def MolFromSmiles(self, smiles):
        '''Read the molecule by smiles string. 
        Save molecule in object. '''
        mol = Chem.MolFromSmiles(smiles, sanitize = True)
        pdbblock = Chem.MolToPDBBlock(mol, flavor = 4)
        self._pmp = self.PMPfromPDBblock(pdbblock)
        self._mol = mol
        return mol

    def MolFromPDBFile(self, fname, refmol2file=None, refsmiles=None):
        ''' Read molecule from PDB/PQR file. 
        You should load the ref mol2 file or smiles at the same time!
        !!Notice: If no reference file/smiles, the bond order may be wrong!
        '''
        mol = Chem.MolFromPDBFile(fname, sanitize=True, removeHs=False)
        if refsmiles:
            mol = self.MolMatchBondBySmiles(refsmiles, mol = mol)
        elif refmol2file:
            mol = self.MolMatchBondByMol2File(refmol2file, mol = mol)
        pdbblock = Chem.MolToPDBBlock(mol, flavor = 4)
        self._pmp = self.PMPfromPDBblock(pdbblock)
        self._mol = mol
        self.SetMolProp(fname, 'Filename', ptype='s', plen=10)
        return mol

    def MolFromPMPFile(self, fname):
        ''' Read the PMP File '''
        mol = self.MolFromPDBFile(fname)
        self._mol = mol
        f = open(fname)
        lines = f.readlines()
        f.close()
        findLabel = False
        labelinfo = {}
        labels = []
        atomcount = 0
        for line in lines:
            ## Read molecular property
            if not findLabel and line[:15]=="REMARK 111 MOLP":
                tmps=line.strip().split(":")
                if len(tmps) >= 4:
                    self.SetMolValue(tmps[3].strip(), tmps[1].strip(), tmps[2].strip())
            ## Read Atom Label
            elif not findLabel and line[:17]=="REMARK 666 LABELS":
                findLabel = True
                tmps=line.split('|')
                if len(tmps) > 2:
                    for label in tmps[2:]:
                        ltmp=label.strip().split(':')
                        labelName=ltmp[0]
                        labelValue=ltmp[1][:-1]
                        labelType=ltmp[1][-1]
                        labels.append(labelName)
                        labelinfo[labelName]=(labelValue,labelType)
            elif line[:6]=="ATOM  " or line[:6]=="HETATM":
                atom = self._mol.GetAtomWithIdx(atomcount)
                if labelinfo:
                    place = 80
                    for label in labels:
                        lvalue=labelinfo[label][0]
                        ltype =labelinfo[label][1]                        
                        self.SetAtomValue(atomcount,line[place:place+int(lvalue)], label, ltype)
                        place += int(lvalue)
                atomcount += 1
        self.SetMolProp(fname, 'Filename', ptype='s', plen=10)
        return self._mol

    def MolToPMPFile(self, fname):
        ''' Write out as PMP file'''
        if self._pmp:
            with open(fname,'w') as f:
                f.write(self._pmp)         
    
    def MolToPMPBlock(self):
        ''' Get the PMP Block'''
        return self._pmp



if __DEBUG:
    os.chdir("/home/hom/Desktop/DailyWork/Cheminfo/20180329_PMPformat/")
    #print os.path.realpath('.')
    pmpf = PMPFormator()

    print "Initial from Mol2 file "
    mol = pmpf.MolFromMol2File('~data/test.mol2')
    print Chem.MolToSmiles(mol)
    #print pmpf._pmp

    print "Initial from Smiles"
    mol2 = pmpf.MolFromSmiles('CC(=O)[O-].[Na+]')
    print Chem.MolToSmiles(mol2)
    #print pmpf._pmp

if __DEBUG:
    print "Initial from PDB based on mol2"
    m1 = pmpf.MolFromPDBFile('~data/test.pqr', refmol2file="~data/test.mol2")
    print Chem.MolToSmiles(m1)
    #print pmpf._pmp

## Explict H in Ref Mol from Smiles will get wrong result 
if __DEBUG:
    print "Initial from PDB based on smiles without explicit H"
    m2 = pmpf.MolFromPDBFile('~data/test.pqr', refsmiles="c1sccc1")
    print Chem.MolToSmiles(m2)
    print "Initial from PDB based on smiles with explicit H"
    m2 = pmpf.MolFromPDBFile('~data/test.pqr', refsmiles="[H]c1sc([H])c([H])c1[H]")
    print Chem.MolToSmiles(m2)
    #print pmpf._pmp

## Test Reading GDMA data
if __DEBUG:
    GDMAdata, MolGDMA=ProcessGDMAdata('~data/test.gdma',multipole=2)
    for item,value in GDMAdata.items():
        pmpf.SetAtomsProp(value, item, ptype='f', plen=10,floatPoint=6)
    for item,value in MolGDMA.items():
        pmpf.SetMolProp(value, item, ptype='f', plen=10,floatPoint=6)    
    print pmpf.MolToPMPBlock()
    print pmpf._mol.GetPropsAsDict()
    print pmpf._mol.GetAtomWithIdx(0).GetPropsAsDict()

## Test Reading Atomic Solvation Energy
if __DEBUG:
    soldata = ProcessAtomSolEng('~data/AtomSoleng.txt')
    pmpf.SetAtomsProp(soldata, "AtomSolEng", ptype='f', plen=10,floatPoint=6)

## Test Reading From PQRTA
if __DEBUG:
    atomtype="AT_gaff"
    AtomDatas, MolDatas= ProcessPQRTA('~data/test.pqrta', atomtype=atomtype,charge="resp", radius="mbondi")
    for item,value in AtomDatas.items():
        if item != atomtype:
            pmpf.SetAtomsProp(value, item, ptype='f', plen=10,floatPoint=4)
        else:
            pmpf.SetAtomsProp(value, item, ptype='s', plen=6)            
    for item,value in MolDatas.items():
        pmpf.SetMolProp(value, item, ptype='f', plen=10,floatPoint=4)  

## Test Reading SYBYL Type and MIBPB result
if __DEBUG:
    sybyltype = ProcessMol2AtomType('~data/test.mol2')
    pmpf.SetAtomsProp(sybyltype, "AT_sybyl", ptype='s', plen=6)
    mibpbOut = ProcessMIBPBresults('~data/mibpb5.log')
    pmpf.SetMolProp(mibpbOut, 'ElecSolvEng', ptype='f', plen=10,floatPoint=6)  

## Test Writing and Reading PMP file
if __DEBUG:
    pmpf.MolToPMPFile('test.pmp')
    pmpf2 = PMPFormator(pmpfile='test.pmp')
    print Chem.MolToSmiles(pmpf2._mol)
    print pmpf2._mol.GetPropsAsDict()
    print pmpf2._mol.GetAtomWithIdx(0).GetPropsAsDict()