from openeye import oechem

# mol = oechem.OEGraphMol()
# # del mol
# # oechem.OESmilesToMol(mol, 'c1ccccc1')
# if oechem.OESmilesToMol(mol, 'c1ccccc1'):
#     pass
# else:
#     print('SMILES string was invalid')
#
# print('Done')

# mol = oechem.OEGraphMol()
# if not oechem.OEParseSmiles(mol, 'C1=CC=CC=C1'):
#     print('SMILES string was invalid')
#
# print('Number of aromatic atoms =', oechem.OECount(mol, oechem.OEIsAromaticAtom()))
# oechem.OEAssignAromaticFlags(mol)
# print('Number of aromatic atoms =', oechem.OECount(mol, oechem.OEIsAromaticAtom()))

# from openeye import oechem
#
# mol = oechem.OEGraphMol()
#
# # oechem.OESmilesToMol(mol, 'c1ccccc1')
# # print('Number of benzene atoms:', mol.NumAtoms())
# #
# # oechem.OESmilesToMol(mol, 'c1ccccc1O')
# # print('Number of phenol atoms:', mol.NumAtoms())
#
# oechem.OEParseSmiles(mol, 'c1ccccc1')
# print(mol.NumAtoms())
# mol.Clear()
#
# oechem.OEParseSmiles(mol, 'c1ccccc1O')
# print(mol.NumAtoms())

# from openeye import oechem
#
# mol = oechem.OEGraphMol()
# print(oechem.OESmilesToMol(mol, 'C1=CC=CC=C1'))
#
# print("Canonical isomeric SMILES is", oechem.OEMolToSmiles(mol))

from openeye import oechem
import sys

# mol = oechem.OEGraphMol()
#
# for smi in sys.stdin:
#     mol.Clear()
#     smi = smi.strip()
#     if oechem.OEParseSmiles(mol, smi):
#         oechem.OEAssignAromaticFlags(mol)
#         print(oechem.OECreateCanSmiString(mol))
#     else:
#         oechem.OEThrow.Warning("%s is an invalid SMILES!" % smi)


mol = oechem.OEGraphMol()
oechem.OESmilesToMol(mol, 'c1ccnc(c1)O')
print(oechem.OEMolToSTDInChI(mol))
print(oechem.OEMolToInChI(mol))

print(oechem.OEMolToSTDInChIKey(mol))


