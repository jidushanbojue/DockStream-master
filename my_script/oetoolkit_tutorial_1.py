from openeye import oechem

def main():
    mol = oechem.OEGraphMol()
    oechem.OEParseSmiles(mol, 'c1ccccc1CCCBr')
    print("mol hjas {0} atoms".format(mol.NumAtoms()))


if __name__ == '__main__':
    main()