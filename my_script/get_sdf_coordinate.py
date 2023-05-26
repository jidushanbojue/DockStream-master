import os
from rdkit import Chem


input_sdf_file = '/data/baiqing/PycharmProjects/Reinvent-master-3.2/result/LINK_invent/BRD9/temp_test/rearranged_h_mini.sdf'
input_pdbqt_file = '/tmp/tmpelyqu88v/0011dxxxn.pdbqt'

output_pdbqt_file = '/tmp/tmpelyqu88v/04amhlv6m.pdbqt'
output_sdf_file = '/tmp/tmpelyqu88v/0w1hn07af.sdf'

import re


def get_input_pdbqt_atom_idx(reference_coord_dic, pdbqt_file):
    """
    reference_coord_list: [(x, y, z), (x, y, z)]
    """
    f = open(pdbqt_file)
    lines = f.readlines()
    f.close()

    input_idx_dic = {}

    for idx, line in enumerate(lines):
        if line.startswith('ATOM'):
            temp_list = re.split('\s+', line)
            # coord = (temp_list[5], temp_list[6], temp_list[7])
            tgt_x, tgt_y, tgt_z = float(temp_list[5]), float(temp_list[6]), float(temp_list[7])
            for ref_coord, isotope in reference_coord_dic.items():
                x, y, z = ref_coord
                if x == tgt_x and y == tgt_y and z == tgt_z:
                    print(ref_coord)
                    # input_idx_list.append(idx)
                    input_idx_dic[idx] = isotope
    # print('Done')
    return input_idx_dic

def get_output_pdbqt_coord(input_idx_dic, pdbqt_file):
    f = open(pdbqt_file)
    lines = f.readlines()
    f.close()

    coord_dic = {}
    for idx, line in enumerate(lines):
        if idx in input_idx_dic:
            temp_list = re.split('\s+', line)
            x, y, z = float(temp_list[5]), float(temp_list[6]), float(temp_list[7])
            # coord_dic.append((x, y, z))
            coord_dic[(x, y, z)] = input_idx_dic[idx]
    return coord_dic

def get_output_sdf_atom_idx(coord_dic, sdf_file):
    suppl = Chem.SDMolSupplier(sdf_file)
    idx_dic = {}
    for mol in suppl:
        for atom in mol.GetAtoms():
            x, y, z = mol.GetConformer().GetAtomPosition(atom.GetIdx())
            x, y, z = round(x, 3), round(y, 3), round(z, 3)
            for coord, isotope in coord_dic.items():
                ref_x, ref_y, ref_z = coord
                if x == ref_x and y == ref_y and z == ref_z:
                    # idx_list.append(atom.GetIdx())
                    idx_dic[atom.GetIdx()] = isotope
    return idx_dic


suppls = Chem.SDMolSupplier(input_sdf_file)
for mol in suppls:
    attachment_atom_coord_dic = {}
    for atom in mol.GetAtoms():
        if atom.GetIsotope() != 0:
            print(atom.GetSymbol())
            x, y, z = mol.GetConformer().GetAtomPosition(atom.GetIdx())
            x, y, z = round(x, 3), round(y, 3), round(z, 3)
            attachment_atom_coord_dic[(x, y, z)] = atom.GetIsotope()

    input_pdbqt_idx_dic = get_input_pdbqt_atom_idx(attachment_atom_coord_dic, input_pdbqt_file)
    output_pdbqt_coord_dic = get_output_pdbqt_coord(input_pdbqt_idx_dic, output_pdbqt_file)
    output_sdf_dic = get_output_sdf_atom_idx(output_pdbqt_coord_dic, output_sdf_file)

    out_sdf_suppl = Chem.SDMolSupplier(output_sdf_file)
    for out_sdf_mol in out_sdf_suppl:
        for atom in out_sdf_mol.GetAtoms():
            if atom.GetIdx() in output_sdf_dic:
                atom.SetIsotope(output_sdf_dic[atom.GetIdx()])

        w = Chem.SDWriter(output_sdf_file)
        w.write(out_sdf_mol)
        w.close()

    # for atom in mol.GetAtoms():
    #     if atom.GetIdx():
    # print('Done')

def add_atommap(sublists, tmp_input_paths, tmp_output_pdbqt_paths, tmp_output_paths):
    for src_sdf_mol, input_pdbqt_path, output_pdbqt_path, output_sdf_path in zip(sublists, tmp_input_paths, tmp_output_pdbqt_paths, tmp_output_paths):
        # print(src_sdf_mol)
        input_sdf_mol = deepcopy(src_sdf_mol[0].get_molecule())

        attachment_atom_coord_dic = {}
        for atom in input_sdf_mol.GetAtoms():
            if atom.GetIsotope() != 0:
                # print(atom.GetSymbol())
                x, y, z = input_sdf_mol.GetConformer().GetAtomPosition(atom.GetIdx())
                x, y, z = round(x, 3), round(y, 3), round(z, 3)
                attachment_atom_coord_dic[(x, y, z)] = atom.GetIsotope()

        input_pdbqt_idx_dic = get_input_pdbqt_atom_idx(attachment_atom_coord_dic, input_pdbqt_path)
        output_pdbqt_coord_dic = get_output_pdbqt_coord(input_pdbqt_idx_dic, output_pdbqt_path)
        output_sdf_dic = get_output_sdf_atom_idx(output_pdbqt_coord_dic, output_sdf_path)

        try:
            out_sdf_suppl = Chem.SDMolSupplier(output_sdf_path)
            for out_sdf_mol in out_sdf_suppl:
                for atom in out_sdf_mol.GetAtoms():
                    if atom.GetIdx() in output_sdf_dic:
                        atom.SetIsotope(output_sdf_dic[atom.GetIdx()])

                w = Chem.SDWriter(output_sdf_path)
                w.write(out_sdf_mol)
                w.close()
                print('Done')
        except AttributeError:
            continue
        except OSError:
            continue



