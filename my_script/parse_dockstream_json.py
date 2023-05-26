import os
import json
from rdkit import Chem

def parse(file):

    sdf_idx_list = []

    with open(file) as f:
        json_input = f.read().replace('\r', '').replace('\n', '')

        content = json.loads(json_input)

        sdf_path = content['docking']['ligand_preparation']['embedding_pools'][0]['input']['input_path']

        sdf_suppl = Chem.SDMolSupplier(sdf_path)

        for mol in sdf_suppl:
            sdf_idx_list.append(mol.GetProp('_Name'))
            # print(mol.GetProp('_Name'))
            # print(Chem.MolToSmiles(mol))



    # try:
    #     return json.loads(json_input)
    # except (ValueError, KeyError, TypeError) as e:
    #     print(f"JSON format error in file ${path}: \n ${e}")

    # content = json.loads(file)
    print('Done')


if __name__ == '__main__':
    base_dir = '/data/baiqing/PycharmProjects/DockStream-master/result/BAF'
    json_file = os.path.join(base_dir, 'BAF_scoring_linkinvent_7_20_scaffold_memory_top100.json')

    parse(json_file)