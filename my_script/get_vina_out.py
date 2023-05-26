import os

if __name__ == '__main__':
    base_dir = '/data/baiqing/PycharmProjects/DockStream-master/result'

    name_list = ['5T35',
                 '6BN7',
                 '6BOY',
                 '6HAY',
                 '6HR2',
                 '6W7O',
                 '6ZHC',
                 '7JTO',
                 '7JTP',
                 '7KHH',
                 '7Q2J',
                 'BRD9',
                 'BAF',
                 'BTK']

    for name in  name_list:
        work_dir = os.path.join(base_dir, name)
        os.chdir(work_dir)
        # cmd = 'obabel -ipdbqt ' + name.lower()+'_crystal_ligand_ligprep_docking_out.pdbqt -osdf -O ' + name.lower()+'_crystal_ligand_ligprep_docking_out.sdf'
        cmd = 'obabel -ipdbqt 0000_ADV_out.pdbqt -osdf -O 0000_ADV_out.sdf'
        print(cmd)

        os.system(cmd)
