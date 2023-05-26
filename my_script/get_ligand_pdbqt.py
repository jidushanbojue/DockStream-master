import os

if __name__ == '__main__':
    reinvent_base_dir = '/data/baiqing/PycharmProjects/Reinvent-master-3.2/result/LINK_invent'
    dockstream_base_dir = '/data/baiqing/PycharmProjects/DockStream-master/result'

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
                 'BTK',
                 'BRD9',
                 'BAF']

    for name in name_list:
        sdf_file = os.path.join(reinvent_base_dir, name, 'initial_pose_generation', '0000.sdf')
        dockstream_work_dir = os.path.join(dockstream_base_dir, name)

        os.chdir(dockstream_work_dir)


        # cmd = 'obabel -isdf ' + name.lower()+'_crystal_ligand_ligprep_docking_out.pdbqt -osdf -O ' + name.lower()+'_crystal_ligand_ligprep_docking_out.sdf'
        cmd = 'obabel -isdf ' + sdf_file + ' -opdbqt -O 0000.pdbqt -h'
        print(cmd)

        os.system(cmd)