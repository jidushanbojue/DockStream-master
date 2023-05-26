import os
import time
from functools import wraps

def timethis(func):
    """
    Decorator that reports the execution time.
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)
        end = time.time()
        print(func.__name__, (end-start)/60)
        return result
    return wrapper

@timethis
def vina_docking(work_dir, cmd):
    os.chdir(work_dir)
    # cmd = '/data/baiqing/PycharmProjects/DockStream-master/DockStreamCommunity/AutoDock-Vina-1.2.3/build/linux/release/vina --receptor brd9_protein_apo.pdbqt --ligand BRD9_cystal_ligand_ligprep.pdbqt --center_x -35.44 --center_y 17.08 --center_z -20.58 --size_x 40 --size_y 40 --size_z 40 --out BRD9_cystal_ligand_ligprep_docking_out.pdbqt --seed 42 --num_modes 2'
    os.system(cmd)


if __name__ == '__main__':
    # work_dir_BRD9 = '/data/baiqing/PycharmProjects/DockStream-master/result/BRD9'
    # cmd = '/data/baiqing/PycharmProjects/DockStream-master/DockStreamCommunity/AutoDock-Vina-1.2.3/build/linux/release/vina --receptor brd9_protein_apo.pdbqt --ligand 0000.pdbqt --center_x -35.44 --center_y 17.08 --center_z -20.58 --size_x 40 --size_y 40 --size_z 40 --out 0000_ADV_out.pdbqt --seed 42 --num_modes 2'
    # vina_docking(work_dir_BRD9, cmd)
    #
    # work_dir_BAF = '/data/baiqing/PycharmProjects/DockStream-master/result/BAF'
    # cmd = '/data/baiqing/PycharmProjects/DockStream-master/DockStreamCommunity/AutoDock-Vina-1.2.3/build/linux/release/vina --receptor 6hax_apo.pdbqt --ligand 0000.pdbqt --center_x -38.83 --center_y 15.67 --center_z -28.64 --size_x 40 --size_y 40 --size_z 40 --out 0000_ADV_out.pdbqt --seed 42 --num_modes 2'
    # vina_docking(work_dir_BAF, cmd)
    #
    # work_dir_BTK = '/data/baiqing/PycharmProjects/DockStream-master/result/BTK'
    # cmd = '/data/baiqing/PycharmProjects/DockStream-master/DockStreamCommunity/AutoDock-Vina-1.2.3/build/linux/release/vina --receptor 6w8i_apo.pdbqt --ligand 0000.pdbqt --center_x -25.12 --center_y -21.29 --center_z 21.71 --size_x 40 --size_y 40 --size_z 40 --out 0000_ADV_out.pdbqt --seed 42 --num_modes 2'
    # vina_docking(work_dir_BTK, cmd)
    #
    #
    # work_dir_5T35 = '/data/baiqing/PycharmProjects/DockStream-master/result/5T35'
    # cmd = '/data/baiqing/PycharmProjects/DockStream-master/DockStreamCommunity/AutoDock-Vina-1.2.3/build/linux/release/vina --receptor 5t35_apo.pdbqt --ligand 0000.pdbqt --center_x 21.89 --center_y -28.24 --center_z -14.53 --size_x 40 --size_y 40 --size_z 40 --out 0000_ADV_out.pdbqt --seed 42 --num_modes 2'
    # vina_docking(work_dir_5T35, cmd)
    #
    # work_dir_6W7O = '/data/baiqing/PycharmProjects/DockStream-master/result/6W7O'
    # cmd = '/data/baiqing/PycharmProjects/DockStream-master/DockStreamCommunity/AutoDock-Vina-1.2.3/build/linux/release/vina --receptor 6w7o_apo.pdbqt --ligand 0000.pdbqt --center_x 22.38 --center_y 42.21 --center_z -5.93 --size_x 40 --size_y 40 --size_z 40 --out 0000_ADV_out.pdbqt --seed 42 --num_modes 2'
    # vina_docking(work_dir_6W7O, cmd)
    #
    # work_dir_6ZHC = '/data/baiqing/PycharmProjects/DockStream-master/result/6ZHC'
    # cmd = '/data/baiqing/PycharmProjects/DockStream-master/DockStreamCommunity/AutoDock-Vina-1.2.3/build/linux/release/vina --receptor 6zhc_apo.pdbqt --ligand 0000.pdbqt --center_x 27.59 --center_y 51.48 --center_z 17.79 --size_x 40 --size_y 40 --size_z 40 --out 0000_ADV_out.pdbqt --seed 42 --num_modes 2'
    # vina_docking(work_dir_6ZHC, cmd)
    #
    # work_dir_7JTP = '/data/baiqing/PycharmProjects/DockStream-master/result/7JTP'
    # cmd = '/data/baiqing/PycharmProjects/DockStream-master/DockStreamCommunity/AutoDock-Vina-1.2.3/build/linux/release/vina --receptor 7jtp_apo.pdbqt --ligand 0000.pdbqt --center_x 7.62 --center_y 7.19 --center_z -10.22 --size_x 40 --size_y 40 --size_z 40 --out 0000_ADV_out.pdbqt --seed 42 --num_modes 2'
    # vina_docking(work_dir_7JTP, cmd)
    #
    # work_dir_7Q2J = '/data/baiqing/PycharmProjects/DockStream-master/result/7Q2J'
    # cmd = '/data/baiqing/PycharmProjects/DockStream-master/DockStreamCommunity/AutoDock-Vina-1.2.3/build/linux/release/vina --receptor 7q2j_apo.pdbqt --ligand 0000.pdbqt --center_x -10.98 --center_y 6.47 --center_z 17.65 --size_x 40 --size_y 40 --size_z 40 --out 0000_ADV_out.pdbqt --seed 42 --num_modes 2'
    # vina_docking(work_dir_7Q2J, cmd)
    #
    # work_dir_6BN7 = '/data/baiqing/PycharmProjects/DockStream-master/result/6BN7'
    # cmd = '/data/baiqing/PycharmProjects/DockStream-master/DockStreamCommunity/AutoDock-Vina-1.2.3/build/linux/release/vina --receptor 6bn7_apo.pdbqt --ligand 0000.pdbqt --center_x 68.22 --center_y 40.96 --center_z 53.03 --size_x 40 --size_y 40 --size_z 40 --out 0000_ADV_out.pdbqt --seed 42 --num_modes 2'
    # vina_docking(work_dir_6BN7, cmd)
    #
    # work_dir_7KHH = '/data/baiqing/PycharmProjects/DockStream-master/result/7KHH'
    # cmd = '/data/baiqing/PycharmProjects/DockStream-master/DockStreamCommunity/AutoDock-Vina-1.2.3/build/linux/release/vina --receptor 7khh_apo.pdbqt --ligand 0000.pdbqt --center_x -11.48 --center_y 8.64 --center_z -16.94 --size_x 40 --size_y 40 --size_z 40 --out 0000_ADV_out.pdbqt --seed 42 --num_modes 2'
    # vina_docking(work_dir_7KHH, cmd)
    #
    # work_dir_6BOY = '/data/baiqing/PycharmProjects/DockStream-master/result/6BOY'
    # cmd = '/data/baiqing/PycharmProjects/DockStream-master/DockStreamCommunity/AutoDock-Vina-1.2.3/build/linux/release/vina --receptor 6boy_apo.pdbqt --ligand 0000.pdbqt --center_x 72.03 --center_y 38.67 --center_z 51.08 --size_x 40 --size_y 40 --size_z 40 --out 0000_ADV_out.pdbqt --seed 42 --num_modes 2'
    # vina_docking(work_dir_6BOY, cmd)
    #
    # work_dir_6HR2 = '/data/baiqing/PycharmProjects/DockStream-master/result/6HR2'
    # cmd = '/data/baiqing/PycharmProjects/DockStream-master/DockStreamCommunity/AutoDock-Vina-1.2.3/build/linux/release/vina --receptor 6hr2_apo.pdbqt --ligand 0000.pdbqt --center_x -27.68 --center_y -9.99 --center_z -22.06 --size_x 40 --size_y 40 --size_z 40 --out 0000_ADV_out.pdbqt --seed 42 --num_modes 2'
    # vina_docking(work_dir_6HR2, cmd)
    #
    # work_dir_6HAY = '/data/baiqing/PycharmProjects/DockStream-master/result/6HAY'
    # cmd = '/data/baiqing/PycharmProjects/DockStream-master/DockStreamCommunity/AutoDock-Vina-1.2.3/build/linux/release/vina --receptor 6hay_apo.pdbqt --ligand 0000.pdbqt --center_x -31.9 --center_y 13.95 --center_z -28.09 --size_x 40 --size_y 40 --size_z 40 --out 0000_ADV_out.pdbqt --seed 42 --num_modes 2'
    # vina_docking(work_dir_6HAY, cmd)
    #
    #
    # work_dir_7JTO = '/data/baiqing/PycharmProjects/DockStream-master/result/7JTO'
    # cmd = '/data/baiqing/PycharmProjects/DockStream-master/DockStreamCommunity/AutoDock-Vina-1.2.3/build/linux/release/vina --receptor 7jto_apo.pdbqt --ligand 0000.pdbqt --center_x -10.15 --center_y 3.09 --center_z 15.65 --size_x 40 --size_y 40 --size_z 40 --out 0000_ADV_out.pdbqt --seed 42 --num_modes 2'
    # vina_docking(work_dir_7JTO, cmd)

    work_dir_BRD9 = '/data/baiqing/PycharmProjects/DockStream-master/result/BRD9'
    cmd = '/data/baiqing/PycharmProjects/DockStream-master/DockStreamCommunity/AutoDock-Vina-1.2.3/build/linux/release/vina --receptor brd9_protein_apo.pdbqt --ligand 0000.pdbqt --center_x -35.44 --center_y 17.08 --center_z -20.58 --size_x 40 --size_y 40 --size_z 40 --out 0000_ADV_out.pdbqt --seed 42 --num_modes 2'
    vina_docking(work_dir_BRD9, cmd)

    work_dir_BAF = '/data/baiqing/PycharmProjects/DockStream-master/result/BAF'
    cmd = '/data/baiqing/PycharmProjects/DockStream-master/DockStreamCommunity/AutoDock-Vina-1.2.3/build/linux/release/vina --receptor 6hax_apo.pdbqt --ligand 0000.pdbqt --center_x -38.83 --center_y 15.67 --center_z -28.64 --size_x 40 --size_y 40 --size_z 40 --out 0000_ADV_out.pdbqt --seed 42 --num_modes 2'
    vina_docking(work_dir_BAF, cmd)

    work_dir_BTK = '/data/baiqing/PycharmProjects/DockStream-master/result/BTK'
    cmd = '/data/baiqing/PycharmProjects/DockStream-master/DockStreamCommunity/AutoDock-Vina-1.2.3/build/linux/release/vina --receptor 6w8i_apo.pdbqt --ligand 0000.pdbqt --center_x -25.12 --center_y -21.29 --center_z 21.71 --size_x 40 --size_y 40 --size_z 40 --out 0000_ADV_out.pdbqt --seed 42 --num_modes 2'
    vina_docking(work_dir_BTK, cmd)


    work_dir_5T35 = '/data/baiqing/PycharmProjects/DockStream-master/result/5T35'
    cmd = '/data/baiqing/PycharmProjects/DockStream-master/DockStreamCommunity/AutoDock-Vina-1.2.3/build/linux/release/vina --receptor 5t35_apo.pdbqt --ligand 0000.pdbqt --center_x 21.89 --center_y -28.24 --center_z -14.53 --size_x 40 --size_y 40 --size_z 40 --out 0000_ADV_out.pdbqt --seed 42 --num_modes 2'
    vina_docking(work_dir_5T35, cmd)

    work_dir_6W7O = '/data/baiqing/PycharmProjects/DockStream-master/result/6W7O'
    cmd = '/data/baiqing/PycharmProjects/DockStream-master/DockStreamCommunity/AutoDock-Vina-1.2.3/build/linux/release/vina --receptor 6w7o_apo.pdbqt --ligand 0000.pdbqt --center_x 22.38 --center_y 42.21 --center_z -5.93 --size_x 40 --size_y 40 --size_z 40 --out 0000_ADV_out.pdbqt --seed 42 --num_modes 2'
    vina_docking(work_dir_6W7O, cmd)

    work_dir_6ZHC = '/data/baiqing/PycharmProjects/DockStream-master/result/6ZHC'
    cmd = '/data/baiqing/PycharmProjects/DockStream-master/DockStreamCommunity/AutoDock-Vina-1.2.3/build/linux/release/vina --receptor 6zhc_apo.pdbqt --ligand 0000.pdbqt --center_x 27.59 --center_y 51.48 --center_z 17.79 --size_x 40 --size_y 40 --size_z 40 --out 0000_ADV_out.pdbqt --seed 42 --num_modes 2'
    vina_docking(work_dir_6ZHC, cmd)

    work_dir_7JTP = '/data/baiqing/PycharmProjects/DockStream-master/result/7JTP'
    cmd = '/data/baiqing/PycharmProjects/DockStream-master/DockStreamCommunity/AutoDock-Vina-1.2.3/build/linux/release/vina --receptor 7jtp_apo.pdbqt --ligand 0000.pdbqt --center_x 7.62 --center_y 7.19 --center_z -10.22 --size_x 40 --size_y 40 --size_z 40 --out 0000_ADV_out.pdbqt --seed 42 --num_modes 2'
    vina_docking(work_dir_7JTP, cmd)

    work_dir_7Q2J = '/data/baiqing/PycharmProjects/DockStream-master/result/7Q2J'
    cmd = '/data/baiqing/PycharmProjects/DockStream-master/DockStreamCommunity/AutoDock-Vina-1.2.3/build/linux/release/vina --receptor 7q2j_apo.pdbqt --ligand 0000.pdbqt --center_x -10.98 --center_y 6.47 --center_z 17.65 --size_x 40 --size_y 40 --size_z 40 --out 0000_ADV_out.pdbqt --seed 42 --num_modes 2'
    vina_docking(work_dir_7Q2J, cmd)

    work_dir_6BN7 = '/data/baiqing/PycharmProjects/DockStream-master/result/6BN7'
    cmd = '/data/baiqing/PycharmProjects/DockStream-master/DockStreamCommunity/AutoDock-Vina-1.2.3/build/linux/release/vina --receptor 6bn7_apo.pdbqt --ligand 0000.pdbqt --center_x 68.22 --center_y 40.96 --center_z 53.03 --size_x 40 --size_y 40 --size_z 40 --out 0000_ADV_out.pdbqt --seed 42 --num_modes 2'
    vina_docking(work_dir_6BN7, cmd)

    work_dir_7KHH = '/data/baiqing/PycharmProjects/DockStream-master/result/7KHH'
    cmd = '/data/baiqing/PycharmProjects/DockStream-master/DockStreamCommunity/AutoDock-Vina-1.2.3/build/linux/release/vina --receptor 7khh_apo.pdbqt --ligand 0000.pdbqt --center_x -11.48 --center_y 8.64 --center_z -16.94 --size_x 40 --size_y 40 --size_z 40 --out 0000_ADV_out.pdbqt --seed 42 --num_modes 2'
    vina_docking(work_dir_7KHH, cmd)

    work_dir_6BOY = '/data/baiqing/PycharmProjects/DockStream-master/result/6BOY'
    cmd = '/data/baiqing/PycharmProjects/DockStream-master/DockStreamCommunity/AutoDock-Vina-1.2.3/build/linux/release/vina --receptor 6boy_apo.pdbqt --ligand 0000.pdbqt --center_x 72.03 --center_y 38.67 --center_z 51.08 --size_x 40 --size_y 40 --size_z 40 --out 0000_ADV_out.pdbqt --seed 42 --num_modes 2'
    vina_docking(work_dir_6BOY, cmd)

    work_dir_6HR2 = '/data/baiqing/PycharmProjects/DockStream-master/result/6HR2'
    cmd = '/data/baiqing/PycharmProjects/DockStream-master/DockStreamCommunity/AutoDock-Vina-1.2.3/build/linux/release/vina --receptor 6hr2_apo.pdbqt --ligand 0000.pdbqt --center_x -27.68 --center_y -9.99 --center_z -22.06 --size_x 40 --size_y 40 --size_z 40 --out 0000_ADV_out.pdbqt --seed 42 --num_modes 2'
    vina_docking(work_dir_6HR2, cmd)

    work_dir_6HAY = '/data/baiqing/PycharmProjects/DockStream-master/result/6HAY'
    cmd = '/data/baiqing/PycharmProjects/DockStream-master/DockStreamCommunity/AutoDock-Vina-1.2.3/build/linux/release/vina --receptor 6hay_apo.pdbqt --ligand 0000.pdbqt --center_x -31.9 --center_y 13.95 --center_z -28.09 --size_x 40 --size_y 40 --size_z 40 --out 0000_ADV_out.pdbqt --seed 42 --num_modes 2'
    vina_docking(work_dir_6HAY, cmd)


    work_dir_7JTO = '/data/baiqing/PycharmProjects/DockStream-master/result/7JTO'
    cmd = '/data/baiqing/PycharmProjects/DockStream-master/DockStreamCommunity/AutoDock-Vina-1.2.3/build/linux/release/vina --receptor 7jto_apo.pdbqt --ligand 0000.pdbqt --center_x -10.15 --center_y 3.09 --center_z 15.65 --size_x 40 --size_y 40 --size_z 40 --out 0000_ADV_out.pdbqt --seed 42 --num_modes 2'
    vina_docking(work_dir_7JTO, cmd)































