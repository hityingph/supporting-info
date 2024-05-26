from pylab import *
from deepmd.calculator import DP
from ase.io import read, write
from tqdm import tqdm

calc = DP(model="./graph-compress.pb")
def dp_dft(file_path):
    atoms = read(file_path, ":")
    dft_energy = []
    dp_energy = []
    dft_force = []
    dp_force = []
    dft_virial = []
    dp_virial = []
    for i in tqdm(range(len(atoms))):
        atom = atoms[i]
        atom.calc = calc
        dft_energy.append(atom.info["energy"] / len(atom.numbers)) #unit in eV/atom
        dp_energy.append(atom.get_potential_energy() / len(atom.numbers)) #unit in eV/atom
        dft_force.append(atom.arrays["forces"].reshape(-1)) #unit in eV/A
        dp_force.append(atom.get_forces().reshape(-1)) #unit in eV/A
        if "virial" in atom.info.keys():
            virial = atom.info["virial"] #unit in eV
            virial = np.array([virial[0, 0], virial[1, 1], virial[2, 2], virial[1, 2], virial[0, 2], virial[0, 1]])  #xx yy zz yz xz xy
            virial = virial / len(atom.numbers) #unit in eV/atom
            dft_virial.append(virial)
            stress = atom.get_stress() #unit in eV/A3
            stress = stress / len(atom.numbers) #unit in eV/atom/A3
            dp_virial.append(stress * atom.get_volume()*(-1)) #unit in eV/atom
        else:
            pass


    dft_energy = np.array(dft_energy)
    dp_energy = np.array(dp_energy)
    dft_force = np.concatenate(dft_force)
    dp_force = np.concatenate(dp_force)
    dft_virial = np.concatenate(dft_virial)
    dp_virial = np.concatenate(dp_virial)
    output_energy = np.c_[dft_energy, dp_energy]
    output_force  = np.c_[dft_force, dp_force]
    output_virial = np.c_[dft_virial, dp_virial]
    # np.save("energy_test.npy", output_energy)
    np.save("force_test.npy", output_force)
    # np.save("virial_test.npy", output_virial)
    return output_energy, output_force, output_virial

energy_test, force_test, virial_test = dp_dft("./test_c-si.xyz")