from pylab import *
# from deepmd.calculator import DP
from ase.io import read, write
from tqdm import tqdm

# calc = DP(model="graph-compress.pb")
# def dp_dft(file_path, label):
#     atoms = read(file_path, ":")
#     dft_energy = []
#     dp_energy = []
#     dft_force = []
#     dp_force = []
#     dft_virial = []
#     dp_virial = []
#     for i in tqdm(range(len(atoms))):
#         atom = atoms[i]
#         atom.calc = calc
#         dft_energy.append(atom.info["energy"] / len(atom.numbers)) #unit in eV/atom
#         dp_energy.append(atom.get_potential_energy() / len(atom.numbers)) #unit in eV/atom
#         dft_force.append(atom.arrays["forces"].reshape(-1)) #unit in eV/A
#         dp_force.append(atom.get_forces().reshape(-1)) #unit in eV/A
#         if "virial" in atom.info.keys():
#             virial = atom.info["virial"] #unit in eV
#             virial = np.array([virial[0, 0], virial[1, 1], virial[2, 2], virial[1, 2], virial[0, 2], virial[0, 1]])  #xx yy zz yz xz xy
#             virial = virial / len(atom.numbers) #unit in eV/atom
#             dft_virial.append(virial)
#             stress = atom.get_stress() #unit in eV/A3
#             stress = stress / len(atom.numbers) #unit in eV/atom/A3
#             dp_virial.append(stress * atom.get_volume()*(-1)) #unit in eV/atom
#         else:
#             pass


#     dft_energy = np.array(dft_energy)
#     dp_energy = np.array(dp_energy)
#     dft_force = np.concatenate(dft_force)
#     dp_force = np.concatenate(dp_force)
#     dft_virial = np.concatenate(dft_virial)
#     dp_virial = np.concatenate(dp_virial)
#     output_energy = np.c_[dft_energy, dp_energy]
#     output_force  = np.c_[dft_force, dp_force]
#     output_virial = np.c_[dft_virial, dp_virial]
#     np.save(f"{label}_energy.npy", output_energy)
#     np.save(f"{label}_force.npy", output_force)
#     return output_energy, output_force, output_virial

# energy_train, force_train, virial_train = dp_dft("./train.xyz", "train")
# energy_test, force_test, virial_test = dp_dft("./Reference_Sampling.xyz", "test")

energy_train = np.load("train_energy.npy")
energy_test = np.load("test_energy.npy")
force_train = np.load("train_force.npy")
force_test = np.load("test_force.npy")

##set figure properties
aw = 1.5
fs = 12
lw = 1.5
font = {'size'   : fs}
matplotlib.rc('font', **font)
matplotlib.rc('axes' , lw=aw)

def set_fig_properties(ax_list):
    tl = 6
    tw = 1.5
    tlm = 3
    
    for ax in ax_list:
        ax.tick_params(which='major', length=tl, width=tw)
        ax.tick_params(which='minor', length=tlm, width=tw)
        ax.tick_params(which='both', axis='both', direction='out', right=False, top=False)


figure(figsize=(12, 5))
subplot(1, 2, 1)
set_fig_properties([gca()])
plot(energy_train[:, 0], energy_train[:, 1], 'o', c="C0", ms = 5, alpha=0.5, label="Train")
e_train_rmse = np.sqrt(sum((energy_train[:,0] - energy_train[:,1]) ** 2) / len(energy_train[:,0])) * 1000
plot(energy_test[:, 0], energy_test[:, 1], 'o', c="C1", ms = 5, alpha=0.5, label="Test")
e_test_rmse = np.sqrt(sum((energy_test[:,0] - energy_test[:,1]) ** 2) / len(energy_test[:,0])) * 1000
begin = np.min(energy_train)-0.1
end = np.max(energy_train)+0.1
plot([begin, end], [begin, end], c = "grey", lw = 1)
xlim([begin, end])
ylim([begin, end])
delta = end - begin
text(begin + delta*0.5, begin + delta*0.2, f"{e_train_rmse:.1f} meV/atom", color="C0")
text(begin + delta*0.5, begin + delta*0.1, f"{e_test_rmse:.1f} meV/atom", color="C1")
xlabel('DFT energy (eV/atom)')
ylabel('DP energy (eV/atom)')
legend(loc="upper left")
title("(a)")


subplot(1, 2, 2)
set_fig_properties([gca()])
plot(force_train[:, 0], force_train[:, 1], 'o', c="C2", ms = 5, alpha=0.5, label="Train")
f_train_rmse = np.sqrt(sum((force_train[:,0] - force_train[:,1]) ** 2) / len(force_train[:,0])) * 1000
plot(force_test[:, 0], force_test[:, 1], 'o', c="C3", ms = 5, alpha=0.5, label="Test")
f_test_rmse = np.sqrt(sum((force_test[:,0] - force_test[:,1]) ** 2) / len(force_test[:,0])) * 1000
begin = np.min(force_train) - 1
end = np.max(force_train) + 1
plot([begin, end], [begin, end], c = "grey", lw = 1)
xlim([begin, end])
ylim([begin, end])
delta = end - begin
text(begin + delta*0.5, begin + delta*0.2, f"{f_train_rmse:.1f} meV/" + r"$\rm{\AA}$", color="C2")
text(begin + delta*0.5, begin + delta*0.1, f"{f_test_rmse:.1f} meV/" + r"$\rm{\AA}$", color="C3")
xlabel(r'DFT force (eV/$\rm{\AA}$)')
ylabel(r'DP force (eV/$\rm{\AA}$)')
legend(loc="upper left")
title("(b)")


# subplot(1, 3, 3)
# set_fig_properties([gca()])
# plot(virial_train[:, 0], virial_train[:, 1], 'o', c="C2", ms = 5, alpha=0.5, label="Train")
# v_train_rmse = np.sqrt(sum((virial_train[:,0] - virial_train[:,1]) ** 2) / len(virial_train[:,0])) * 1000
# begin = np.min(virial_train) - 0.1
# end = np.max(virial_train) + 0.1
# plot([begin, end], [begin, end], c = "grey", lw = 1)
# xlim([begin, end])
# ylim([begin, end])
# delta = end - begin
# text(begin + delta*0.5, begin + delta*0.2, f"{v_train_rmse:.1f} meV/atom", color="C2")
# xlabel('DFT virial (eV/atom)')
# ylabel('DP virial (eV/atom)')
# legend(loc="upper left")
# title("(c)")

subplots_adjust(wspace=0.3, hspace=0.3)
savefig("RMSE.png", bbox_inches='tight')
