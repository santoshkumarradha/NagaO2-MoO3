import pickle
from lmf import lmf  # load the lmf calculator
from ase.io import read
from ase.eos import calculate_eos


print("Calculating Pn21 GGA \n =========================")
atoms = read("POSCAR_p2n1")
atoms.calc = lmf(nkabc=[4, 4, 4], gmax=12, p=12,
                 pbe=True)  # p=num of processors
eos_p2n1_gga = calculate_eos(atoms, trajectory='nagao2-p2n1_gga.traj')
eos_p2n1_gga.plot('nagao2-p2n1-eos_gga.png')

print("Calculating Pn21 LDA \n =========================")
atoms.calc = lmf(nkabc=[4, 4, 4], gmax=12, p=12,
                 pbe=False)  # p=num of processors
eos_p2n1_lda = calculate_eos(atoms, trajectory='nagao2-p2n1_lda.traj')
eos_p2n1_lda.plot('nagao2-p2n1-eos_lda.png')

print("Calculating r3m GGA \n =========================")
atoms = read("POSCAR")
atoms.calc = lmf(nkabc=[6, 6, 6], gmax=12, p=12,
                 pbe=True)  # p=num of processors
eos_r3m_gga = calculate_eos(atoms, trajectory='nagao2-r3m_gga.traj')
eos_r3m_gga.plot('nagao2-r3m-eos_gga.png')

print("Calculating r3m LDA \n =========================")
atoms.calc = lmf(nkabc=[6, 6, 6], gmax=12, p=12,
                 pbe=False)  # p=num of processors
eos_r3m_lda = calculate_eos(atoms, trajectory='nagao2-r3m_lda.traj')
eos_r3m_lda.plot('nagao2-r3m-eos_lda.png')

print("daving calc \n =========================")
data = [eos_p2n1_gga, eos_p2n1_lda, eos_r3m_gga, eos_r3m_lda]
with open('eos_data.pickle', 'wb') as handle:
    pickle.dump(a, handle, protocol=pickle.HIGHEST_PROTOCOL)
print("done :) \n =========================")
