import dill as pickle
from ase.build import bulk
from gpaw import GPAW, FermiDirac
from ase.phonons import Phonons
import ase
import numpy as np
atoms = ase.io.read("moo3_bulk.cif")

calc = GPAW(symmetry={'point_group': False},
            mode='lcao',
            basis='szp(dzp)',
            kpts=(6, 6, 3),
            convergence={'density': 1e-7},
            xc='PBE',  # No PBEsol :(
            occupations=FermiDirac(0.01))

ph = Phonons(atoms, calc, supercell=(3, 3, 2), delta=0.05)
ph.run()

ph.read(acoustic=True)
force_constants = ph.get_force_constant()
ph.acoustic(force_constants)
force_constants = ph.symmetrize(force_constants)
np.save("force_constant.npy", force_constants)
with open('phonons.pkl', 'wb') as file:
    pickle.dump(ph, file)
