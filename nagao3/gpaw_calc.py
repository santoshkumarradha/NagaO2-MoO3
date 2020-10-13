from gpaw import GPAW, PW
import pickle
import numpy as np
from ase.io import read
from ase.parallel import parprint


def calc(fname="POSCAR_li_p"):
    vol = []
    e = []
    for x in np.linspace(.95, 1.05, 6):
        parprint(f"GGA {fname} @ {x}")
        parprint("="*50)
        a = read(fname)
        pos = a.get_scaled_positions()
        a.set_cell(a.cell*x)
        a.set_scaled_positions(pos)
        calc = GPAW(mode=PW(800),
                    xc='PBE',
                    kpts=(15, 15, 15),
                    basis='dzp',
                    txt=f'pn21-{x}.txt')
        a.calc = calc
        e.append(a.get_potential_energy())
        vol.append(a.get_volume())

    data = [vol, e]
    with open(f'data_{fname}.pickle', 'wb') as handle:
        pickle.dump(data, handle, protocol=pickle.HIGHEST_PROTOCOL)
    parprint("done")
    parprint("-"*50)


calc(fname="POSCAR_li_p")
calc(fname="POSCAR_li_r")
