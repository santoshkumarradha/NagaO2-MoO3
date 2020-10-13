from lmf import lmf  # load the lmf calculator
import numpy as np
from ase.io import read
from ase.eos import calculate_eos
from ase.io.trajectory import Trajectory
import matplotlib.pyplot as plt
from ase.eos import EquationOfState as eos
from ase.units import Bohr, Rydberg
from ase.io.trajectory import TrajectoryReader as trread


def plot_gga(ax, ax1):
    def murnaghan(V, E0, B0, BP, V0):
        'From PRB 28,5480 (1983'

        E = E0 + B0 * V / BP * \
            (((V0 / V)**BP) / (BP - 1) + 1) - V0 * B0 / (BP - 1)
        return E
    traj = Trajectory('nagao2-r3m-gga.traj')
    v = [i.get_volume() for i in traj]
    e1 = np.array([i.get_total_energy() * Rydberg for i in traj])
    eos_data = eos(v, e1, eos='murnaghan')
    eos_data.fit()
    V = np.linspace(35, 54, 100)
    E = murnaghan(V, eos_data.eos_parameters[0], eos_data.eos_parameters[1],
                  eos_data.eos_parameters[2], eos_data.eos_parameters[3])
    r3m_e = E
    r3m_v = V

    traj = Trajectory('nagao2-p2n1_gga.traj')
    v = [i.get_volume()*0.25 for i in traj]
    e1 = np.array([i.get_total_energy() * Rydberg*0.25 for i in traj])
    eos_data = eos(v, e1, eos='murnaghan')
    eos_data.fit()
    V = np.linspace(45, 65, 100)
    E = murnaghan(V, eos_data.eos_parameters[0], eos_data.eos_parameters[1],
                  eos_data.eos_parameters[2], eos_data.eos_parameters[3])
    p2n1_e = E
    p2n1_v = V

    minimum = np.min([r3m_e, p2n1_e])
    c = "k"
    label = "R$\\bar{3}$m"
    ax.plot(r3m_v, r3m_e-minimum, c=c, label=label)
    c = "r"
    label = "Pna$2_1$"
    ax.plot(p2n1_v, p2n1_e-minimum, c=c, label=label)
    ax.legend()

    c = "k"
    label = "R$\\bar{3}$m"
    p = -1*np.gradient(r3m_e, np.diff(r3m_v)[0])
    h = r3m_e-p*r3m_v
    h1 = h
    p1 = p
    ax1.plot(p, h, c=c, label=label)

    c = "r"
    label = "Pna$2_1$"
    p = -1*np.gradient(p2n1_e, np.diff(p2n1_v)[0])
    h = p2n1_e-p*p2n1_v
    h2 = h
    p2 = p
    ax1.plot(p, h, c=c, label=label)
    ax1.legend()
    return h1, h2, p1, p2


fig, ax = plt.subplots(1, 2)
h1, h2, p1, p2 = plot_gga(ax[0], ax[1])
plt.tight_layout()
