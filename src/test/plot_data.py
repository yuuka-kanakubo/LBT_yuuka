import numpy as np
import matplotlib.pyplot as plt
import os

from scipy.special import zeta

DebyeMass2 = 0.5089380098815464**0.5

#T_med = 0.3 is assumed
def qhat_q(E):
    return (4/(3*np.pi))*42.0*zeta(3)*0.3*0.3*0.3*0.3*0.3*np.log(5.8*E*0.3/DebyeMass2)

def Gamma_inel_q(E):
    return (1.0 * (0.3 * 3.0 * qhat_q(E)) / (4.0*np.pi)) * np.log(E/np.sqrt(DebyeMass2))

# List your data files and their plot settings
data_files = [
    {"filename": "DATA_Gamma_inel_q.dat", "xlabel": "$E_0$", "ylabel": "$Gamma_{inel}$", "title": "",  "func": Gamma_inel_q},
    {"filename": "DATA_qhat_q.dat",       "xlabel": "$E_0$", "ylabel": "$\hat{q}$",      "title": "",  "func": qhat_q},
]

n_files = len(data_files)

fig, axes = plt.subplots(n_files, 1, figsize=(6, 4*n_files))  # One column, n_files rows

if n_files == 1:
    axes = [axes]  # make iterable

for ax, data_info in zip(axes, data_files):
    data = np.loadtxt(data_info["filename"])
    x = data[:,0]
    y = data[:,1]
    ax.plot(x, [data_info["func"](x_) for x_ in x], marker='x', linestyle='-', ms=1.5, color="gray", alpha=0.5, label="semi-analytical")
    ax.plot(x, y, marker='o', linestyle='-', ms=1.5, color="goldenrod", alpha=0.5, label="LBT")
    ax.set_xlabel(data_info["xlabel"])
    ax.set_ylabel(data_info["ylabel"])
    ax.set_title(data_info["title"])
    ax.grid(True)
    ax.legend(fontsize='small', loc='best')

plt.tight_layout()
plt.savefig('TestPlots.pdf')
plt.close()
