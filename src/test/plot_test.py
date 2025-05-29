import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Constants
alpha_s = 0.3
CA = 3.0
t_minus_ti = 1.0/0.197  # [fm/c]
pi = np.pi
zeta_3 = 1.202
C2_g = 3.0
c_a = 5.6

# Updated mu_D
def mu_D(T):
    return np.sqrt((3.0/2.0) * 4.0 * pi * alpha_s * T**2)

# Updated qhat_a
def qhat_a(E, T):
    muD = mu_D(T)
    ln_arg = c_a * E * T / (muD**2)
    return (C2_g * (42 * zeta_3) / pi) * alpha_s**2 * T**3 * np.log(ln_arg)

# Inelastic rate function
def Gamma_inel(E, T):
    muD = mu_D(T)
    qhat = qhat_a(E, T)
    factor = (t_minus_ti * alpha_s * CA * qhat) / (2 * muD)
    term2 = (muD / E) * np.log(muD / E)
    return factor * (1 + term2)

# Create E and T grid
E_vals = np.linspace(1, 200, 100)  # E in GeV
T_vals = np.linspace(0.2, 0.5, 100)  # T in GeV
E_grid, T_grid = np.meshgrid(E_vals, T_vals)

# Compute Gamma_inel
Gamma_grid = Gamma_inel(E_grid, T_grid)

# Plotting
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')

surf = ax.plot_surface(T_grid, E_grid, Gamma_grid, cmap='coolwarm', edgecolor='none')

ax.set_xlabel('T (GeV)')
ax.set_ylabel('E (GeV)')
ax.set_zlabel(r'$\Gamma_{\mathrm{inel}}$')
ax.set_title(r'Inelastic Rate $\Gamma_{\mathrm{inel}}(E, T)$')

# Set viewing angle similar to FIG.5
ax.view_init(elev=30, azim=-45)

# Color bar
fig.colorbar(surf, shrink=0.5, aspect=5, label=r'$\Gamma_{\mathrm{inel}}$')

plt.tight_layout()
plt.savefig('test.pdf')

