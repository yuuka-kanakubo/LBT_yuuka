import numpy as np
import matplotlib.pyplot as plt
import os

from scipy.special import zeta

AnalyticalShow = True

C_2g = 3
C_2q = 4/3
c_g = 5.6
c_q = 5.8
alphas=0.3
N_C = 3
N_F = 3
C_A = 3
t_ti = 1.0/0.1973 #GEV^-1


def DebyeMass2(T):
	return (N_C + N_F/2)*4*np.pi*alphas*T*T/3

def qhat_g(E, T):
	"""
	Eq.(25) in arXiv: 2306.13742
        """ 
	return (C_2g/np.pi)*42.0*zeta(3)*alphas*alphas*T*T*T*np.log(c_g*E*T/(4.0*DebyeMass2(T)))

def qhat_q(E, T):
	""" 
	Eq.(25) in arXiv: 2306.13742
	""" 
	return (C_2q/np.pi)*42.0*zeta(3)*alphas*alphas*T*T*T*np.log(c_q*E*T/(4.0*DebyeMass2(T)))

def Gamma_inel_q(E, T):
	above_4 = (t_ti * (alphas * C_A * qhat_q(E, T)) / (2.0 * np.sqrt(DebyeMass2(T)))) * \
		(1.0 + (np.sqrt(DebyeMass2(T))/E) * np.log(np.sqrt(DebyeMass2(T))/E))
	below_4 = (t_ti * t_ti * (alphas * C_A * qhat_q(E, T)) / (4.0 * np.pi)) * \
		np.log(E/np.sqrt(DebyeMass2(T)))
	# Use np.where to choose
	return np.where(E > 4, above_4/0.1973, below_4/0.1973)

def Gamma_inel_g(E, T):
	above_4 = (t_ti * (alphas * C_A * qhat_g(E, T)) / (2.0 * np.sqrt(DebyeMass2(T)))) * \
		(1.0 + (np.sqrt(DebyeMass2(T))/E) * np.log(np.sqrt(DebyeMass2(T))/E))
	below_4 = (t_ti * t_ti * (alphas * C_A * qhat_g(E, T)) / (4.0 * np.pi)) * \
		np.log(E/np.sqrt(DebyeMass2(T)))
	# Use np.where to choose
	return np.where(E > 4, above_4/0.1973, below_4/0.1973)

# List your data files and their plot settings
data_files = [
    {"filename": "DATA_Gamma_inel_q.dat", "xlabel": "$E_0$", "ylabel": "$\Gamma_{inel} \ \mathrm{[fm^{-1}]}$", "title": "Inelastic scat. rate of quarks",  "func": Gamma_inel_q},
    {"filename": "DATA_qhat_q.dat",       "xlabel": "$E_0$", "ylabel": "$\hat{q} \ \mathrm{[GeV^{3}}]$",      "title": "qhat of quarks",                  "func": qhat_q},
    {"filename": "DATA_Gamma_inel_g.dat", "xlabel": "$E_0$", "ylabel": "$\Gamma_{inel} \ \mathrm{[fm^{-1}]}$", "title": "Inelastic scat. rate of gluons",  "func": Gamma_inel_g},
    {"filename": "DATA_qhat_g.dat",       "xlabel": "$E_0$", "ylabel": "$\hat{q} \ \mathrm{[GeV^{3}]}$",      "title": "qhat of gluons",                  "func": qhat_g},
]


if AnalyticalShow:
	# Vectorize the function
	Gamma_inel_q_vec = np.vectorize(Gamma_inel_q)
	
	# Create E and T mesh
	E_vals = np.linspace(1.0, 200.0, 100)
	T_vals = np.linspace(0.1, 0.6, 100)
	E_grid, T_grid = np.meshgrid(E_vals, T_vals)
	
	# Compute Gamma_inel_q on the grid
	Gamma_grid = Gamma_inel_q_vec(E_grid, T_grid)
	
	# Plotting
	fig = plt.figure(figsize=(10, 7))
	ax = fig.add_subplot(111, projection='3d')
	surf = ax.plot_surface(E_grid, T_grid, Gamma_grid, cmap='viridis', edgecolor='none')
	
	ax.set_xlabel('E')
	ax.set_ylabel('T')
	ax.set_zlabel(r'$\Gamma_{\mathrm{inel}}^q(E, T) \mathrm{fm^{-1}}$')
	ax.set_title(r'3D Plot of $\Gamma_{\mathrm{inel}}^q(E, T)$')
	fig.colorbar(surf, shrink=0.5, aspect=5)
	ax.invert_xaxis()
	ax.view_init(elev=30, azim=45)
	
	plt.savefig('AnalyticalGammaInel_q.pdf')
	plt.clf()
	
	
	# Vectorize the function
	Gamma_inel_g_vec = np.vectorize(Gamma_inel_g)
	
	# Create E and T mesh
	E_vals = np.linspace(1.0, 200.0, 100)
	T_vals = np.linspace(0.1, 0.6, 100)
	E_grid, T_grid = np.meshgrid(E_vals, T_vals)
	
	# Compute Gamma_inel_g on the grid
	Gamma_grid = Gamma_inel_g_vec(E_grid, T_grid)
	
	# Plotting
	fig = plt.figure(figsize=(10, 7))
	ax = fig.add_subplot(111, projection='3d')
	surf = ax.plot_surface(E_grid, T_grid, Gamma_grid, cmap='viridis', edgecolor='none')
	
	ax.set_xlabel('E')
	ax.set_ylabel('T')
	ax.set_zlabel(r'$\Gamma_{\mathrm{inel}}^g(E, T) \mathrm{fm^{-1}}$')
	ax.set_title(r'3D Plot of $\Gamma_{\mathrm{inel}}^g(E, T)$')
	fig.colorbar(surf, shrink=0.5, aspect=5)
	ax.invert_xaxis()
	ax.view_init(elev=30, azim=45)
	
	plt.savefig('AnalyticalGammaInel_g.pdf')
	plt.clf()

##
n_files = len(data_files)

fig, axes = plt.subplots(n_files, 1, figsize=(6, 4*n_files))  # One column, n_files rows
T_fixed = 0.3

if n_files == 1:
    axes = [axes]  # make iterable

for ax, data_info in zip(axes, data_files):
    data = np.loadtxt(data_info["filename"])
    x = data[:,0]
    y = data[:,1]
    ax.plot(x, y, marker='o', linestyle='-', ms=1.5, color="goldenrod", alpha=0.5, label="LBT")
    ax.plot(x, [data_info["func"](x_, T_fixed) for x_ in x], linestyle='--', lw=1, color="black", alpha=0.75, label="analytical")
    ax.set_xlabel(data_info["xlabel"])
    ax.set_ylabel(data_info["ylabel"])
    ax.set_title(data_info["title"])
    ax.grid(True)
    ax.legend(fontsize='small', loc='best')

plt.tight_layout()
plt.savefig('TestPlots.pdf')
plt.close()
