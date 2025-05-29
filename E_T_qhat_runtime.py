import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

# Load your data
data = np.loadtxt('E_T_qhat.dat')  # replace with your filename

E0 = data[:, 0]
T = data[:, 1]
qhat = data[:, 2]

fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(E0, T, qhat, c=qhat, cmap='viridis', marker='o')
ax.set_xlabel('E₀')
ax.set_ylabel('T')
ax.set_zlabel('q̂')
ax.set_title('3D plot of q̂(E₀, T)')
plt.tight_layout()
plt.savefig('qhat_in3D.pdf')
plt.show()

# Choose a tolerance for T (since floats might not match exactly)
T_slice = 0.3
tol = 1e-4
mask = np.abs(T - T_slice) < tol

E0_slice = E0[mask]
qhat_slice = qhat[mask]

plt.figure(figsize=(7,5))
plt.plot(E0_slice, qhat_slice, marker='o', linestyle='-', color='blue', label=f'T = {T_slice}')
plt.xlabel('E₀')
plt.ylabel('q̂')
plt.title(f'q̂ vs E₀ at T={T_slice}')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig('qhat_at_T0.3.pdf')
plt.show()
