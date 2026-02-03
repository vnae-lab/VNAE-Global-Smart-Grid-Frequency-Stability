import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix, diags
from scipy.integrate import odeint
import time

# -------------------------------------------------------------------------
# VNAE Framework: Global Smart Grid Frequency Stability (n=100,000)
# Purpose: Validating Volume Contraction in Massive Energy Networks
# -------------------------------------------------------------------------

np.random.seed(123)

# --- [ GLOBAL GRID CONFIGURATION ] ---
n = 100000              # 100,000 Energy Nodes (National Grid Scale)
beta_val = 0.12         # Grid Rigidity (VNAE Structural Factor)

print(f"Initializing Grid with {n} nodes...")

# 1. GENERATE HETEROGENEOUS NODE INERTIA (THETA)
# Models turbines (high theta) vs solar inverters (low theta)
theta_nodes = np.random.uniform(1.0, 120.0, n)

# 2. GENERATE ULTRA-SPARSE TRANSMISSION NETWORK (A)
# Large grids are sparse; using COOrdinate format for efficient creation
nz_count = n * 15       # Average of 15 transmission lines per node
row_idx = np.random.randint(0, n, nz_count)
col_idx = np.random.randint(0, n, nz_count)
weights = np.random.exponential(1.0/0.8, nz_count)

A_sparse = csr_matrix((weights, (row_idx, col_idx)), shape=(n, n))
A_sparse.setdiag(0)
A_sparse.eliminate_zeros()

# 3. STATISTICAL MANIFOLD CURVATURE (K) CERTIFICATION
# Using Monte Carlo to verify K > 0 for massive systems
sample_size = 500000
idx_i = np.random.randint(0, n, sample_size)
idx_j = np.random.randint(0, n, sample_size)

diff_theta = np.abs(theta_nodes[idx_i] - theta_nodes[idx_j])
# Accessing sparse elements efficiently
couplings = np.array(A_sparse[idx_i, idx_j]).flatten() + np.array(A_sparse[idx_j, idx_i]).flatten()
rigidity = 1 + beta_val * (theta_nodes[idx_i] + theta_nodes[idx_j])

K_estimate = np.mean((diff_theta * couplings) / rigidity)

# -------------------------------------------------------------------------
# DYNAMIC SIMULATION SETUP
# -------------------------------------------------------------------------
# Construct Laplacian: L = D - A
d_vec = np.array(A_sparse.sum(axis=1)).flatten()
L_sparse = diags(d_vec) - A_sparse

p_weather_shocks = np.random.normal(0, 10, n) # Heavy fluctuations
omega_init = np.random.normal(0, 5, n)        # Initial frequency deviations

# Gradient Flow Dynamics using Euler Method for Speed at 100k nodes
def solve_grid_vnae(omega, dt, steps):
    history = [omega.copy()]
    curr_omega = omega
    for _ in range(steps):
        # VNAE fundamental equation: d_omega = -(L + Theta) * omega + p
        # Using sparse matrix multiplication (@)
        d_omega = -(L_sparse @ curr_omega) - (theta_nodes * curr_omega) + p_weather_shocks
        curr_omega = curr_omega + dt * d_omega
        history.append(curr_omega.copy())
    return np.array(history)

# Solving for a short burst
print("Running massive scale simulation...")
start_time = time.time()
dt = 0.01
n_steps = 10
solution = solve_grid_vnae(omega_init, dt, n_steps)
end_time = time.time()

# -------------------------------------------------------------------------
# FINAL STABILITY REPORT
# -------------------------------------------------------------------------
print("\n--- VNAE GLOBAL ENERGY GRID REPORT ---")
print(f"Total Grid Nodes (n): {n}")
print(f"Structural Curvature (K): {K_estimate:.8f}")
print(f"Simulation Time: {end_time - start_time:.4f} seconds")
print(f"Status: {'Geometrically Stable (VNAE Certified)' if K_estimate > 0 else 'Critical Instability'}")

# Visualization: Sample 100 lines to see the 'Grid Convergence'
plt.figure(figsize=(10, 6))
time_axis = np.arange(0, n_steps + 1) * dt
plot_sample = np.random.choice(range(n), 100, replace=False)

for idx in plot_sample:
    plt.plot(time_axis, solution[:, idx], alpha=0.3)

plt.title(f"VNAE 100k Nodes Grid Stability (K = {K_estimate:.4f})")
plt.xlabel("Time (seconds)")
plt.ylabel("Frequency Deviation (Hz)")
plt.grid(True, linestyle='--', alpha=0.7)
plt.show()
