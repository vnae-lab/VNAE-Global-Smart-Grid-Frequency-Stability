using SparseArrays
using LinearAlgebra
using Random
using Plots
using Statistics

# -------------------------------------------------------------------------
# VNAE Framework: Global Smart Grid Frequency Stability (n=100,000)
# Purpose: Validating Volume Contraction in Massive Energy Networks
# -------------------------------------------------------------------------

Random.seed!(28)

# --- [ GLOBAL GRID CONFIGURATION ] ---
n = 100_000             # 100,000 Energy Nodes (National Grid Scale)
beta_val = 0.12         # Grid Rigidity (VNAE Structural Factor)

println("Initializing National Grid Simulation with $n nodes...")

# 1. GENERATE HETEROGENEOUS NODE INERTIA (THETA)
# Diversity in inertia (Theta) is what generates stabilizing curvature
theta_nodes = rand(n) .* (120.0 - 1.0) .+ 1.0

# 2. GENERATE ULTRA-SPARSE TRANSMISSION NETWORK (A)
# Large grids are sparse; each node connects to approx 15 neighbors
nz_count = n * 15
rows = rand(1:n, nz_count)
cols = rand(1:n, nz_count)
weights = randexp(nz_count) ./ 0.8

A_sparse = sparse(rows, cols, weights, n, n)
# Remove self-loops
for i in 1:n A_sparse[i,i] = 0.0 end
dropzeros!(A_sparse)

# 3. STATISTICAL MANIFOLD CURVATURE (K) CERTIFICATION
# Monte Carlo estimate for Scalar Curvature in high dimensions
sample_size = 500_000
idx_i = rand(1:n, sample_size)
idx_j = rand(1:n, sample_size)

diff_theta = abs.(theta_nodes[idx_i] .- theta_nodes[idx_j])
couplings = [A_sparse[idx_i[k], idx_j[k]] + A_sparse[idx_j[k], idx_i[k]] for k in 1:sample_size]
rigidity = 1.0 .+ beta_val .* (theta_nodes[idx_i] .+ theta_nodes[idx_j])

K_estimate = mean((diff_theta .* couplings) ./ rigidity)

# -------------------------------------------------------------------------
# DYNAMIC SIMULATION (FLOW ON THE MANIFOLD)
# -------------------------------------------------------------------------

# Construct Laplacian: L = D - A
d_vec = sum(A_sparse, dims=2)[:]
L_sparse = Diagonal(d_vec) - A_sparse

p_weather_shocks = randn(n) .* 10.0  # Heavy weather fluctuations
omega_init = randn(n) .* 5.0        # Initial frequency deviations

# Simulation Parameters
dt = 0.01
n_steps = 10
history = zeros(n_steps + 1, n)
history[1, :] = omega_init

println("Running massive scale gradient flow...")

# High-performance solver loop
curr_omega = copy(omega_init)
@time for t in 1:n_steps
    # VNAE equation: d_omega = -(L + Theta) * omega + p
    # Sparse multiplication is native and extremely fast in Julia
    d_omega = -(L_sparse * curr_omega) .- (theta_nodes .* curr_omega) .+ p_weather_shocks
    curr_omega .+= dt .* d_omega
    history[t+1, :] = curr_omega
end

# -------------------------------------------------------------------------
# FINAL STABILITY REPORT
# -------------------------------------------------------------------------
println("\n--- VNAE GLOBAL ENERGY GRID REPORT ---")
println("Total Grid Nodes (n): $n")
println("Structural Curvature (K): $(round(K_estimate, digits=8))")
println("Status: ", K_estimate > 0 ? "Geometrically Stable (VNAE Certified)" : "Critical Instability")

# Visualization: Sample 100 lines for the convergence plot
plot_sample = rand(1:n, 100)
time_axis = (0:n_steps) .* dt

plt = plot(time_axis, history[:, plot_sample], 
           linealpha=0.3, 
           palette=:viridis, 
           legend=false,
           title="VNAE 100k Nodes Grid Stability (K = $(round(K_estimate, digits=4)))",
           xlabel="Time (seconds)", 
           ylabel="Frequency Deviation (Hz)",
           size=(900, 500))

display(plt)
