# Global Smart Grid Frequency Stability

This repository presents a large-scale stability validation of power grid frequency dynamics using the Victoria-Nash Asymmetric Equilibrium (VNAE) and/or Victoria-Nash Geometry framework.

The model demonstrates how **global stability can emerge from asymmetric dissipation**, even in massive (with 100,000 agents), heterogeneous, and sparsely connected energy networks without requiring synchronization or eigenvalue-based tuning.

---

## Motivation

Modern power grids operate under:
- high renewable penetration,
- heterogeneous inertia (turbines, inverters, storage),
- sparse and directed transmission topologies,
- persistent external disturbances (weather, demand shocks).

Classical small-signal or symmetric consensus models struggle to certify global stability at scale.

We can say that the **VNAE provides a geometric alternative** in which stability is governed by **volume contraction on a curved manifold** rather than pointwise convergence.

---

## Core Model

The grid dynamics are modeled as a dissipative flow on a networked manifold:

dω/dt = − (L + Θ) · ω + p

where:

ω ∈ Rⁿ = Vector of node frequency deviations (one component per grid node).

L = Directed graph Laplacian encoding transmission-line couplings and network topology.

Θ = diag(θ₁, θ₂, …, θₙ)
Diagonal matrix of heterogeneous asymmetric dissipation parameters
(e.g. generator inertia, inverter damping, load responsiveness).

p = Persistent external forcing vector
(e.g. weather-induced renewable intermittency or demand shocks).

---

## Large-Scale Configuration

- **Grid size:** 100,000 nodes  
- **Topology:** ultra-sparse directed network  
- **Node heterogeneity:** wide inertia distribution  
- **Dynamics:** nonlinear dissipative flow  
- **Numerics:** sparse matrices + Monte Carlo geometry

This scale moves the framework well beyond toy models and into realistic infrastructure regimes.

### Network Sparsity Assumption

The transmission network is modeled as a sparse directed graph, where each node is connected on average to approximately **15 transmission lines**.  
This reflects realistic large-scale power grids, in which physical, economic, and geographical constraints lead to low-degree connectivity even at national or continental scales.  
The sparse structure ensures computational scalability while preserving the essential propagation pathways of frequency disturbances.


---

## Curvature Proxy (Scalable Stability Certification)

Global stability is assessed via a statistical curvature proxy:

K = E [ |θ_i − θ_j| · |A_ij| / ( 1 + β · (θ_i + θ_j) ) ]

Estimated via Monte Carlo sampling:

(i, j) ~ Uniform( {1,…,n} × {1,…,n} )


### Interpretation

- **\(K > 0\)** → positively curved effective manifold  
- **Positive curvature** → volume contraction  
- **Volume contraction** → global stability (VNAE criterion)

No eigenvalue alignment or synchronization is required.

---

## Simulation Results

**Report excerpt:**

- Total Grid Nodes (n): 100,000
- Structural Curvature (K): 0.00121638
- Metric Determinant Proxy: 1.4e+120
- Status: Geometrically Stable (VNAE Certified)


### Plot Interpretation

- Only a **random subset of nodes (n = 100)** is visualized for readability.
- All trajectories rapidly contract into a narrow band.
- Late-time dispersion reflects:
  - heterogeneous dissipation,
  - sparse directed coupling,
  - numerical integration noise.

**Importantly:**  
Pointwise convergence is not required. Volume contraction is the stability mechanism.

---

## Why This Matters

- Scales to **national-grid-sized systems**
- Robust to heterogeneity and asymmetry
- Avoids fragile spectral assumptions
- Applicable to:
  - smart grids
  - inverter-dominated systems
  - systemic risk networks
  - large-scale cyber-physical infrastructure

---

## Reference

Pereira, D. H. (2025). Riemannian Manifolds of Asymmetric Equilibria: The Victoria-Nash Geometry.

