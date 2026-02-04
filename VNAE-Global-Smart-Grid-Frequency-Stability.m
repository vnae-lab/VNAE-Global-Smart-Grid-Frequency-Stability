% -------------------------------------------------------------------------
% VNAE Framework: Global Smart Grid Frequency Stability (n=100,000)
% Purpose: Validating Volume Contraction in Massive Energy Networks
% -------------------------------------------------------------------------
clear; clc;

% Set seed for reproducibility
rng(42);

% --- [ GLOBAL GRID CONFIGURATION ] ---
n = 100000;             % 100,000 Energy Nodes (National Grid Scale)
beta_val = 0.12;        % Grid Rigidity (VNAE Structural Factor)

fprintf('Initializing National Grid Simulation with %d nodes...\n', n);

% 1. GENERATE HETEROGENEOUS NODE INERTIA (THETA)
% Asymmetry in theta creates the stabilizing manifold
theta_nodes = 1.0 + (120.0 - 1.0) .* rand(n, 1);

% 2. GENERATE ULTRA-SPARSE TRANSMISSION NETWORK (A)
% Average of 15 transmission lines per node
nz_count = n * 15;
i = randsample(n, nz_count, true);
j = randsample(n, nz_count, true);
weights = exprnd(1/0.8, [nz_count, 1]);

% Construct Sparse Adjacency Matrix
A_sparse = sparse(i, j, weights, n, n);
A_sparse = (A_sparse + A_sparse') / 2; % Symmetrize grid
A_sparse = spdiags(zeros(n,1), 0, A_sparse); % Remove self-loops

% 3. STATISTICAL MANIFOLD CURVATURE (K) CERTIFICATION
% Monte Carlo estimation of Scalar Curvature for massive N
sample_size = 500000;
idx_i = randi(n, sample_size, 1);
idx_j = randi(n, sample_size, 1);

diff_theta = abs(theta_nodes(idx_i) - theta_nodes(idx_j));
% Linear indexing for fast sparse access
couplings = diag(A_sparse(idx_i, idx_j)) + diag(A_sparse(idx_j, idx_i));
rigidity = 1 + beta_val .* (theta_nodes(idx_i) + theta_nodes(idx_j));

K_estimate = mean((diff_theta .* couplings) ./ rigidity);

% -------------------------------------------------------------------------
% DYNAMIC SIMULATION (FLOW ON THE ENERGY MANIFOLD)
% -------------------------------------------------------------------------

% Construct Laplacian: L = D - A
d_vec = sum(A_sparse, 2);
L_sparse = spdiags(d_vec, 0, n, n) - A_sparse;

p_weather_shocks = normrnd(0, 10, [n, 1]); % Weather noise
omega_init = normrnd(0, 5, [n, 1]);       % Initial state

% Simulation parameters
dt = 0.01;
n_steps = 10;
history = zeros(n_steps + 1, n);
history(1, :) = omega_init';

fprintf('Running massive scale simulation (National Grid Scale)...\n');
tic;
curr_omega = omega_init;
for t = 1:n_steps
    % VNAE fundamental equation: d_omega = -(L + Theta) * omega + p
    % MATLAB's sparse engine is highly optimized for this multiplication
    d_omega = -(L_sparse * curr_omega) - (theta_nodes .* curr_omega) + p_weather_shocks;
    curr_omega = curr_omega + dt * d_omega;
    history(t+1, :) = curr_omega';
end
sim_time = toc;

% -------------------------------------------------------------------------
% FINAL STABILITY REPORT
% -------------------------------------------------------------------------
fprintf('\n--- VNAE GLOBAL ENERGY GRID REPORT ---\n');
fprintf('Total Grid Nodes (n): %d\n', n);
fprintf('Structural Curvature (K): %.8f\n', K_estimate);
fprintf('Simulation Time: %.4f seconds\n', sim_time);
fprintf('Status: %s\n', char(datetime('now')));

% Visualization: Sample 100 nodes
plot_sample = randsample(n, 100);
time_axis = (0:n_steps) * dt;

figure('Color', 'w', 'Position', [100, 100, 900, 500]);
plot(time_axis, history(:, plot_sample), 'LineWidth', 1);
title(sprintf('VNAE 100k Nodes Grid Stability (K = %.4f)', K_estimate));
xlabel('Time (seconds)');
ylabel('Frequency Deviation (Hz)');
grid on;
