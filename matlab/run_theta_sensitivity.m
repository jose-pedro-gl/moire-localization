clear variables;
format short g;

cache_path = "cache/";

sigma = "largestreal";
tol_str = "1e-6";

eigs_impl = @eigs_cupy;
eigs_impl = @(M, N) eigs_matlab(M, N, sigma, eval(tol_str));

pde_impl = @(A, V, bc, Nb) fin_diff_fn_2(A, V, bc, Nb, eigs_impl);

e = exp(1);

Nu = 40 * 5;
Nv = 40 * 5;

TU = linspace(-0.5, 0.5, Nu);
TV = linspace(-0.5, 0.5, Nv);

UV = space_grid_2(TU, TV);

phi_step_str = "pi / 1800";
phi_n_str = "(0:6)";

phi_step = eval(phi_step_str);
phi_n = eval(phi_n_str);

Nb = 1;

L = pi;
LL = 20 * L;

modes = zeros(size(phi_n, 2), Nb);
q_largest_b = zeros(Nv, Nu, size(phi_n, 2));

A = comp2mat(LL);

% 2^2 - 1^2 = 3; 2 * 2 * 1 = 4
% tan = 3 / 4
A1 = square(L);
A2_periodic = comp2mat(e^(j * atan2(3, 4))) * A1;

bc = struct("type", "zero");

for idx = 1:size(phi_n, 2)
    phi = phi_n(idx) * phi_step;

    A2 = comp2mat(e^(j * phi)) * A2_periodic;

    [~, UU1, VV1] = change_basis(UV, A1^(-1) * A);
    [~, UU2, VV2] = change_basis(UV, A2^(-1) * A);

    p1 = 1;
    p2 = 0.4;

    V = V_2016(UU1, VV1, UU2, VV2, p1, p2);
    V = V_2019(UU1, VV1, UU2, VV2, p1, p2, 7);

    V = V_cut(UV(:, :, 1), UV(:, :, 2), 1, 2, -dot([1, 2], [0.48, 0.48]), V);

    start = tic;
    [b, q] = compute_modes(A, V, bc, Nb, 1, pde_impl, cache_path, {});
    toc(start)

    modes(idx, :) = b;
    q_largest_b(:, :, idx) = q;
end

show_figures_theta_sensitivity(UV, A, phi_n, modes, Nb, q_largest_b);
    
%figure(7);
%plot_cells_from = [-3, -3];
%plot_cells_to = [2, 2];
%plot_cells(A1, plot_cells_from, plot_cells_to);
%hold on;
%plot_cells(A2, plot_cells_from, plot_cells_to);
%plot_cells(A, [-1, -1], [0, 0]);
%hold off;
