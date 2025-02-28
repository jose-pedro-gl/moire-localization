clear variables;
format short g;

cache_path = "cache/";

method = 1; % 1: finite differences; 2: plane wave

sigma = "largestreal";
tol_str = "1e-6";

eigs_impl = @(M, N) eigs_matlab(M, N, sigma, eval(tol_str));

Nu = 60;
Nv = 60;
Nb = 120;

Nb2 = 5;

Nb = min(Nb, Nv * Nu);
Nb2 = min(Nb2, Nb);

TU = linspace_ex(-0.5, 0.5, Nu);
TV = linspace_ex(-0.5, 0.5, Nv);

UV = space_grid_2(TU, TV);

Nk = 2^2;

L_str = "pi"; L = eval(L_str);

% 2^2 - 1^2 = 3; 2 * 2 * 1 = 4
% tan = 3 / 4
n1_str = "1"; n1 = eval(n1_str);
n2_str = "2"; n2 = eval(n2_str);

% sqrt(3) approximation
%
% 4^2 - 1^2 = 15; 2 * 4 * 1 = 8
% tan = 15 / 8
%
% 15^2 - 4^2 = 209; 2 * 15 * 4 = 120
% tan = 209 / 120
%
% 56^2 - 15^2 = 2911; 2 * 56 * 15 = 1680
% tan = 2911 / 1680
%
%[A1, A2, A] = square(L, 1, 4);
%[A1, A2, A] = square(L, 15, 4);
%[A1, A2, A] = square(L, 15, 56);

lattice_type = "square_moire";

if lattice_type == "square_moire"
    [A1, A2, A] = square_moire(L, n1, n2);
elseif lattice_type == "hex_moire"
    [A1, A2, A] = hex_moire(L / sqrt(sqrt(3) / 2), n1, n2);
end

[~, UU1, VV1] = change_basis(UV, A1^(-1) * A);
[~, UU2, VV2] = change_basis(UV, A2^(-1) * A);

p1 = 1;
p2 = 0;
E0 = 7;

V = V_2016(UU1, VV1, UU2, VV2, p1, p2);
%V = V_2019(UU1, VV1, UU2, VV2, p1, p2, E0);

% K in normalized coordinates in reciprocal space
K = build_path(1 / 2 * [0, 0; 1, 0; 1, 1; 0, 0], Nk);
K_ticks = 1 + Nk * (0:3);
K_tick_labels = ["$\Gamma$", "$M$", "$X$", "$\Gamma$"];

if method == 1
    pde_impl = @(A, V, bc, Nb) fin_diff_fn_2(A, V, bc, Nb, eigs_impl);
else
    pde_impl = @(A, V, bc, Nb) plane_wave_fn_2(A, V, bc, Nb, eigs_impl);
end

start = tic;
[bands, q_largest_b] = compute_bands(A, V, K, Nb, 1, pde_impl, cache_path, {});
toc(start)

LL = 40 * pi;
show_figures(UV, A, V, K, K_ticks, K_tick_labels, bands, Nb2, q_largest_b, LL * LL);

%figure(7);
%plot_cells_from = [-3, -3];
%plot_cells_to = [2, 2];
%plot_cells(A1, plot_cells_from, plot_cells_to);
%hold on;
%plot_cells(A2, plot_cells_from, plot_cells_to);
%plot_cells(A, [-1, -1], [0, 0]);
%hold off;
