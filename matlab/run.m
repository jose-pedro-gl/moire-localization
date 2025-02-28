clear variables;
format short g;

e = exp(1);

cache_path = "cache/";
export_path = {};

sigma = "largestreal";
tol = 1e-6;

% eigs_impl = @eigs_cupy;
eigs_impl = @(M, N) eigs_matlab(M, N, sigma, tol);


% 2016
% p1 = 1, p2 = 0.4
%

% % L_cell = pi;
% % N_cell = 20;
% % L = N_cell * L_cell;
% % 
% % p1 = 1;
% % p2 = 0.4;
% % 
% % phi = atan(3 / 4); %+ 6 * pi / 1800;
% % 
% % A = comp2mat(L);
% % A1 = comp2mat(L_cell);
% % A2 = comp2mat(e^(j * phi)) * A1;
% % 
% % V_fn = @(UU, VV) V_2016(UU, VV, UU, VV, 1, 0);
% % V_fn = @(UV) V_moire(UV, A1, A2, A, V_fn, p1, p2);
% % 
% % Nb = 150;
% % Nu = N_cell * 10;
% % Nv = N_cell * 10;
% % 
% % Bn0 = [1];


% % Nk = 2^3;
% % 
% % % K in normalized coordinates in reciprocal space
% % K = build_path(1 / 2 * [0, 0; 1, 0; 1, 1; 0, 0], Nk);
% % K_ticks = 1 + Nk * (0:3);
% % K_tick_labels = ["$\Gamma$", "$M$", "$X$", "$\Gamma$"];
% % 
% % Nb = 5;
% % Nu = 30;
% % Nv = 30;
% % 
% % [A1, A2, A] = square_moire(L_cell, 1, 2);
% % 
% % V_fn = @(UU, VV) V_2016(UU, VV, UU, VV, 1, 0);
% % V_fn = @(UV) V_moire(UV, A1, A2, A, V_fn, p1, p2);


% 2019
% p1 = 1, p2 = 0.2
%

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

L_cell = pi;

Nk = 2^2;

% K in normalized coordinates in reciprocal space
K = build_path(1 / 2 * [0, 0; 1, 0; 1, 1; 0, 0], Nk);
K_ticks = 1 + Nk * (0:3);
K_tick_labels = ["$\Gamma$", "$M$", "$X$", "$\Gamma$"];

p1 = 1;
p2 = 0.1;
E0 = 7;

Nb = 15;
Nu = ceil(sqrt(15^2 + 56^2) * 10) + 1;
Nv = Nu;

%[A1, A2, A] = square_moire(L_cell, 1, 4);
%[A1, A2, A] = square_moire(L_cell, 15, 4);
[A1, A2, A] = square_moire(L_cell, 15, 56);

V_fn = @(UU, VV) V_2016(UU, VV, UU, VV, 1, 0);
V_fn = @(UV) V_moire(UV, A1, A2, A, V_fn, p1, p2);
V_fn = @(UV) -E0 ./ (1 + V_fn(UV) .* conj(V_fn(UV)));


% 2023
% ws = 3.3, wx = 0.25, wy = 0.75
% p1 = p2 = 1.95
% square(ws, 1, 2) -> (3, 4, 5) -> phi = atan(3 / 4)
%
% case 1
% L_lattice = 11 * ws
% L = 50
%
% case 2
% L_lattice = 21 * ws
% L = 84
%

% % ws = 3.3;
% % wx = 0.25;
% % wy = 0.75;
% % 
% % p1 = 1.95;
% % p2 = 1.95;
% % 
% % phi = atan(3 / 4);
% % 
% % % case 1
% % % 
% % N_lattice = 11;
% % L = 50;
% % Nb = 71;
% % Bn0 = [1; 21; 18];
% % 
% % % case 2
% % % 
% % % N_lattice = 21;
% % % L = 84;
% % % Nb = 248;
% % % Bn0 = [1; 79; 71];
% % 
% % L_lattice = N_lattice * ws;
% % 
% % A = comp2mat(L);
% % A1 = square(ws);
% % A2 = comp2mat(e^(j * phi)) * A1;
% % 
% % Aw = diag([1 / wx, 1 / wy]);
% % 
% % lattice0 = space_grid_2(-N_lattice:N_lattice, -N_lattice:N_lattice);
% % 
% % lattice1 = change_basis(lattice0, A1);
% % mask1 = V_trunc(ones(size(lattice1)), lattice1(:, :, 1), lattice1(:, :, 2), L_lattice / 2, L_lattice / 2);
% % 
% % lattice2 = change_basis(lattice0, A2);
% % mask2 = V_trunc(ones(size(lattice2)), lattice2(:, :, 1), lattice2(:, :, 2), L_lattice / 2, L_lattice / 2);
% % 
% % V_fn = @(XX, YY, dx, dy) V_gaussian(cat(3, XX - dx, YY - dy), Aw);
% % V_fn = @(XX, YY, lattice, mask)...
% %     V_lattice(lattice, mask, @(dx, dy) V_fn(XX, YY, dx, dy), size(XX, 2), size(XX, 1));
% % V_fn = @(XX, YY) p1 * V_fn(XX, YY, lattice1, mask1) + p2 * V_fn(XX, YY, lattice2, mask2);
% % V_fn = @(XY, UV) V_fn(XY(:, :, 1), XY(:, :, 2));
% % V_fn = @(UV) V_fn(change_basis(UV, A), UV);
% % 
% % Nu = round(L / ws * 15) + 1;
% % Nv = round(L / ws * 15) + 1;

% workbench

% % ws = 3.3;
% % wx = 0.25;
% % wy = 0.75;
% % 
% % f = 1;
% % 
% % p1 = f * 1.95;
% % p2 = f * 1.95;
% % 
% % phi = atan(3 / 4);
% % 
% % N_lattice_u = 51;
% % N_lattice_v = 21;
% % Lu = 180;
% % Lv = 84;
% % Nb = 200;
% % Bn0 = [];
% % 
% % N_lattice_max = max(N_lattice_u, N_lattice_v);
% % 
% % L_lattice_u = N_lattice_u * ws;
% % L_lattice_v = N_lattice_v * ws;
% % 
% % R = comp2mat(e^(j * (pi/4-pi/30)));
% % 
% % A = diag([Lu, Lv]);
% % A1 = R * square(ws);
% % A2 = comp2mat(e^(j * phi)) * A1;
% % 
% % Aw = diag([1 / wx, 1 / wy]);
% % 
% % % u = [1, 1];
% % % p = [-1, -1]*ws*1e-6;
% % % c = dot(u, p);
% % % u_ = [-1, 1];
% % % p_ = [1, -1]*ws*1e-6;
% % % c_ = dot(u_, p_);
% % 
% % u = [1, 1];
% % p = [-11, -4] * ws * (1 + 1e-6);
% % c = dot(u, p);
% % 
% % lattice0 = space_grid_2(-N_lattice_max:N_lattice_max, -N_lattice_max:N_lattice_max);
% % 
% % lattice1 = change_basis(lattice0, A1);
% % mask1 = V_trunc(ones(size(lattice1)), lattice1(:, :, 1), lattice1(:, :, 2), L_lattice_u / 2, L_lattice_v / 2);
% % %mask1 = V_cut(lattice1(:, :, 1), lattice1(:, :, 2), u(1), u(2), c, mask1);
% % %mask1 = V_cut(lattice1(:, :, 1), lattice1(:, :, 2), u_(1), u_(2), c_, mask1);
% % 
% % lattice2 = change_basis(lattice0, A2);
% % mask2 = V_trunc(ones(size(lattice2)), lattice2(:, :, 1), lattice2(:, :, 2), L_lattice_u / 2, L_lattice_v / 2);
% % %mask2 = V_cut(lattice2(:, :, 1), lattice2(:, :, 2), u(1), u(2), c, mask2);
% % %mask2 = V_cut(lattice2(:, :, 1), lattice2(:, :, 2), u_(1), u_(2), c_, mask2);
% % 
% % V_fn = @(XX, YY, dx, dy) V_gaussian(cat(3, XX - dx, YY - dy), Aw);
% % V_fn = @(XX, YY, lattice, mask)...
% %     V_lattice(lattice, mask, @(dx, dy) V_fn(XX, YY, dx, dy), size(XX, 2), size(XX, 1));
% % V_fn = @(XX, YY) p1 * V_fn(XX, YY, lattice1, mask1) + p2 * V_fn(XX, YY, lattice2, mask2);
% % V_fn = @(XY, UV) V_fn(XY(:, :, 1), XY(:, :, 2));
% % V_fn = @(UV) V_fn(change_basis(UV, A), UV);
% % 
% % % V_fn = @(UU, VV) V_2016(UU, VV, UU, VV, 1, 0);
% % % V_fn = @(UV, UU, VV) V_cut(UU, VV, u(1), u(2), c, V_moire(UV, A1, A2, A, V_fn, p1, p2));
% % % V_fn = @(UV) V_fn(UV, UV(:, :, 1), UV(:, :, 2));
% % 
% % Nu = round(Lu / ws * 15) + 1;
% % Nv = round(Lv / ws * 15) + 1;




% % L_cell = pi;
% % N_cell = 20;
% % L = N_cell * L_cell;
% % 
% % p1 = 1;
% % p2 = 0.4;
% % 
% % % alpha = pi;
% % % n1 = 1;
% % % n2 = 2;
% % % phi = atan(-2 * alpha * n1 * n2 / (n1^2 - (alpha * n2)^2)) + pi / 300;
% % 
% % phi = atan(3 / 4);
% % 
% % A = comp2vec([L, j * L]);
% % %A1 = comp2vec([L_cell / sqrt(alpha), j * sqrt(alpha) * L_cell]);
% % A1 = comp2vec([L_cell, j * L_cell]);
% % A2 = comp2mat(e^(j * phi)) * A1;
% % 
% % V_fn = @(UU, VV) V_2016(UU, VV, UU, VV, 1, 0);
% % V_fn = @(UV) V_moire(UV, A1, A2, A, V_fn, p1, p2);
% % %V_fn = @(UV) V_fn(UV) + 0.05*UV(:, :, 1) + 0.05*UV(:,:,2);
% % 
% % Nb = 150;
% % Nu = N_cell * 10;
% % Nv = N_cell * 10;
% % 
% % Bn0 = [1];


%export_path = "export/V_2023__triangle_right__11_ws/";

TU = linspace(-0.5, 0.5, Nu);
TV = linspace(-0.5, 0.5, Nv);

[UV, UU, VV] = space_grid_2(TU, TV);

V = V_fn(UV);

[b, q] = run_periodic_fn(A, V, K, Nb, 1, eigs_impl, cache_path, export_path);
show_figures_periodic(UV, A, V, K, K_ticks, K_tick_labels, b, Nb, q, (20 * pi)^2);

% [b, q] = run_aperiodic_fn(A, V, Nb, eigs_impl, cache_path, export_path);
% show_figures_aperiodic(UV, A, V, b, Nb, q, Bn0, false);
