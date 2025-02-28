clear variables;
format short g;

cache_path = "cache/";

sigma = "largestreal";
tol_str = "1e-6";

% eigs_impl = @eigs_cupy;
eigs_impl = @(M, N) eigs_matlab(M, N, sigma, eval(tol_str));

pde_impl = @(A, V, bc, Nb) fin_diff_fn_2(A, V, bc, Nb, eigs_impl);

e = exp(1);

Nu = 11 * 20 + 1;
Nv = 11 * 20 + 1;

TU = linspace(-0.5, 0.5, Nu);
TV = linspace(-0.5, 0.5, Nv);

[UV, UU, VV] = space_grid_2(TU, TV);

Nb = 72;

L = pi;
LL = 14 * L;

phi = atan(sqrt(3));

A = comp2mat(LL);

% 2^2 - 1^2 = 3; 2 * 2 * 1 = 4
% tan = 3 / 4
A1 = square(L);
A2 = comp2mat(e^(j * phi)) * A1;

bc = struct("type", "zero");

[UV1, UU1, VV1] = change_basis(UV, A1^(-1) * A);
[UV2, UU2, VV2] = change_basis(UV, A2^(-1) * A);

p1 = 1;
p2 = 0.4;

[XY, XX, YY] = change_basis(UV, A);

V = V_2016(UU1, VV1, UU2, VV2, p1, p2);
%V = V_2019(UU1, VV1, UU2, VV2, p1, p2, 7);

%Aw = [1 / 0.25, 0; 0, 1 / 0.75];

%XY = change_basis(UV, A);

%V = V_gaussian(XY, A1, A2, p1, p2, Aw, 11, 11);
%V = V .* (abs(UU) < 11/14/2) .* (abs(VV) < 11/14/2);

%V = V_cut(UU, VV, 1, 2, dot([1, 2], [-0.2, -0.2]), V);

start = tic;
[b, q] = compute_modes(A, V, bc, Nb, Nb, pde_impl, cache_path, {});
toc(start)

show_figures_aperiodic(UV, A, V, b, Nb, q);
