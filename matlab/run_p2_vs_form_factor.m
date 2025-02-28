clear variables;
format short g;

cache_path = "cache/";

sigma = "largestreal";
tol = 1e-6;

% eigs_impl = @eigs_cupy;
eigs_impl = @(M, N) eigs_matlab(M, N, sigma, tol);

pde_impl = @(A, V, bc, Nb) fin_diff_fn_2(A, V, bc, Nb, eigs_impl);

e = exp(1);

Nu = 40 * 5;
Nv = 40 * 5;

TU = linspace(-0.5, 0.5, Nu);
TV = linspace(-0.5, 0.5, Nv);

UV = space_grid_2(TU, TV);

p2_values = linspace(0, 0.6, 2^3 + 1);

L = pi;
LL = 20 * L;

q_largest_b = zeros(Nv, Nu, size(p2_values, 2), 2);

A1 = comp2mat(L);
A = LL / L * A1;

bc = struct("type", "zero");

for periodic = 0:1
    if periodic
        A2 = comp2mat(e^(j * atan2(3, 4))) * A1;
    else
        A2 = comp2mat(e^(j * atan2(sqrt(3), 1))) * A1;
    end

    for idx = 1:size(p2_values, 2)
        [~, UU1, VV1] = change_basis(UV, A1^(-1) * A);
        [~, UU2, VV2] = change_basis(UV, A2^(-1) * A);
    
        p1 = 1;
        p2 = p2_values(idx);

        V = V_2016(UU1, VV1, UU2, VV2, p1, p2);
        %V = V_2019(UU1, VV1, UU2, VV2, p1, p2, 7);
    
        start = tic;
        [b, q] = compute_modes(A, V, bc, 1, 1, pde_impl, cache_path, {});
        toc(start)
    
        q_largest_b(:, :, idx, periodic + 1) = q;
    end
end

show_figures_p2_vs_form_factor(UV, A, p2_values, q_largest_b);

%figure(7);
%plot_cells_from = [-3, -3];
%plot_cells_to = [2, 2];
%plot_cells(A1, plot_cells_from, plot_cells_to);
%hold on;
%plot_cells(A2, plot_cells_from, plot_cells_to);
%plot_cells(A, [-1, -1], [0, 0]);
%hold off;
