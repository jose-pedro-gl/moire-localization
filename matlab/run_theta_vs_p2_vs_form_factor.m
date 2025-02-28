clear variables;
format short g;

cache_path = "cache/";

sigma = "largestreal";
tol = 1e-6;

% eigs_impl = @eigs_cupy;
eigs_impl = @(M, N) eigs_matlab(M, N, sigma, tol);

pde_impl = @(A, V, bc, Nb) fin_diff_fn_2(A, V, bc, Nb, eigs_impl);

e = exp(1);

N_cell = 20;

Nu = N_cell * 10;
Nv = N_cell * 10;

TU = linspace(-0.5, 0.5, Nu);
TV = linspace(-0.5, 0.5, Nv);

UV = space_grid_2(TU, TV);

% square
A1 = square(pi);
theta_min = 0;
theta_max = pi / 2;
theta_max_label = "\pi / 2";
theta_target = pi / 4;
theta_target = atan(3 / 4);

% hex
A1 = hex(1);
A1 = A1 / sqrt(abs(det(A1))) * pi;
theta_min = 0;
theta_max = pi / 3;
theta_max_label = "\pi / 3";
theta_target = pi / 6;

theta_values = linspace_quantized(theta_min, theta_max, theta_target, 8, 8);

p2_values = linspace(0, 0.6, 2^5 + 1);

Nb = 1;

LL = N_cell * pi;

q_largest_b = zeros(Nv, Nu, size(theta_values, 2), size(p2_values, 2));

A = square(LL);

[~, UU1, VV1] = change_basis(UV, A1^(-1) * A);

bc = struct("type", "zero");

for idx_theta = 1:size(theta_values, 2)
    theta = theta_values(idx_theta);
    A2 = comp2mat(e^(j * theta)) * A1;

    [~, UU2, VV2] = change_basis(UV, A2^(-1) * A);
    idx_theta
    for idx_p2 = 1:size(p2_values, 2)
        p1 = 1;
        p2 = p2_values(idx_p2);

        V = V_2016(UU1, VV1, UU2, VV2, p1, p2);
    
        start = tic;
        [b, q] = compute_modes(A, V, bc, 2, 2, pde_impl, cache_path, {});
        toc(start)
    
        q_largest_b(:, :, idx_theta, idx_p2) = q(:, :, 1);
    end
end

form_factors = form_factor(q_largest_b, A, abs(det(A)));

figure_latex(1);

[theta_grid, p2_grid] = meshgrid(theta_values, p2_values);
surf(theta_grid ./ theta_max, p2_grid, form_factors.');
xlabel("$\theta / \left(" + theta_max_label + "\right)$"); ylabel("$p_2$"); zlabel("$\chi$");
surf_style(gca());

%figure(7);
%plot_cells_from = [-3, -3];
%plot_cells_to = [2, 2];
%plot_cells(A1, plot_cells_from, plot_cells_to);
%hold on;
%plot_cells(A2, plot_cells_from, plot_cells_to);
%plot_cells(A, [-1, -1], [0, 0]);
%hold off;
