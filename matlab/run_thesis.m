clear variables;
format short g;

e = exp(1);

cache_path = "cache/";
export_path = {};

sigma = "largestreal";
tol = 1e-6;

eigs_impl = @(M, N) eigs_matlab(M, N, sigma, tol);

analysis = struct("type", "");

tan_theta_str = "$\tan\left(\theta\right) = ";

% 2016
% p1 = 1, p2 = 0.4
%

L_cell = pi;
N_cell = 20;
L = N_cell * L_cell;

p1 = 1;
p2 = 0.4;

A = comp2mat(L);

% % % square
% % phi = atan(3 / 4); %+ 6 * pi / 1800;
% % phi = atan(sqrt(3));
% % 
% % A1 = comp2mat(L_cell);
% % A2 = comp2mat(e^(j * phi)) * A1;
% % p2 = 0.4;
% % 
% % % p2_values = linspace(0, 0.6, 2^3 + 1);
% % % theta_values = [atan(3 / 4), atan(sqrt(3))];
% % % theta_values = [atan(3 / 4) + (0:6) * pi / 1800, atan(sqrt(3))];
% % 
% % p2_values = linspace(0, 0.6, 2^5 + 1);
% % p2_values(2 * [(1:2), (8:16)]) = [];
% % 
% % [A1, A2, A] = square_moire(L_cell, 1, 2);
% % 
% % analysis = struct("type", "periodic");
% % 
% % % analysis = struct(...
% % %     "type", "p2_vs_form_factor",...
% % %     "p2_values", p2_values, ...
% % %     "sweep_values", [atan(3 / 4), atan(sqrt(3))], ...
% % %     "sweep_labels", [tan_theta_str + "3 / 4$", tan_theta_str + "\sqrt{3}$"]);


% % % hex
% % phi = atan(sqrt(3) * (-8) / (2 * 3 - 8));
% % % phi = acos(11 / 14);
% % % phi = atan(1);
% % % phi = atan(sqrt(3) * (-21) / (2 * 16 - 21));
% % 
% % A1 = comp2vec([1, e^(j * pi / 3)]); A1 = L_cell / sqrt(abs(det(A1))) * A1;
% % A2 = comp2mat(e^(j * phi)) * A1;
% % 
% % p2_values = linspace(0, 0.6, 2^5 + 1);
% % p2_values(2 * [(1:2), (8:16)]) = [];
% % 
% % [A1, A2, A] = hex_moire(L_cell, 1, 2);
% % p2 = 0.4;
% % 
% % analysis = struct("type", "periodic");
% % 
% % % analysis = struct(...
% % %     "type", "p2_vs_form_factor",...
% % %     "p2_values", p2_values, ...
% % %     "sweep_values", [atan(sqrt(3) * (-8) / (2 * 3 - 8)), atan(1)], ...
% % %     "sweep_labels", [tan_theta_str + "4 \sqrt{3}$", tan_theta_str + "1$"]);


% % % rect
% % % % alpha = pi;
% % % % n1 = 2; n2 = -1;
% % % % % n1 = 2; n2 = 3;
% % % % phi = atan(-2 * alpha * n1 * n2 / (n1^2 - (alpha * n2)^2));
% % % % p2 = 0.8;
% % 
% % alpha_sq_num = 2;
% % alpha_sq_den = 1;
% % alpha = sqrt(alpha_sq_num / alpha_sq_den);
% % 
% % % equation: a^2 + num * den * b^2 = c^2 => a^2 + 2b^2 = c^2
% % % solution for a, b, c != 0 => 1^2 + 2 * 2^2 = 3^2
% % a = 1;
% % b = 2;
% % phi = atan(sqrt(alpha_sq_num * alpha_sq_den) * b / a);
% % 
% % %alpha = pi;
% % %phi = atan(1);
% % 
% % A1 = comp2vec([L_cell / sqrt(alpha), j * sqrt(alpha) * L_cell]);
% % A2 = comp2mat(e^(j * phi)) * A1;
% % 
% % V_fn = @(UU, VV) V_2016(UU, VV, UU, VV, 1, 0);
% % 
% % % % n1 = -1;
% % % % n2 = 1;
% % % % A = A1 * [n1, -alpha_sq_num * n2; n2, alpha_sq_den * n1];
% % % % p2 = 0.4;
% % % % 
% % % % analysis = struct("type", "periodic");
% % 
% % p2_values = linspace(0, 0.6, 2^5 + 1);
% % p2_values(2 * [(1:2), (8:16)]) = [];
% % p2_values_2 = linspace(0.6, 0.6 + 1 * 0.6, 1 * 2^3 + 1);
% % p2_values = [p2_values, p2_values_2(2:end)];
% % 
% % alpha_values = [sqrt(2), pi, sqrt(2), pi];
% % alpha_labels = ["\sqrt{2}", "\pi", "\sqrt{2}", "\pi"];
% % 
% % tan_theta_values = [2 * sqrt(2), 4 * pi / (4 - pi^2), 1, 1];
% % tan_theta_labels = ["2 \sqrt{2}", "4 \pi / \left(4 - \pi^2\right)", "1", "1"];
% % 
% % A1_fn = @(alpha)...
% %     comp2vec([L_cell / sqrt(alpha), j * sqrt(alpha) * L_cell]);
% % 
% % V_fn = @(UV, A1, theta, p2)...
% %     V_moire(UV, A1, comp2mat(e^(j * theta)) * A1, A, V_fn, p1, p2);
% % % 
% % % analysis = struct(...
% % %     "type", "p2_vs_form_factor",...
% % %     "p2_values", p2_values, ...
% % %     "sweep_values", struct(...
% % %         "alpha", num2cell(alpha_values), ...
% % %         "theta", num2cell(atan(tan_theta_values))), ...
% % %     "sweep_labels", tan_theta_str + tan_theta_labels + ", \alpha = " + alpha_labels + "$", ...
% % %     "V_fn", @(UV, sweep, p2) V_fn(UV, A1_fn(sweep.alpha), sweep.theta, p2));


% % % centered rect
% % 
% % cos_alpha_num = 1;
% % cos_alpha_den = 3;
% % alpha = acos(cos_alpha_num / cos_alpha_den);
% % % equation: a^2 + (den^2 - num^2) b^2 = c^2 => a^2 + 8b^2 = c^2
% % % solution for a, b, c != 0 => 1^2 + 8 * 1^2 = 3^2
% % a = 1;
% % b = 1;
% % phi = atan(sqrt(cos_alpha_den^2 - cos_alpha_num^2) * b / a);
% % %phi = atan(1);
% % 
% % % % n1 = 1;
% % % % n2 = -2;
% % % % alpha = acos(sqrt(3) / 2);
% % % % phi = atan((n1^2 - n2^2) * sin(alpha) / ((n1^2 + n2^2) * cos(alpha) + 2 * n1 * n2));
% % % % phi = atan(1);
% % % % p2 = 0.8;
% % % % p2 = 1.6;
% % 
% % A1 = comp2vec([1, e^(j * alpha)]); A1 = L_cell / sqrt(abs(det(A1))) * A1;
% % A2 = comp2mat(e^(j * phi)) * A1;
% % 
% % V_fn = @(UU, VV) V_2016(UU, VV, UU, VV, 1, 0);
% % 
% % n1 = -1;
% % n2 = 1;
% % A = A1 * [3, 0; 0, 1];
% % p2 = 0.1;
% % 
% % analysis = struct("type", "periodic");
% % 
% % p2_values = linspace(0, 0.6, 2^5 + 1);
% % p2_values(2 * [(1:2), (8:16)]) = [];
% % p2_values_2 = linspace(0.6, 0.6 + 2 * 0.6, 2 * 2^3 + 1);
% % p2_values = [p2_values, p2_values_2(2:end)];
% % 
% % cos_alpha_values = [1 / 3, sqrt(3) / 2, 1 / 3, sqrt(3) / 2];
% % cos_alpha_labels = ["1 / 3", "\sqrt{3} / 2", "1 / 3", "\sqrt{3} / 2"];
% % 
% % tan_theta_values = [2 * sqrt(2), -3 / 11 * (5 * sqrt(3) + 8), 1, 1];
% % tan_theta_labels = ["2 \sqrt{2}", "-3 / 11\left(5 \sqrt{3} + 8\right)", "1", "1"];
% % 
% % A1_fn = @(A1) L_cell / sqrt(abs(det(A1))) * A1;
% % A1_fn = @(alpha) A1_fn(comp2vec([1, e^(j * alpha)]));
% % 
% % V_fn = @(UV, A1, theta, p2)...
% %     V_moire(UV, A1, comp2mat(e^(j * theta)) * A1, A, V_fn, p1, p2);
% % 
% % % analysis = struct(...
% % %     "type", "p2_vs_form_factor",...
% % %     "p2_values", p2_values, ...
% % %     "sweep_values", struct(...
% % %         "alpha", num2cell(acos(cos_alpha_values)), ...
% % %         "theta", num2cell(atan(tan_theta_values))), ...
% % %     "sweep_labels", tan_theta_str + tan_theta_labels + ...
% % %         ", \cos\left(\varphi\right) = " + cos_alpha_labels + "$", ...
% % %     "V_fn", @(UV, sweep, p2) V_fn(UV, A1_fn(sweep.alpha), sweep.theta, p2));


V_fn = @(UU, VV) V_2016(UU, VV, UU, VV, 1, 0);

if analysis.type ~= "p2_vs_form_factor"
    V_fn = @(UV) V_moire(UV, A1, A2, A, V_fn, p1, p2);
elseif isfield(analysis, "V_fn")
    V_fn = analysis.V_fn;
else
    V_fn = @(UV, theta, p2) V_moire(UV, A1, comp2mat(e^(j * theta)) * A1, A, V_fn, p1, p2);
end

Nb = 150;
Nu = N_cell * 10;
Nv = N_cell * 10;

Bn0 = [1];

if analysis.type == "periodic"
    Nk = 2^3;
    
    % K in normalized coordinates in reciprocal space
    K = build_path(1 / 2 * [0, 0; 1, 0; 1, 1; 0, 0], Nk);
    K_ticks = 1 + Nk * (0:3);
    K_tick_labels = ["$\Gamma$", "$M$", "$X$", "$\Gamma$"];

    N_cell = ceil(sqrt(abs(det(A) / det(A1))));
    Nu = N_cell * 10;
    Nv = N_cell * 10;
    Nb = 5;
end



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

ws = 3.3;
wx = 0.25;
wy = 0.75;

p1 = 1.95;
p2 = 1.95;

phi = atan(3 / 4);

% case 1
% 
% N_lattice = 11;
% L = 50;
% Nb = 71;
% Bn0 = [1; 21; 18];

% case 2
% 
% N_lattice = 21;
% L = 84;
% Nb = 248;
% Bn0 = [1; 79; 71];

% N_lattice_u = N_lattice;
% N_lattice_v = N_lattice;
% Lu = L
% Lv = L

N_lattice_u = 31;
N_lattice_v = 11;
Lu = 116;
Lv = 50;
Nb = 150;
Bn0 = [];

L_lattice_u = N_lattice_u * ws;
L_lattice_v = N_lattice_v * ws;

R = comp2mat(e^(j * (pi / 4 - pi / 30)));

R = comp2mat(e^(j * pi / 4));
phi = phi + pi / 30;

A = comp2vec([Lu, j * Lv]);
A1 = R * square(ws);
A2 = comp2mat(e^(j * phi)) * A1;

Aw = diag([1 / wx, 1 / wy]);

N_lattice_max = max(N_lattice_u, N_lattice_v);
lattice0 = space_grid_2(-N_lattice_max:N_lattice_max, -N_lattice_max:N_lattice_max);

lattice1 = change_basis(lattice0, A1);
mask1 = V_trunc(ones(size(lattice1)), lattice1(:, :, 1), lattice1(:, :, 2), L_lattice_u / 2, L_lattice_v / 2);

lattice2 = change_basis(lattice0, A2);
mask2 = V_trunc(ones(size(lattice2)), lattice2(:, :, 1), lattice2(:, :, 2), L_lattice_u / 2, L_lattice_v / 2);

V_fn = @(XX, YY, dx, dy) V_gaussian(cat(3, XX - dx, YY - dy), Aw);
V_fn = @(XX, YY, lattice, mask)...
    V_lattice(lattice, mask, @(dx, dy) V_fn(XX, YY, dx, dy), size(XX, 2), size(XX, 1));
V_fn = @(XX, YY) p1 * V_fn(XX, YY, lattice1, mask1) + p2 * V_fn(XX, YY, lattice2, mask2);
V_fn = @(XY, UV) V_fn(XY(:, :, 1), XY(:, :, 2));
V_fn = @(UV) V_fn(change_basis(UV, A), UV);

Nu = round(Lu / ws * 15) + 1;
Nv = round(Lv / ws * 15) + 1;



TU = linspace(-0.5, 0.5, Nu);
TV = linspace(-0.5, 0.5, Nv);

[UV, UU, VV] = space_grid_2(TU, TV);

% export_path = "export/thesis_finite_square/";

if analysis.type == ""
    % V = V_fn(UV);
    % [b, q] = run_periodic_fn(A, V, K, Nb, 1, eigs_impl, cache_path, export_path);
    % show_figures_periodic(UV, A, V, K, K_ticks, K_tick_labels, b, Nb, q, (20 * pi)^2);
    
    V = V_fn(UV);
    [b, q] = run_aperiodic_fn(A, V, Nb, eigs_impl, cache_path, export_path);
    show_figures_aperiodic(UV, A, V, b, Nb, q, Bn0, true);
end

if analysis.type == "periodic"
    V = V_fn(UV);
    [b, q] = run_periodic_fn(A, V, K, Nb, 2, eigs_impl, cache_path, export_path);
    show_figures_periodic(UV, A, V, K, K_ticks, K_tick_labels, b, Nb, q, L^2);
end

if analysis.type == "p2_vs_form_factor"
    p2_values = analysis.p2_values;
    sweep_values = analysis.sweep_values;
    sweep_labels = analysis.sweep_labels;

    V_fn = @(sweep, p2) V_fn(UV, sweep, p2);

    q = run_p2_vs_form_factor_fn(A, V_fn, p2_values, sweep_values, eigs_impl, cache_path, export_path);
    show_figures_p2_vs_form_factor(UV, A, p2_values, sweep_labels, q);
end
