function show_figures_aperiodic(UV, A, V, b, Nb, q, Bn0, highlight_modes)
    area = abs(det(A));
    
    form_factors = form_factor(q, A, area);

    [~, Bn] = sort(form_factors, "descend");
    Bn = setdiff(Bn, Bn0, "stable");
    Bn = [Bn0; Bn(1:(4 - length(Bn0)))];

    [~, XX, YY] = change_basis(UV, A);

    matlab_blue = [0, 0.4470, 0.7410];
    matlab_orange = [0.8500, 0.3250, 0.0980];

    cm_scatter = repmat(matlab_blue, Nb, 1);

    if highlight_modes
        cm_scatter(Bn, :) = repmat(matlab_orange, size(Bn));
    end

    figure_latex(1);
    surf(XX, YY, V);
    xlabel("$x$"); ylabel("$y$"); zlabel("$V$");
    title("$V$");
    surf_style(gca());

    figure_latex(2);

    scatter(0:(Nb - 1), b, [], cm_scatter, ...
        "MarkerFaceColor", "white");

    xlabel("$n$");
    ylabel("$\beta_n$");

    figure_latex(3);

    scatter(0:(Nb - 1), form_factors, [], cm_scatter, ...
        "MarkerFaceColor", "white");
    
    xlabel("$n$"); ylabel("$\chi_n$");

    figure_latex(4);

    for idx = 1:length(Bn)
        bn = Bn(idx);
        q_name = "$\left|\psi_{n = " + (bn - 1) + "}\right|^2$";

        subplot(2, 2, idx);
        surf(XX, YY, q(:, :, bn) .* conj(q(:, :, bn)));
        xlabel("$x$"); ylabel("$y$"); zlabel(q_name);
        title(q_name);
        surf_style(gca());

        % cm = hot();
        % cm = cm(17:end, :);
        % colormap(cm);
    end

    figure_latex(5);
    surf(XX, YY, q(:, :, 1) .* conj(q(:, :, 1)));
    xlabel("$x$"); ylabel("$y$"); zlabel("$\left|\psi\right|^2$");
    title("$\left|\psi\right|^2$");
    surf_style(gca());
end
