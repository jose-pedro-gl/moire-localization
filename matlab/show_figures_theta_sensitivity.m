function show_figures_theta_sensitivity(UV, A, phi_n, modes, Nb, q_largest_b)
    area = abs(det(A));

    [~, XX, YY] = change_basis(UV, A);
    % 
    % figure_latex(1);
    % surf(XX, YY, V);
    % xlabel("$x$"); ylabel("$y$"); zlabel("$V$");
    % surf_style(gca());
    
    figure_latex(2);
    hold off;
    
    for idx = 1:size(phi_n, 2)
        scatter(0:(Nb - 1), modes(idx, :), [], ...
            "DisplayName", "n = " + phi_n(idx), ...
            "MarkerFaceColor", "white", ...
            "LineWidth", 0.25);
    
        hold on;
    end
    
    hold off;
    
    xlabel("$n$");
    ylabel("$\beta_n$");
    legend();
    
    form_factors = form_factor(q_largest_b, A, area);
    
    figure_latex(3);
    plot(phi_n, form_factors, "LineWidth", 2);
    xlabel("$n$"); ylabel("$\chi_n$");
    
    q_s = q_largest_b(:, :, 1);

    figure_latex(4);
    surf(XX, YY, q_s .* conj(q_s));
    xlabel("$x$"); ylabel("$y$"); zlabel("$\left|\psi\right|^2$");
    surf_style(gca());

    figure_latex(5);
    surf(XX, YY, log(q_s .* conj(q_s)));
    xlabel("$x$"); ylabel("$y$"); zlabel("$\ln\left|\psi\right|^2$");
    
    q_e = q_largest_b(:, :, end);

    figure_latex(6);
    surf(XX, YY, q_e .* conj(q_e));
    xlabel("$x$"); ylabel("$y$"); zlabel("$\left|\psi\right|^2$");
    surf_style(gca());
end
