function show_figures_p2_vs_form_factor(UV, A, p2_values, sweep_labels, q_largest_b)
    area = abs(det(A));

    [~, XX, YY] = change_basis(UV, A);
    
    form_factors = form_factor(q_largest_b, A, area);
    
    figure_latex(1);

    for sweep_idx = 1:size(sweep_labels, 2)
        plot(p2_values, form_factors(:, sweep_idx),...
            "DisplayName", sweep_labels(sweep_idx),...
            "LineWidth", 2);
        hold on;
    end

    hold off;
    xlabel("$p_2$"); ylabel("$\chi$");
    legend("Location", "northwest");
    
    q_s = q_largest_b(:, :, end, 1);

    figure_latex(2);
    surf(XX, YY, q_s .* conj(q_s));
    xlabel("$x$"); ylabel("$y$"); zlabel("$\left|\psi\right|^2$");
    surf_style(gca());
    
    q_e = q_largest_b(:, :, end, end);

    figure_latex(3);
    surf(XX, YY, q_e .* conj(q_e));
    xlabel("$x$"); ylabel("$y$"); zlabel("$\left|\psi\right|^2$");
    surf_style(gca());
end
