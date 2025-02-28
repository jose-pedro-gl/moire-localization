function show_figures_periodic(UV, A, V, K, K_ticks, K_tick_labels, bands, Nb2, q_largest_b, area)
    [~, XX, YY] = change_basis(UV, A);

    figure_latex(1);
    surf(XX, YY, V);
    xlabel("$x$"); ylabel("$y$"); zlabel("$V$");
    title("$V$");
    surf_style(gca());

    K_dist = sqrt(sum((K(2:end, :) - K(1:(end - 1), :)).^2, 2));
    K_plot = cumsum([0; K_dist]);

    figure_latex(2);
    hold off;
    
    for b = 1:Nb2
        plot(K_plot, bands(b, :), "LineWidth", 2);
        xticks(K_plot(K_ticks));
        xticklabels(K_tick_labels);
        hold on;
    end
    
    hold off;
    
    xlabel("$k$"); ylabel("$\beta$");
    xlim("tight");
    
    form_factors = form_factor(q_largest_b, A, area);
    
    figure_latex(3);
    plot(1:size(K, 1), form_factors, "LineWidth", 2);
    xticks(K_ticks); xticklabels(K_tick_labels);
    xlabel("$k$"); ylabel("$\chi$");
    xlim("tight");
    
    q = q_largest_b(:, :, 1);

    figure_latex(4);
    surf(XX, YY, q .* conj(q));
    xlabel("$x$"); ylabel("$y$"); zlabel("$\left|\psi\right|^2$");
    title("$\left|\psi\right|^2$");
    surf_style(gca());
    
    figure_latex(5);
    hist_bands = bands(:);
    hist_bands = sort(hist_bands, "descend");
    histogram(hist_bands, size(bands, 1) * 10, "Normalization", "pdf");
    xlabel("$\beta$"); ylabel("pdf");
end
