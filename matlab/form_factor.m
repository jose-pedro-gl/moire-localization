function form_factor = form_factor(q, A, area)
    q_size = size(q);

    dv = 1 / (q_size(1) - 1);
    du = 1 / (q_size(2) - 1);

    abs_q_sq = q .* conj(q);

    I4 = trapz(du, trapz(dv, abs_q_sq.^2));
    I2 = trapz(du, trapz(dv, abs_q_sq));

    % cell_area = abs(det(A));
    % N = area / cell_area;

    form_factor = sqrt(1 / area * I4) ./ I2;
    form_factor = reshape(form_factor, [q_size(3:end), 1, 1]);
end
