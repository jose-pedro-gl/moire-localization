function q = run_p2_vs_form_factor_fn(A, V_fn, p2_values, sweep_values, eigs_impl, cache_path, export_path)
    [Nv, Nu] = size(V_fn(sweep_values(1), p2_values(1)));

    p2_length = size(p2_values, 2);
    sweep_length = size(sweep_values, 2);

    q = zeros(Nv, Nu, size(p2_values, 2), size(sweep_values, 2));

    for sweep_idx = 1:sweep_length
        sweep = sweep_values(sweep_idx);
    
        for p2_idx = 1:p2_length
            p2 = p2_values(p2_idx);

            V = V_fn(sweep, p2);

            [~, q_largest_b] = run_aperiodic_fn(A, V, 1, eigs_impl, cache_path, export_path);

            q(:, :, p2_idx, sweep_idx) = q_largest_b;
        end
    end
end
