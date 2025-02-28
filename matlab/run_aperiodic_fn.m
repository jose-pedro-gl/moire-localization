function [b, q] = run_aperiodic_fn(A, V, Nb, eigs_impl, cache_path, export_path)
    bc = struct("type", "zero");
    pde_impl = @(A, V, bc, Nb) fin_diff_fn_2(A, V, bc, Nb, eigs_impl);

    V = V(2:(end - 1), 2:(end - 1));

    start = tic;
    [b, q] = compute_modes(A, V, bc, Nb, Nb, pde_impl, cache_path, export_path);
    toc(start)

    q_size = size(q);
    q_old = q(:, :, :);
    q = zeros(size(q_old) + [2, 2, zeros(1, size(size(q_old), 2) - 2)]);

    for idx = 1:size(q_old, 3)
        q(2:(q_size(1) + 1), 2:(q_size(2) + 1), idx) = q_old(:, :, idx);
    end

    q = reshape(q, [size(q, 1), size(q, 2), q_size(3:end)]);
end
