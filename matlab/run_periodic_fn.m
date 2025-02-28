function [b, q] = run_periodic_fn(A, V, K, Nb, method, eigs_impl, cache_path, export_path)
    % method
    % 1: finite differences; 2: plane wave

    e = exp(1);

    % compat
    Nq = max(1, -Nb);
    Nb = abs(Nb);

    V = V(1:(end - 1), 1:(end - 1));

    if method == 1
        pde_impl = @(A, V, bc, Nb) fin_diff_fn_2(A, V, bc, Nb, eigs_impl);
    else
        pde_impl = @(A, V, bc, Nb) plane_wave_fn_2(A, V, bc, Nb, eigs_impl);
    end

    start = tic;
    [b, q] = compute_bands(A, V, K, Nb, Nq, pde_impl, cache_path, export_path);
    toc(start)

    %q_size = size(q);
    % q_old = q(:, :, :);
    % q = zeros(size(q_old) + [1, 1, 0]);
    % 
    % for idx = 1:size(q_old, 3)
    %     q(1:q_size(1), 1:q_size(2), idx) = q_old(:, :, idx);
    %     q(:, (q_size(2) + 1), idx) = q(:, 1, idx) .* e.^(j * 2 * pi * K(idx, 1));
    %     q((q_size(1) + 1), :, idx) = q(1, :, idx) .* e.^(j * 2 * pi * K(idx, 2));
    % end

    %q = reshape(q, [size(q, 1), size(q, 2), q_size(3:end)]);

    q_size = size(q);
    Kx = reshape(K(:, 1), [1, 1, 1, size(K, 1)]);
    q_ex1 = reshape(q(:, 1, :), [q_size(1), 1, q_size(3:end)]);
    q = cat(2, q, q_ex1 .* e.^(j * 2 * pi * Kx));

    q_size = size(q);
    Ky = reshape(K(:, 2), [1, 1, 1, size(K, 1)]);
    q_ex2 = reshape(q(1, :, :), [1, q_size(2), q_size(3:end)]);
    q = cat(1, q, q_ex2 .* e.^(j * 2 * pi * Ky));
end
