function [bands, qs] = compute_bands(A, V, K, Nb, Nq, pde_impl, cache_path, export_path)
    [Nv, Nu] = size(V);
    Nk = size(K, 1);

    bands = zeros(Nb, Nk);
    qs = zeros(Nv, Nu, Nq, Nk);

    for k_idx = 1:Nk
        k = K(k_idx, :);
        bc = struct("type", "bloch", "k", k);

        start = tic;
        [b, q] = compute_modes(A, V, bc, Nb, Nq, pde_impl, cache_path, export_path);

        disp("k = " + mat2str(k) + ":");
        toc(start);

        bands(:, k_idx) = b;
        qs(:, :, :, k_idx) = q;
    end
end
