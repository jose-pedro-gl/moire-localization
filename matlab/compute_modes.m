function [b, q] = compute_modes(A, V, bc, Nb, Nq, pde_impl, cache_path, export_path)
    A_str = string(mat2str(A, 6));
    V_str = string(mat2str(V, 6));

    if bc.type == "zero"
        bc_str = "zero";
    elseif bc.type == "bloch"
        bc_str = "bloch;k=" + string(mat2str(bc.k, 6));
    end

    pde_impl_str = func2str(pde_impl);

    hash = hash_string(strjoin([A_str, V_str, bc_str, pde_impl_str], ";"));

    % delete(cache_path + hash + ".mat");

    if isempty(cache_path)
        Nb_saved = 0;
        Nq_saved = 0;
    else
        [b, q, Nb_saved, Nq_saved] = load_modes(cache_path, hash, Nb, Nq);
    end

    if Nb_saved < Nb || Nq_saved < Nq
        [b, q] = pde_impl(A, V, bc, max(Nb, Nq));

        if ~isempty(cache_path)
            save_modes(cache_path, hash, b(1:max(Nb_saved, Nb)), q(:, :, 1:max(Nq_saved, Nq)));
        end

        b = b(1:Nb);
        q = q(:, :, 1:Nq);
    end

    if ~isempty(export_path)
        export_modes(cache_path, export_path, hash, Nb, Nq);
    end
end
