function [b, q] = fin_diff_fn_2(A, V, bc, Nb, eigs_impl)
    if bc.type ~= "bloch" && bc.type ~= "zero"
        ex = MException(...
            "fin_diff_fn_2:invalid_bc",...
            "Finite differences method requires either zero or Bloch" +...
                " boundary conditions");

        throw(ex);
    end

    [Nv, Nu] = size(V);

    M = fin_diff_matrix_sparse_2(A, V, bc);

    [q, b] = eigs_impl(M, Nb);

    q = reshape(q, [Nv, Nu, Nb]);
end
