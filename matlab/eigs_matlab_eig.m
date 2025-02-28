function [Q, b] = eigs_matlab_eig(M, N)
    [Q, b] = eig(M, "vector");

    Q = Q(:, 1:N);
    b = b(1:N);

    [b, eig_idx] = sort(b, "descend");
    Q = Q(:, eig_idx);
end
