function [Q, b] = eigs_matlab(M, N, sigma, tol)
    % eigs has a bug:
    %
    % apparently, when there are eigenvalues with a large
    % multiplicity, some instances are skipped
    %
    % (this is likely related to the fact that the algorithm
    % assumes a gap between eigenvalues)
    %
    % computing twice the eigenvalues apparently works
    if N == 1
        N_safe = 1;
    else
        N_safe = min(2 * N, size(M, 1));
    end

    subspace_dim = min(max(2 * N_safe, 20), size(M, 1));
    max_iterations = 300;

    flag = 1;

    while flag
        options = struct();
        options.tol = tol;
        options.maxit = max_iterations;

        if subspace_dim <= size(M, 1)
            options.p = subspace_dim;
        end

        [Q, B, flag] = eigs(M, N_safe, sigma, options);

        b = diag(B);

        if flag
            old_subspace_dim = subspace_dim;
            old_max_iterations = max_iterations;

            subspace_dim = min(2 * subspace_dim, size(M, 1));
            max_iterations = max_iterations + 20;

            warning("eigs failed to converge (%d/%d converged) with " + ...
                "MaxIterations = %d, SubspaceDimension = %d" + ...
                newline() + ...
                "retrying with " + ...
                "MaxIterations = %d, SubspaceDimension = %d", ...
                (N_safe - sum(isnan(b))), N_safe, ...
                old_max_iterations, old_subspace_dim, ...
                max_iterations, subspace_dim);
        end

        Q = Q(:, 1:N);
        b = b(1:N);
    end
end
