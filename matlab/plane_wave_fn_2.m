function [b, q] = plane_wave_fn_2(A, V, bc, Nb, eigs_impl)
    if bc.type ~= "bloch"
        ex = MException(...
            "plane_wave_fn_2:invalid_bc",...
            "Plane wave method requires Bloch boundary conditions");

        throw(ex);
    end

    e = exp(1);

    [Nv, Nu] = size(V);
    N = Nv * Nu;

    % VF = 1 / N * fftshift(fft2(V));
    VF = fft_series_approx_2(V, Nu, Nv, 100);

    M = plane_wave_matrix_sparse_2(A, VF, bc);

    [Q, b] = eigs_impl(M, Nb);

    % a normalization is convenient because:
    % - eigenvectors are normalized to unitary L2 norm
    % - these are Fourier series coefficients
    % - to compare to other methods, we want unitary vectors in space coordinates
    %
    % so, assuming q is the vector in space coordinates, and Q is the vector in
    % Fourier coordinates, we know, according to Parseval's theorem, that:
    %
    % |q|^2 = |Q|^2 / N
    %
    % taking |q| = 1, as reference:
    %
    % |Q|^2 = N <=> |Q| = sqrt(N)
    %
    Q = sqrt(N) * Q;
    Q = reshape(Q, [Nv, Nu, Nb]);

    [UU, VV] = meshgrid(linspace_ex(0, 1, Nu), linspace_ex(0, 1, Nv));
    UV = cat(3, UU, VV);
    KK = tensorprod_shim(UV, 2 * pi * bc.k, 3, 2);

    Q = ifftshift(ifftshift(Q, 1), 2);
    q = ifft2(Q) .* e.^(j * KK);
end
