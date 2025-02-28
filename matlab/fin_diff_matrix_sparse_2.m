function M = fin_diff_matrix_sparse_2(A, V, bc)
    e = exp(1);

    [Nv, Nu] = size(V);
    V_size = size(V);
    N = Nv * Nu;

    A_inv = A^(-1);
    A2 = A_inv * A_inv.';

    if bc.type == "bloch"
        du_inv = Nu;
        dv_inv = Nv;
    else
        du_inv = Nu + 1;
        dv_inv = Nv + 1;
    end

    opu = A2(1, 1) * du_inv^2;
    opv = A2(2, 2) * dv_inv^2;
    opu2 = opu / 2;
    opv2 = opv / 2;
    opuv2 = (A2(1, 2) + A2(2, 1)) * dv_inv * du_inv / 4 / 2;

    I = zeros(1, 9 * N);
    J = zeros(1, 9 * N);
    P = zeros(1, 9 * N);

    [du, dv] = meshgrid(-1:1, -1:1);

    du = reshape(du, [1, 9]);
    dv = reshape(dv, [1, 9]);

    x = [
        opuv2,      opv2,       -opuv2
        opu2,   -(opu + opv),   opu2
        -opuv2,     opv2,       opuv2
    ];

    x = reshape(x, [1, 9]);

    for vn = 1:Nv
        for un = 1:Nu
            m = vec2idx([vn - 1, un - 1], V_size);
            r = (9 * m + 1):(9 * m + 9);

            [um_q, um_r] = idiv(un + du - 1, Nu);
            [vm_q, vm_r] = idiv(vn + dv - 1, Nv);

            I(r) = m + 1;
            J(r) = Nv * um_r + vm_r + 1;

            if bc.type == "zero"
                P(r) = x .* (um_q == 0) .* (vm_q == 0);
            elseif bc.type == "bloch"
                P(r) = x .* e.^(j * 2 * pi * bc.k * [um_q; vm_q]);
            end

            P(9 * m + 5) = P(9 * m + 5) + V(vn, un);
        end
    end

    idx = find(abs(P) > 1e-9);

    I = I(idx);
    J = J(idx);
    P = P(idx);

    M = sparse(I, J, P, N, N);
end
