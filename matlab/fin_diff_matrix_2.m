function M = fin_diff_matrix_2(A, V, bc)
    [Nv, Nu] = size(V);
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

    M = zeros(Nv, Nu, Nv, Nu);

    for vn = 1:Nv
        for un = 1:Nu
            M = fill(M, Nv, Nu, vn, un, vn, un, bc, -(opu + opv) + V(vn, un));

            M = fill(M, Nv, Nu, vn, un, vn, un + 1, bc, opu2);
            M = fill(M, Nv, Nu, vn, un, vn, un - 1, bc, opu2);

            M = fill(M, Nv, Nu, vn, un, vn + 1, un, bc, opv2);
            M = fill(M, Nv, Nu, vn, un, vn - 1, un, bc, opv2);

            M = fill(M, Nv, Nu, vn, un, vn + 1, un + 1, bc, opuv2);
            M = fill(M, Nv, Nu, vn, un, vn - 1, un + 1, bc, -opuv2);
            M = fill(M, Nv, Nu, vn, un, vn + 1, un - 1, bc, -opuv2);
            M = fill(M, Nv, Nu, vn, un, vn - 1, un - 1, bc, opuv2);
        end
    end

    M = reshape(M, [N, N]);
end

function M = fill(M, Nv, Nu, vn, un, vm, um, bc, value)
    e = exp(1);

    [um_q, um_r] = idiv(um - 1, Nu);
    [vm_q, vm_r] = idiv(vm - 1, Nv);

    if bc.type == "zero"
        op = (um_q == 0) .* (vm_q == 0);
    elseif bc.type == "bloch"
        op = e^(j * 2 * pi * bc.k * [um_q; vm_q]);
    end

    M(vn, un, vm_r + 1, um_r + 1) = value * op;
end
