function M = plane_wave_matrix_sparse_2(A, VF, bc)
    [Nv, Nu] = size(VF);

    Nu2 = (Nu - 1 - mod(Nu - 1, 2)) / 2;
    Nv2 = (Nv - 1 - mod(Nv - 1, 2)) / 2;

    Nu = 2 * Nu2 + 1;
    Nv = 2 * Nv2 + 1;

    N = Nv * Nu;

    BT = 2 * pi * A^(-1);

    M = sparse(N, N);
    [VF_R, VF_C] = find(VF);

    for vn = 1:Nv
        for un = 1:Nu
            gv = vn - 1 - Nv2;
            gu = un - 1 - Nu2;

            kxy = (bc.k + [gu, gv]) * BT;

            m = vec2idx([vn - 1, un - 1], [Nv, Nu]);

            for idx = 1:length(VF_R)
                tn = VF_R(idx);
                sn = VF_C(idx);

                if tn >= max(1, vn - Nv2) && tn <= min(2 * Nv2 + 1, vn + Nv2) ...
                    && sn >= max(1, un - Nu2) && sn <= min(2 * Nu2 + 1, un + Nu2) ...
                    && abs(VF(tn, sn)) > 1e-9
                    l = vec2idx([vn - tn + Nv2, un - sn + Nu2], [Nv, Nu]);
                    M(m + 1, l + 1) = VF(tn, sn);
                end
            end

            M(m + 1, m + 1) = M(m + 1, m + 1) - 1 / 2 * kxy * kxy.';
        end
    end
end
