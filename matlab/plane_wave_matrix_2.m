function M = plane_wave_matrix_2(A, VF, bc)
    [Nv, Nu] = size(VF);

    Nu2 = (Nu - 1 - mod(Nu - 1, 2)) / 2;
    Nv2 = (Nv - 1 - mod(Nv - 1, 2)) / 2;

    Nu = 2 * Nu2 + 1;
    Nv = 2 * Nv2 + 1;

    N = Nv * Nu;

    BT = 2 * pi * A^(-1);

    M = zeros(Nv, Nu, Nv, Nu);

    for vn = 1:Nv
        for un = 1:Nu
            gv = vn - 1 - Nv2;
            gu = un - 1 - Nu2;

            kxy = (bc.k + [gu, gv]) * BT;

            for tn = max(1, vn - Nv2):min(2 * Nv2 + 1, vn + Nv2)
                for sn = max(1, un - Nu2):min(2 * Nu2 + 1, un + Nu2)
                    M(vn, un, vn - tn + Nv2 + 1, un - sn + Nu2 + 1) = VF(tn, sn);
                end
            end

            M(vn, un, vn, un) = M(vn, un, vn, un) - 1 / 2 * kxy * kxy.';
        end
    end

    M = reshape(M, [N, N]);
end
