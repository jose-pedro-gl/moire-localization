function V = V_gaussian(UV, Aw)
    e = exp(1);

    [~, UUW, VVW] = change_basis(UV, Aw);
    V = e.^(-(UUW.^2 + VVW.^2));
end

% function V = V_gaussian(UVw, A1, A2, p1, p2, Aw, Nu, Nv)
%     V = p1 * V_gaussian_simple(UVw, A1, Aw, Nu, Nv) +...
%         p2 * V_gaussian_simple(UVw, A2, Aw, Nu, Nv);
% end

function V = V_gaussian_simple(UVw, A, Aw, Nu, Nv)
    e = exp(1);

    V = zeros(size(UVw, [1, 2]));

    for nv = -Nv:Nv
        for nu = -Nu:Nu
            UVN = UVw;
            N = A * [nu; nv];

            if abs(N(1)) > 3.3*5.5 || abs(N(2)) > 3.3*5.5
                continue;
            end

            UVN(:, :, 1) = UVN(:, :, 1) - N(1);
            UVN(:, :, 2) = UVN(:, :, 2) - N(2);

            [~, UUW, VVW] = change_basis(UVN, Aw);
            V = V + e.^(-(UUW.^2 + VVW.^2));
        end
    end
end
