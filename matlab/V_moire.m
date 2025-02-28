function V = V_moire(UV, A1, A2, A, V_fn, p1, p2)
    [~, UU1, VV1] = change_basis(UV, A1^(-1) * A);
    [~, UU2, VV2] = change_basis(UV, A2^(-1) * A);

    V = p1 * V_fn(UU1, VV1) + p2 * V_fn(UU2, VV2);
end
