function V = V_2019(UU1, VV1, UU2, VV2, p1, p2, E0)
    V_I = V_2016(UU1, VV1, UU2, VV2, p1, p2);
    V = -E0 ./ (1 + V_I .* conj(V_I));
end
