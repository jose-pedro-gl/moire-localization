function V = V_cut(UU, VV, mu, mv, c, V)
    V = ((mu * UU + mv * VV) > c) .* V;
end
