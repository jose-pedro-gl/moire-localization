function V = V_trunc(V, UU, VV, maxU, maxV)
    V = V .* (abs(UU) < maxU) .* (abs(VV) < maxV);
end
