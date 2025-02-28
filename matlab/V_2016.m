function V = V_2016(UU1, VV1, UU2, VV2, p1, p2)
    V1 = cos(2 * pi * UU1) + cos(2 * pi * VV1);
    V2 = cos(2 * pi * UU2) + cos(2 * pi * VV2);

    V = p1 * V1 + p2 * V2;
end
