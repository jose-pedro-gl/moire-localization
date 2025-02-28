function [A1, A2, A] = hex_moire(L, n1, n2)
    % a = n2^2 - n1^2;
    % b = -2 * n1 * n2 - n2^2;
    % c = n1^2 + n1 * n2 + n2^2;

    e = exp(1);

    [A1, A1_complex] = hex(L);

    zm1 = dot([n1, n2], A1_complex);
    zm2 = zm1 * e^(j * pi/3);

    r = (n2^2 - n1^2 + e^(j * pi / 3) * (-2 * n1 * n2 - n2^2)) / (n1^2 + n1 * n2 + n2^2);

    A2 = comp2vec(r * A1_complex);
    A = comp2vec([zm1, zm2]);
end
