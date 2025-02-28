function [A1, A2, A] = square_moire(L, n1, n2)
    [A1, A1_complex] = square(L);

    zm1 = dot([n1, n2], A1_complex);
    zm2 = dot([-n2, n1], A1_complex);

    r = (-2 * n1 * n2 + j * (n1^2 - n2^2)) / (n1^2 + n2^2);
    %r = ((n1^2 - n2^2) + j * 2 * n1 * n2) / (n1^2 + n2^2);

    A2 = comp2vec(r * A1_complex);
    A = comp2vec([zm1, zm2]);
end
