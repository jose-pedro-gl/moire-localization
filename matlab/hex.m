function [A, A_complex] = hex(L)
    e = exp(1);

    A_complex = L * [1, e^(j * pi / 3)];
    A = comp2vec(A_complex);
end
