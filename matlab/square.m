function [A, A_complex] = square(L)
    A_complex = L * [1, j];
    A = comp2vec(A_complex);
end
