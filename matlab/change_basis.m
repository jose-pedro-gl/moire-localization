function [to, to1, to2] = change_basis(from, M)
    if exist("tensorprod", "builtin") ~= 0
        to = tensorprod_shim(from, M, 3, 2);
    else
        to = zeros(size(from));

        for n = 1:size(M, 2)
            for m = 1:size(M, 1)
                to(:, :, m) = to(:, :, m) + M(m, n) * from(:, :, n);
            end
        end
    end
    
    to1 = to(:, :, 1);
    to2 = to(:, :, 2);
end
