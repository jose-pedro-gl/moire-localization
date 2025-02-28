function Z = tensorprod_shim(X, Y, from, to)
    if exist("tensorprod", "builtin") ~= 0
        Z = tensorprod(X, Y, from, to);
        return;
    end

    sX = size(X);
    sY = size(Y);

    X = permute(X, [from, 1:(from - 1), (from + 1):length(sX)]);
    Y = permute(Y, [to, 1:(to - 1), (to + 1):length(sY)]);

    sX = size(X);
    sY = size(Y);

    N = sX(1);

    if (N ~= sY(1))
        throw MException();
    end

    sX2 = sX(2:end);
    sY2 = sY(2:end);

    sX2_lin = prod(sX2);
    sY2_lin = prod(sY2);

    sZ = [sX2, sY2];
    Z = zeros(sX2_lin, sY2_lin);

    for nx = 1:sX2_lin
        for ny = 1:sY2_lin
            for n = 1:N
                Z(nx, ny) = Z(nx, ny) + X(n, nx) * Y(n, ny);
            end
        end
    end

    Z = reshape(Z, sZ);
end
