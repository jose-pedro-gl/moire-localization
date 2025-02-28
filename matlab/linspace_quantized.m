function X = linspace_quantized(x1, x2, y, N, M)
    a = x2 - x1;
    n = max(0, floor(y / a * 2^N) - 2^(M - 1));
    X = (n:(n + 2^M)) / (2^N) * a;
end
