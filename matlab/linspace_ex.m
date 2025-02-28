function X = linspace_ex(x1, x2, n)
    X = linspace(x1, x2 - (x2 - x1) / n, n);
end
