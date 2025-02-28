function path = build_path(points, N)
    [Np, Dp] = size(points);
    a = linspace_ex(0, 1, N);
    A = cat(3, 1 - a, a);
    start_end = cat(3, points(1:(end - 1), :), points(2:end, :));
    path = reshape(tensorprod_shim(A, start_end, 3, 3), [(Np - 1) * N, Dp]);
    path(end + 1, :) = points(end, :);
end
