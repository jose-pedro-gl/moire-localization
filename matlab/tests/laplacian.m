function laplacian = laplacian(A, k, x, f)
    e = exp(1);

    [Nv, Nu] = size(x);

    A_inv = A^(-1);
    A2 = A_inv * A_inv.';

    opufw = f * e^(j * 2 * pi * k(1));
    opubw = f * e^(-j * 2 * pi * k(1));
    opvfw = f * e^(j * 2 * pi * k(2));
    opvbw = f * e^(-j * 2 * pi * k(2));

    vm = 1:Nv;
    v0 = 2:(Nv + 1);
    vp = 3:(Nv + 2);

    um = 1:Nu;
    u0 = 2:(Nu + 1);
    up = 3:(Nu + 2);

    y = zeros(Nv + 2, Nu + 2);
    y(v0, u0) = x;
    y(1, u0) = opvbw * y(Nv + 1, u0);
    y(Nv + 2, u0) = opvfw * y(2, u0);
    y(:, 1) = opubw * y(:, Nu + 1);
    y(:, Nu + 2) = opufw * y(:, 2);

    d2y_du2 = Nu^2 * (y(v0, up) - 2 * y(v0, u0) + y(v0, um));
    d2y_dv2 = Nv^2 * (y(vp, u0) - 2 * y(v0, u0) + y(vm, u0));
    d2y_duv = Nv * Nu / 4 * (y(vp, up) - y(vm, up) - y(vp, um) + y(vm, um));

    laplacian = A2(1, 1) * d2y_du2 + (A2(1, 2) + A2(2, 1)) * d2y_duv + A2(2, 2) * d2y_dv2;
end

