function V = V_lattice(DUV, mask, V_fn, Nu, Nv)
    V = zeros(Nv, Nu);

    for nv = 1:size(DUV, 1)
        for nu = 1:size(DUV, 2)
            if mask(nv, nu) ~= 0
                V = V + mask(nv, nu) * V_fn(DUV(nv, nu, 1), DUV(nv, nu, 2));
            end
        end
    end
end
