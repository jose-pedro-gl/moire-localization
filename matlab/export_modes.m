function export_modes(source_path, target_path, hash, Nb, Nq)
    source = matfile(source_path + hash + ".mat");
    target = create_modes_matfile(target_path, hash);

    if isempty(who(target, "b"))
        target.b = zeros(0, 1);
    end

    if isempty(who(target, "q"))
        target.q = zeros(0, 0, 0);
    end

    if Nb > size(target.b, 1)
        target.b = source.b(1:Nb, 1);
    end

    if Nq > size(target.q, 3)
        target.q = source.q(:, :, 1:Nq);
    end
end
