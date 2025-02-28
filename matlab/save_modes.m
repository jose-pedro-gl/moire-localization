function save_modes(path, hash, b, q)
    file = create_modes_matfile(path, hash);

    file.b = b;
    file.q = q;
end
