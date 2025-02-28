function file = create_modes_matfile(path, hash)
    results_path = path + hash + ".mat";

    [~, ~] = mkdir(path);
    file = matfile(results_path, "Writable", true);
end

