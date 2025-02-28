function [b, q, Nb_saved, Nq_saved] = load_modes(path, hash, Nb, Nq)
    results_path = path + hash + ".mat";

    if ~isfile(results_path)
        b = [];
        q = [];
        Nb_saved = 0;
        Nq_saved = 0;
    else
        file = matfile(results_path);

        Nb_saved = size(file.b, 1);
        b = file.b(1:min(Nb, Nb_saved), 1);

        if isempty(who(file, "q"))
            Nq_saved = size(file.q_largest_b, 3);
            q = file.q_largest_b(:, :, 1:min(Nq, Nq_saved));
        else
            Nq_saved = size(file.q, 3);
            q = file.q(:, :, 1:min(Nq, Nq_saved));
        end
    end
end
