function m = vec2idx(n, s)
    if (size(s, 2) == 0)
        m = 0;
        return;
    end

    m = n(1) + s(1) * vec2idx(n(2:end), s(2:end));
end
