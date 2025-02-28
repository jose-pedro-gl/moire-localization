function [q, r] = idiv(a, b)
    r = mod(a, b);
    q = (a - r) / b;
end
