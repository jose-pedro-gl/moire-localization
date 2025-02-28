function v = comp2mat(c)
    v = [ real(c), -imag(c); imag(c), real(c) ];
end
