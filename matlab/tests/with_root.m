function x = with_root(f)
    addpath("..");
    x = f();
    rmpath("..");
end
