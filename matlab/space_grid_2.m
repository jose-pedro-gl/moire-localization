function [UV, UU, VV] = space_grid_2(TU, TV)
    [UU, VV] = meshgrid(TU, TV);
    UV = cat(3, UU, VV);
end
