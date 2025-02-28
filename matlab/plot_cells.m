function [] = plot_cells(A, N1, N2)
    rU = N1(1):1:N2(1);
    rV = N1(2):1:N2(2);
    [UU, VV] = meshgrid(rU, rV);
    UV = cat(3, UU, VV);
    [~, XX, YY] = change_basis(UV, A);
    XX = repelem(XX(:), 2, 1);
    YY = repelem(YY(:), 2, 1);
    DD = repmat(A.', size(XX, 1) / 2, 1);
    quiver(XX, YY, DD(:, 1), DD(:, 2), "off",...
        "LineWidth", 0.5, "ShowArrowHead", false);
end
