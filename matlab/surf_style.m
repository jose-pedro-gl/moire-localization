function surf_style(ax)
    maxX = max(abs(ax.Children(1).XData), [], "all");
    maxY = max(abs(ax.Children(1).YData), [], "all");
    maxZ = max(abs(ax.Children(1).ZData), [], "all");
    maxXY = sqrt(maxX * maxY);

    turbo_cm_file = matfile("turbo_cm.mat");

    axis(ax, "tight");
    colorbar(ax);
    colormap(ax, turbo_cm_file.turbo_cm);
    daspect(ax, [1, 1, maxZ / maxXY]);
    shading(ax, "interp");
    view(ax, 0, 90);
end
