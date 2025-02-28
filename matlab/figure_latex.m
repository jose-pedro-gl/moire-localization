function f = figure_latex(n)
    f = figure(n);
    clf(f, "reset");

    set(f, "defaultAxesTickLabelInterpreter", "latex");
    set(f, "defaultLegendInterpreter", "latex");
    set(f, "defaultTextInterpreter", "latex");

    set(f, "defaultAxesFontSize", 18);
end
