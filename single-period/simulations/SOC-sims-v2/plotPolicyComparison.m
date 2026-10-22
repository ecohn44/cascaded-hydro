function plotPolicyComparison(IVI_bon, IVI_ssh, energy_bon, energy_ssh, alpha_set, baseline_set)

    %{
    IVI_bon_plot    = IVI_bon;    IVI_bon_plot(isnan(IVI_bon_plot))       = 0;
    IVI_ssh_plot    = IVI_ssh;    IVI_ssh_plot(isnan(IVI_ssh_plot))       = 0;
    %}
    energy_bon_plot = energy_bon; energy_bon_plot(isnan(energy_bon_plot)) = 0;
    energy_ssh_plot = energy_ssh; energy_ssh_plot(isnan(energy_ssh_plot)) = 0;

    xLabels = arrayfun(@(x) sprintf('%.2f',x), alpha_set,    'UniformOutput', false);
    yLabels = arrayfun(@(x) sprintf('%.3f',x), baseline_set, 'UniformOutput', false);

    figure;
    %{
    tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

    nexttile;
    b1 = bar3(IVI_bon_plot); hold on;
    for k = 1:length(b1)
        b1(k).FaceColor = 'flat';
        b1(k).CData = b1(k).ZData;
    end
    b2 = bar3(IVI_ssh_plot);
    for k = 1:length(b2)
        b2(k).FaceColor = 'flat';
        b2(k).CData = b2(k).ZData;
        b2(k).FaceAlpha = 0.4;
    end
    set(gca, 'XTickLabel', xLabels, 'YTickLabel', yLabels);
    xlabel('Drought Intensity (\alpha)');
    ylabel('Baseline Flow (q_0)');
    zlabel('IVI (m^3)');
    title('Constraint Violation Severity');
    colormap parula; colorbar;
    legend([b1(1) b2(1)], 'BON','SSH','Location','northwest');

    nexttile;
    %}
    b3 = bar3(energy_bon_plot/1e3); hold on;
    for k = 1:length(b3)
        b3(k).FaceColor = 'flat';
        b3(k).CData = b3(k).ZData;
    end
    b4 = bar3(energy_ssh_plot/1e3);
    for k = 1:length(b4)
        b4(k).FaceColor = 'flat';
        b4(k).CData = b4(k).ZData;
        b4(k).FaceAlpha = 0.4;
    end
    set(gca, 'XTickLabel', xLabels, 'YTickLabel', yLabels);
    xlabel('Drought Intensity (\alpha)');
    ylabel('Baseline Flow (q_0)');
    zlabel('Total Energy (GWh)');
    colormap parula; %colorbar;
    %legend([b3(1) b4(1)], 'BON','SSH','Location','northeast');

    bon_color = [0.2 0.5 0.8];  % may return 'flat', so hardcode if needed
    ssh_color = [0.85 0.85 0.85];
    
    h_bon = patch(nan, nan, nan, 'FaceColor', bon_color,   'EdgeColor', 'k');
    h_ssh = patch(nan, nan, nan, 'FaceColor', ssh_color, 'EdgeColor', 'k', 'FaceAlpha', 0.4);
    
    legend([h_bon h_ssh], 'BON', 'SSH', 'Location', 'northwest');
    

end