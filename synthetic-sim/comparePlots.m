function comparePlots(path, tag1, tag2, printplot)
% comparePlots: Overlay DIU and DDU simulation results for all units
%
% Inputs:
%   path       - folder where results are saved (e.g. './results/')
%   tag1,tag2  - framework names ('diu','ddu')
%   printplot  - true/false to save figures

    if nargin < 4, printplot = false; end
    font = 10; xfont = 8;

    % Load DIU and DDU files
    D1 = load(fullfile(path, sprintf('results_unit1_%s.mat', lower(tag1))));
    D2 = load(fullfile(path, sprintf('results_unit1_%s.mat', lower(tag2))));

    n_units = numel(D1.sysparams);

    for i = 1:n_units
        sp = D1.sysparams(i);
        base = (i-1)*5;

        % Extract trajectories
        V1 = D1.X(:, base+1); p1 = D1.X(:, base+2);
        u1 = D1.X(:, base+3); s1 = D1.X(:, base+4); q1 = D1.X(:, base+5);

        V2 = D2.X(:, base+1); p2 = D2.X(:, base+2);
        u2 = D2.X(:, base+3); s2 = D2.X(:, base+4); q2 = D2.X(:, base+5);

        head1 = sp.a .* (V1.^sp.b);
        head2 = sp.a .* (V2.^sp.b);

        % Plot overlays
        fig = figure('Position',[100 100 1200 600]);
        sgtitle(sprintf('DIU vs DDU – Unit %d', sp.unit), ...
            'FontSize', font+2, 'FontWeight','bold');

        subplot(2,3,1)
        plot(u1,'-b','LineWidth',2); hold on;
        plot(u2,'--r','LineWidth',2);
        yline(sp.max_ut,'--k'); yline(sp.min_ut,'--k');
        title('Generation Outflow'); xlabel('Hour'); ylabel('Flow (m^3/hr)');
        legend('DIU','DDU','Location','best');

        subplot(2,3,2)
        plot(p1,'-b','LineWidth',2); hold on;
        plot(p2,'--r','LineWidth',2);
        yline(sp.F,'--k');
        title('Hydropower Generation'); xlabel('Hour'); ylabel('MWh');

        subplot(2,3,3)
        plot(s1,'-b','LineWidth',2); hold on;
        plot(s2,'--r','LineWidth',2);
        title('Spill Outflow'); xlabel('Hour'); ylabel('Flow (m^3/hr)');

        subplot(2,3,4)
        plot(V1,'-b','LineWidth',2); hold on;
        plot(V2,'--r','LineWidth',2);
        yline(sp.max_V,'--k'); yline(sp.min_V,'--k');
        title('Reservoir Volume'); xlabel('Hour'); ylabel('m^3');

        subplot(2,3,5)
        plot(q1,'-b','LineWidth',2); hold on;
        plot(q2,'--r','LineWidth',2);
        title('Inflow'); xlabel('Hour'); ylabel('Flow (m^3/hr)');

        subplot(2,3,6)
        plot(head1,'-b','LineWidth',2); hold on;
        plot(head2,'--r','LineWidth',2);
        yline(sp.max_h,'--k'); yline(sp.min_h,'--k');
        title('Forebay Elevation'); xlabel('Hour'); ylabel('Elevation (m)');

        % Save figure if requested
        if printplot
            outFile = fullfile(path, sprintf('Unit%d_DIUvsDDU.png', sp.unit));
            saveas(fig, outFile);
        end
    end
end
