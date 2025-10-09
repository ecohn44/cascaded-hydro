function simPlots(path, X, q, q_hist, sysparams, T, c, lag, printplot)
    % simPlots: Create one figure per unit and save as PNG
    % X columns: 1=V1, 2=p1, 3=u1, ...

    font = 10;
    xfont = 8;

    n_units = numel(sysparams);

    for i = 1:n_units
        sp = sysparams(i);
        base = (i-1)*3; % offset into X matrix

        % Extract decision variables for this unit
        V = X(:, base+1);
        p = X(:, base+2);
        u = X(:, base+3);
        % Extract streamflow time series for this unit 
        q_pred = q(:, i);
        q_hist = q_hist((1+lag):end); % temp: bonneville as ref

        % Compute head and max power
        head = sp.a .* (V.^sp.b);
        p_max = (c * u .* head);

        % ---- Create new figure for this unit ----
        fig = figure('Position',[100 100 1200 600]);

        % Add a banner title
        sgtitle(['Simulation Results – Unit ', num2str(sp.unit)], 'FontSize', font+2, 'FontWeight','bold');

        % Subplot 1: Outflow
        subplot(2,3,1);
        plot(1:T, u, 'LineWidth', 2, 'DisplayName', 'Outflow'); hold on;
        yline(sp.max_ut, '--r', 'LineWidth', 1.5);
        yline(sp.min_ut, '--r', 'LineWidth', 1.5);
        xlabel('Hour','FontSize',xfont); ylabel('Flow (m^3/hr)');
        xlim([1, T]);
        title('Generation Outflow','FontSize',font); 

        % Subplot 2: Generation
        subplot(2,3,2);
        plot(1:T, p, 'LineWidth', 2, 'DisplayName','Generation'); hold on;
        yline(sp.F, '--r', 'LineWidth', 1.5);
        plot(1:T, p_max, '--k', 'DisplayName','P_{max}');
        xlabel('Hour','FontSize',xfont); ylabel('MWh');
        title('Hydropower Generation','FontSize',font); 
        xlim([1, T]);

        % Subplot 4: Volume
        subplot(2,3,4);
        plot(1:T, V, 'LineWidth', 2, 'DisplayName', 'Volume'); hold on;
        % yline(sp.max_V, '--r');
        % yline(sp.min_V, '--r');
        xlabel('Hour','FontSize',xfont); ylabel('m^3');
        title('Reservoir Volume','FontSize',font);
        xlim([1, T]);

        % Subplot 5: Inflow
        subplot(2,3,5);
        plot(1:T, q_pred, 'LineWidth', 2, 'DisplayName', 'Pred Inflow'); hold on
        plot(1:T, q_hist, '--g', 'LineWidth', 2, 'DisplayName', 'Hist Inflow')
        xlabel('Hour','FontSize',xfont); ylabel('Flow (m^3/hr)');
        title('Inflow','FontSize',font);
        xlim([1, T]);

        % Subplot 6: Forebay Elevation 
        subplot(2,3,6)
        plot(1:T, head, 'LineWidth', 2, 'DisplayName', 'Elevation'); hold on;
        yline(sp.max_h, '--r','LineWidth', 1.5);
        yline(sp.min_h, '--r','LineWidth', 1.5);
        xlabel('Hour','FontSize',xfont); ylabel('Elevation (m)');
        title('Forebay Elevation','FontSize',font); 
        xlim([1, T]);

        if printplot
            % ---- Save figure with unit number ----
            filename = fullfile(path, sprintf("Unit%d_results.png", sp.unit));
            saveas(fig, filename);
            % close(fig); % close to avoid clutter
        end
    end
end
