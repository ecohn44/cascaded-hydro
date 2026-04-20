%% Plot Driver for Risk Allocation Convergence
% ========================================================================
% Author: Eliza Cohn
% Description: Plot the risk allocation convergence for SSH to BON
% ========================================================================

close all;

% Set default font to match other figures
set(groot,'defaultAxesFontName','Helvetica');
set(groot,'defaultTextFontName','Helvetica');

% Parameters
eps = 0.05;
T   = 20;
n   = 3;
tt  = 1:T;

t1 = 5;
t2 = 10;
t3 = 15;

% Font / style settings
font_axis = 18;    % axis labels
font_tick = 16;    % tick labels
font_leg  = 18;    % legend
lw_main   = 2.5;
lw_ref    = 2.0;

% Initialize alpha vals 
alpha_vals = (1/n) * ones(T,n);

% Set up allocation values
alpha_blocks = [
    1      0      0
    1/2    1/2    0
    1/n    1/n    1/n
];

row_starts = [t1, t2, t3];
row_ends   = [t2-1, t3-1, T];

for k = 1:length(row_starts)
    rows = row_ends(k) - row_starts(k) + 1;
    alpha_vals(row_starts(k):row_ends(k), :) = repmat(alpha_blocks(k,:), rows, 1);
end

% Plot Figure 
figure('Name','SSH Diagnostics','NumberTitle','off', ...
       'Position',[100 100 900 500]);   % make figure wider/taller like paper figs

plot(tt, alpha_vals(:,1)*eps, '-',  'LineWidth', lw_main); hold on
plot(tt, alpha_vals(:,2)*eps, '--', 'LineWidth', lw_main);
plot(tt, alpha_vals(:,3)*eps, '-.', 'LineWidth', lw_main);
yline(eps/n, 'k:', 'LineWidth', lw_ref);

xlabel('Time (hours)', 'FontSize', font_axis, 'FontName','Helvetica');
ylabel('Allocated Risk', 'FontSize', font_axis, 'FontName','Helvetica'); 

xlim([1, T]);
ylim([-0.0025, 0.055]);

% Create a cell array of strings for the legend labels
legendLabels = cell(1, n+1);
for i = 1:n
    legendLabels{i} = sprintf('$\\epsilon_{%d}$', i);
end
legendLabels{n+1} = '$\frac{\epsilon}{n}$';

lg = legend(legendLabels, ...
   'Location', 'southoutside', ...  
   'Orientation', 'horizontal', ...  %'NumColumns', 1, ...                
   'Box', 'off', ...
   'FontSize', font_leg, ...
   'Interpreter','latex');         

yticks([0, 0.01, 0.02, 0.03, 0.04, 0.05]);

grid on
set(gca, ...
    'YGrid', 'on', ...
    'XGrid', 'off', ...
    'GridAlpha', 0.15, ...
    'GridColor', [0 0 0], ...
    'FontSize', font_tick, ...
    'FontName', 'Helvetica', ...
    'LineWidth', 1.2, ...
    'Box', 'on');

