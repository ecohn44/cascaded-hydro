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
font_axis = 20;
font_tick = 21;
font_leg  = 24;
lw_main   = 3.5;
lw_ref    = 1.5;

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
       'Position',[100 100 900 500], ...
       'Color','w');

plot(tt, alpha_vals(:,1)*eps, '-',  'LineWidth', lw_main); hold on
plot(tt, alpha_vals(:,2)*eps, '--', 'LineWidth', lw_main);
plot(tt, alpha_vals(:,3)*eps, '-.', 'LineWidth', lw_main);

yline(eps/n, ':', ...
    'Color', [0.35 0.35 0.35], ...
    'LineWidth', lw_ref);

xlabel('Time (hours)', 'FontSize', font_axis, 'FontName','Helvetica');
ylabel('Allocated Risk', 'FontSize', font_axis, 'FontName','Helvetica'); 

xlim([1, T]);
ylim([-0.0025, 0.052]);

% Symbolic y-axis labels
yticks([0, eps/n, eps/2, eps])

%{
yticklabels({ ...
    '0', '$\frac{\varepsilon}{6}$', ...
    '$\frac{\varepsilon}{4}$', ...
    '$\frac{\varepsilon}{2}$'})
%}

yticklabels({ ...
    '0', '$\frac{\varepsilon}{3}$', ...
    '$\frac{\varepsilon}{2}$', ...
    '$\varepsilon$'})

% Legend labels
legendLabels = cell(1, n);
for i = 1:n
    legendLabels{i} = sprintf('$\\epsilon_{%d}$', i);
end
%legendLabels{n+1} = '$\epsilon / n$';

lg = legend(legendLabels, ...
   'Location', 'northeast', ...
   'Orientation', 'horizontal', ...
   'Box', 'off', ...
   'FontSize', font_leg, ...
   'Interpreter','latex');

grid on

set(gca, ...
    'Color', 'w', ...
    'YGrid', 'on', ...
    'XGrid', 'off', ...
    'GridAlpha', 0.12, ...
    'GridColor', [0 0 0], ...
    'FontSize', font_tick, ...
    'FontName', 'Times New Roman', ...
    'TickLabelInterpreter','latex', ...
    'TickDir','out', ...
    'LineWidth', 1.1, ...
    'Box', 'on');