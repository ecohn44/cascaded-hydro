%% Plot Driver: hydraulic head + PWL approximation (N intervals)
% ========================================================================
% Author: Eliza Cohn
% Description: Visualize hydraulic head fit phi(V)=a V^b and its piecewise-
%              linear approximation with N intervals.
% ========================================================================

close all; clc;

% Parameters
a = 5;
b = 0.45;

N = 10;                      % number of PWL intervals
Vmin = 0.0;
Vmax = 1.0;

% Dense curve for true head function
V = linspace(Vmin, Vmax, 400);
phi = a .* V.^b;

% PWL breakpoints
Vk = linspace(Vmin, Vmax, N+1);
phik = a .* Vk.^b;

% Evaluate PWL approximation on dense grid
phi_pwl = interp1(Vk, phik, V, 'linear');

% Plot
figure('Color','w','Position',[200 200 780 360]); hold on;

% Shaded PWL feasible region
x_fill = [V, fliplr(V)];
y_fill = [phi_pwl, zeros(size(phi_pwl))];
hFill = fill(x_fill, y_fill, [0.80 0.95 1.00], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.85);
set(hFill,'HandleVisibility','off');

% True nonlinear curve
plot(V, phi, '-', ...
    'LineWidth', 3, ...
    'Color', 'k', ...
    'DisplayName', 'True Head');

% PWL approximation
plot(V, phi_pwl, '--', ...
    'LineWidth', 2.6, ...
    'Color', [0.10 0.45 1.00], ...
    'DisplayName', 'PWL Approx.');

% Breakpoints (red X)
plot(Vk, phik, 'x', ...
    'MarkerSize', 12, ...
    'LineWidth', 2.4, ...
    'Color', [0.85 0.10 0.10], ...
    'DisplayName', 'Intervals');

% Formatting
xlim([Vmin Vmax]);
ylim([0 1.08*max(phi)]);

idx = [1 2 3 (length(Vk)-1), length(Vk)];   % V1, V2, V3, ..., VN+1

% Set tick locations
set(gca,'XTick', Vk(idx));

% Symbolic LaTeX labels
xticklabels({ ...
    '$V_1$', ...
    '$V_2$', ...
    '$V_3$', ...
    '$V_N$', ...
    '$V_{N+1}$'
});


set(gca,'FontSize',20,'LineWidth',1.2,'TickLabelInterpreter','latex');
set(gca,'YTick',[]);

% Labels
xlabel('Volume (m^3)', 'FontSize',16);
ylabel('Hydraulic Head (m)', 'FontSize',16);

grid off;

box on
set(gca,'XAxisLocation','bottom','YAxisLocation','left');

yl = ylim;
text(0.5*Vmax, yl(1) - 0.1*(yl(2)-yl(1)), '$\cdots$', ...
    'Interpreter','latex', 'FontSize',18, ...
    'HorizontalAlignment','center', 'VerticalAlignment','bottom');

% Legend
lgd = legend( ...
    'Orientation','horizontal', ...
    'Location','northoutside', ...
    'Box','off');

lgd.FontSize = 14;

% Export
exportgraphics(gcf,'figures/head_function_pwl.png','Resolution',300);
