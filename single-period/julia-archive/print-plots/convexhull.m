
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot_power_surface.m
% 
% Visualizes the convex hull approximation grid for hydropower production.
% Plots 3D lines over turbine release and volume breakpoints with power.
% Ideal for verifying lambda-based convex combinations in hydro models.
%
% Author: Eliza Cohn
% Date:   July 24, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function m3s = cfs_to_m3s(cfs)

conversion_factor = 0.0283168;  % 1 cfs = 0.0283168 m³/s
m3s = cfs * conversion_factor;

end

function meters = ft_to_m(feet)

conversion_factor = 0.3048;  % 1 foot = 0.3048 meters
meters = feet * conversion_factor;

end

% Parameters
s2hr = 3600;
eta = 0.9;
g = 9.81;      % m/s^2
rho = 1000;    % kg/m^3
scale = 3.6e9;

% Forebay elevation relationship
a = 0.0026;
b = 0.4404;

% Linear constraints
min_u = s2hr*cfs_to_m3s(30600);    % min water release rate  [m3/hr] 
max_u = s2hr*cfs_to_m3s(134100);   % max water release rate [m3/hr]
min_h = ft_to_m(70);                 % min forebay elevation levels [m]
max_h = ft_to_m(77);                 % max forebay elevation levels [m]
min_v = (min_h/a)^(1/b);           % min volume [m3]
max_v = (max_h/a)^(1/b);           % max volume [m3]

% Breakpoints
N = 50;
M = 50;
U_vals = linspace(min_u, max_u, N);    % N turbine release points
V_vals = linspace(min_v, max_v, M);    % M volume breakpoints (Mm^3)

% Pre-allocate power grid
P = zeros(length(V_vals), length(U_vals)); % M x N

% Compute power values at each (n,m)
for m = 1:length(V_vals)
    for n = 1:length(U_vals)
        y = U_vals(n) * (V_vals(m)^b); 
        P(m,n) = eta * g * rho * a * y / scale;
    end
end

%% Create 3D wireframe plot
%{
figure;
hold on;

% Plot lines across N (varying U for fixed V)
for m = 1:length(V_vals)
    plot3(U_vals, V_vals(m)*ones(size(U_vals)), P(m,:), 'b', 'LineWidth', 1.5);
end

% Plot lines across M (varying V for fixed U)
for n = 1:length(U_vals)
    plot3(U_vals(n)*ones(size(V_vals)), V_vals, P(:,n), 'b', 'LineWidth', 1.5);
end

xlabel('Turbine Release U (m^3/s)');
ylabel('Reservoir Volume V (m^3)');
zlabel('Power (MW)');
title('Convex Hull Support Grid');
grid on;
view(45, 30); % angled 3D view
%}

%% Plot heatmap with color representing power
fontsize = 16;

figure;
contourf(U_grid, V_grid, P, 20);
colormap(parula)
cb = colorbar;
cb.Label.String = 'Power (MWh)';      
cb.Label.FontSize = fontsize; 
hold on;

%{
% Draw vertical grid lines (constant u)
for n = 1:length(U_vals)
    x = U_vals(n) * ones(size(V_vals));
    y = V_vals;
    line(x, y, 'Color', [0.2 0.2 0.2], 'LineStyle', '-', 'LineWidth', 0.5);
end

% Draw horizontal grid lines (constant v)
for m = 1:length(V_vals)
    x = U_vals;
    y = V_vals(m) * ones(size(U_vals));
    line(x, y, 'Color', [0.2 0.2 0.2], 'LineStyle', '-', 'LineWidth', 0.5);
end
%}

xlabel('Turbine Release u (m^3/s)', 'FontSize', fontsize);
ylabel('Reservoir Volume v (m^3)', 'FontSize', fontsize);

title('Power Production Heatmap', 'FontSize', fontsize);
set(gca, 'YDir', 'normal'); % flip y-axis so low volume at bottom