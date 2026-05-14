%% Visualize drought parameter PDFs
clear; clc; close all;

figure('Position',[200 200 300 900])

%% 1. Baseline flow q0 ~ N(mu, sigma)
mu_q = 0.0425;
sigma_q = 0.004;

x = linspace(mu_q-4*sigma_q, mu_q+4*sigma_q, 500);
pdf_q = normpdf(x, mu_q, sigma_q);

subplot(3,1,1)
plot(x, pdf_q, 'LineWidth',2)
xlabel('q_0'); ylabel('PDF')
title('Baseline Flow: q_0 ~ N(\mu,\sigma)')
grid on

%% 2. Drought amplitude alpha ~ Beta(1.5, 8.5)
a = 1.5;
b = 8.5;

x = linspace(0,1,500);
pdf_a = betapdf(x, a, b);

subplot(3,1,2)
plot(x, pdf_a, 'LineWidth',2)
xlabel('\alpha'); ylabel('PDF')
title('\alpha ~ Beta(1.5,8.5)')
grid on

%% 3. Drought duration D ~ Gamma(4,0.2)
k = 4;
theta = 0.2;

x = linspace(0,2,500);
pdf_D = gampdf(x, k, theta);

subplot(3,1,3)
plot(x, pdf_D, 'LineWidth',2)
xlabel('D (days)'); ylabel('PDF')
title('D ~ \Gamma(4,0.2)')
grid on

%% Plot #2: Two Gaussian CDFs

% Means
mu = [0; 0];

% Variances / standard deviations
sigma1 = 1.5;   % larger uncertainty
sigma2 = 0.6;   % smaller uncertainty

% Diagonal covariance matrix
Sigma_DDU = [sigma1^2 0;
             0        sigma2^2];

% Grid
x1 = linspace(-5, 5, 200);
x2 = linspace(-2, 2, 200);
[X1, X2] = meshgrid(x1, x2);

% 2D Gaussian PDF
Z = zeros(size(X1));
for i = 1:numel(X1)
    x = [X1(i); X2(i)];
    Z(i) = 1/(2*pi*sqrt(det(Sigma_DDU))) * ...
           exp(-0.5*(x-mu)'*(Sigma_DDU\(x-mu)));
end

% Plot
figure;
set(gcf, 'Color', 'w');

contour(X1, X2, Z, 8, 'LineWidth', 1.5);
hold on;

% Mean point
plot(mu(1), mu(2), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5);

% Axes labels
% xlabel('Forecast error of reservoir 1');
% ylabel('Forecast error of reservoir 2');

% title('\Sigma_{DDU} = diag(\sigma_1^2,\sigma_2^2),  \sigma_1 > \sigma_2');

% Use the outer contour scale
r = sqrt(5.991);   % 95% confidence ellipse for 2D Gaussian

% Arrow endpoints from mean to contour edge
xEnd = mu(1) + r*sigma1;
yEnd = mu(2) + r*sigma2;

% Arrows from mu to contour boundary
% quiver(mu(1), mu(2), r*sigma1, 0, 0, ...
%     'k', 'LineWidth', 1.5, 'MaxHeadSize', 0.25);

% quiver(mu(1), mu(2), 0, r*sigma2, 0, ...
%     'k', 'LineWidth', 1.5, 'MaxHeadSize', 0.25);

% Labels outside contour
text(xEnd + 0.05, mu(2), '\sigma_1', ...
    'FontName','Times', 'FontSize',24, ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','middle');

text(mu(1), yEnd + 0.05, '\sigma_2', ...
    'FontName','Times', 'FontSize',24, ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','bottom');

% IEEE formatting
grid off;
box off;
axis off;

set(gca, 'FontName', 'Times', ...
         'FontSize', 10, ...
         'LineWidth', 1);
