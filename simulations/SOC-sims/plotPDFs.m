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