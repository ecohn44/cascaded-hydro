function plotDroughtParameterPDFs(scale)
% Plot drought parameter PDFs in physical units
%
% Inputs:
%   scale : struct with fields
%       dry_muQ   = physical mean flow (m^3/s)
%       dry_stdQ  = physical std flow (m^3/s)

    figure('Color','w','Position',[200 200 600 750]);

    % Style settings
    fs_ax  = 18;
    fs_lab = 18;
    fs_txt = 18;
    lw     = 2.5;
    col    = [0 0.2 0.6];
    red    = [0.85 0.1 0.1];

    % -----------------------------
    % 1) Baseline flow q0 ~ N(mu, sigma)
    % -----------------------------
    mu_q    = 0.0425;
    sigma_q = 0.004;

    mu_Q    = scale.dry_muQ;
    sigma_Q = scale.dry_stdQ;

    xQ = linspace(mu_Q - 4*sigma_Q, mu_Q + 4*sigma_Q, 500);
    pdf_Q = normpdf(xQ, mu_Q, sigma_Q);

    subplot(3,1,1); hold on
    plot(xQ, pdf_Q, 'LineWidth', lw, 'Color', col)

    q0_phys = scale.dry_muQ + (scale.dry_stdQ/sigma_q) * (mu_q - mu_q);  % = dry_muQ
    xline(q0_phys, '--', 'Color', red, 'LineWidth', 1.8)
    scatter(q0_phys, normpdf(q0_phys, mu_Q, sigma_Q), 35, red, 'filled')

    xlabel('$q_0$ (m$^3$/s)', 'Interpreter','latex', 'FontSize', fs_lab)
    ylabel('PDF', 'Interpreter','latex', 'FontSize', fs_lab)

    text(0.95, 0.85, '$q_0 \sim \mathcal{N}(\mu_q^{dry},\sigma_q^{dry})$', ...
        'Units','normalized', ...
        'Interpreter','latex', ...
        'FontSize', fs_txt, ...
        'HorizontalAlignment','right');

    grid on; box on
    set(gca, 'FontSize', fs_ax, 'LineWidth', 0.8)

    % -----------------------------
    % 2) Drought amplitude alpha ~ Beta(1.5, 8.5)
    % -----------------------------
    a = 1.5;
    b = 8.5;

    xa = linspace(0,1,500);
    pdf_a = betapdf(xa, a, b);

    subplot(3,1,2); hold on
    plot(xa, pdf_a, 'LineWidth', lw, 'Color', col)

    alpha_mean = a/(a+b);
    xline(alpha_mean, '--', 'Color', red, 'LineWidth', 1.8)
    scatter(alpha_mean, betapdf(alpha_mean, a, b), 35, red, 'filled')

    xlabel('$\alpha$', 'Interpreter','latex', 'FontSize', fs_lab)
    ylabel('PDF', 'Interpreter','latex', 'FontSize', fs_lab)

    text(0.95, 0.85, '$\alpha \sim \mathrm{Beta}(1.5,8.5)$', ...
        'Units','normalized', ...
        'Interpreter','latex', ...
        'FontSize', fs_txt, ...
        'HorizontalAlignment','right');

    grid on; box on
    set(gca, 'FontSize', fs_ax, 'LineWidth', 0.8)

    % -----------------------------
    % 3) Drought duration D ~ Gamma(4, 0.2)
    % -----------------------------
    k = 4;
    theta = 0.2;

    xD = linspace(0,2,500);
    pdf_D = gampdf(xD, k, theta);

    subplot(3,1,3); hold on
    plot(xD, pdf_D, 'LineWidth', lw, 'Color', col)

    D_mean = k*theta;
    xline(D_mean, '--', 'Color', red, 'LineWidth', 1.8)
    scatter(D_mean, gampdf(D_mean, k, theta), 35, red, 'filled')

    xlabel('$D$ (days)', 'Interpreter','latex', 'FontSize', fs_lab)
    ylabel('PDF', 'Interpreter','latex', 'FontSize', fs_lab)

    text(0.95, 0.85, '$D \sim \Gamma(4,0.2)$', ...
        'Units','normalized', ...
        'Interpreter','latex', ...
        'FontSize', fs_txt, ...
        'HorizontalAlignment','right');

    grid on; box on
    set(gca, 'FontSize', fs_ax, 'LineWidth', 0.8)
end