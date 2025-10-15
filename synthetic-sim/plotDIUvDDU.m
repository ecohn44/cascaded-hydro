function plotDIUvDDU(X_diu, std_diu, X_ddu, std_ddu, s, savePath, opts, Ueff_diu, Ueff_ddu)
% plotDIUvDDU  Minimal DIU vs DDU plots (Unit 2 only) with optional corridors
% Required:
%   X_diu, X_ddu   : [T x 10] solver outputs
%                    cols: 1=V1,2=p1,3=u1,4=s1,5=q1, 6=V2,7=p2,8=u2,9=s2,10=q2
%   std_diu/_ddu   : [T x 2]  forecast std; col 2 is downstream (can be [])
%   s              : sysparams struct array (uses s(2))
%   savePath       : '' to skip saving, else folder path
%   opts.titleTag  : (optional) text for titles
% Optional (for shaded SoC corridors based on post-clip Ueff):
%   Ueff_diu/_ddu  : [T x 4] effective outflow bands; cols 3=U_hi (unit 2), 4=U_lo (unit 2)
%
% Produces two figures:
%   (A) σ₂(t): DIU vs DDU (+ event marker at DDU spike, if available)
%   (B) SoC(t): DIU vs DDU (+ same event marker)
%       with two corridor types (if available):
%         - Ueff corridors (what the solver exported)   [green/red, semi-transparent]
%         - Reconstructed "clamped" corridors from storage+ramp logic (post-process)
%           (shown hatched outlines to distinguish from Ueff)

    if nargin < 6 || isempty(savePath), savePath = ''; end
    if nargin < 7 || ~isstruct(opts), opts = struct; end
    if ~isfield(opts,'titleTag'), opts.titleTag = ''; end
    if ~isempty(savePath) && ~exist(savePath,'dir'), mkdir(savePath); end

    T = size(X_ddu,1); t = (1:T)';

    % ============ Figure A: sigma2(t) ============
    haveSigma_diu = ~isempty(std_diu) && size(std_diu,2) >= 2;
    haveSigma_ddu = ~isempty(std_ddu) && size(std_ddu,2) >= 2;

    % Find DDU event time (max σ2), if available
    t_evt = NaN;
    if haveSigma_ddu
        [~, t_evt] = max(std_ddu(:,2));
    end

    fh1 = figure('Position',[100 100 900 360]); hold on; grid on;
    if haveSigma_diu, plot(t, std_diu(:,2), 'LineWidth',1.8, 'DisplayName','DIU'); end
    if haveSigma_ddu, plot(t, std_ddu(:,2), 'LineWidth',1.8, 'DisplayName','DDU'); end
    if ~isnan(t_evt)
        xline(t_evt, 'k--', '\sigma spike', 'LineWidth',1, ...
            'LabelOrientation','horizontal','LabelVerticalAlignment','bottom');
    end
    xlabel('Hour'); ylabel('\sigma_{q2} (m^3/hr)');
    if haveSigma_diu || haveSigma_ddu
        legend('Location','best');
    end
    title(sprintf('Unit 2 Inflow Uncertainty — DIU vs DDU %s', opts.titleTag));

    % ============ Figure B: SoC(t) with corridors ============
    % Extract Unit 2 series
    V2_diu = X_diu(:,6);  V2_ddu = X_ddu(:,6);
    u2_diu = X_diu(:,8);  u2_ddu = X_ddu(:,8);
    q2_diu = X_diu(:,10); q2_ddu = X_ddu(:,10);

    % SoC normalization
    Vmin = s(2).min_V; Vmax = s(2).max_V; span = Vmax - Vmin;
    soc_diu = (V2_diu - Vmin)/span;
    soc_ddu = (V2_ddu - Vmin)/span;

    fh2 = figure('Position',[100 480 1200 500]); hold on; grid on;
    yline(0,'b:','SoC=0'); yline(1,'b:','SoC=1'); % visual SoC bounds

    % ---- Corridor A: from Ueff (post-clip from your saved results) ----
    haveUeff = (nargin >= 9) && ~isempty(Ueff_diu) && ~isempty(Ueff_ddu) ...
               && size(Ueff_diu,2) >= 4 && size(Ueff_ddu,2) >= 4;

    if haveUeff
        V0 = getfield_def(s(2), 'V0', V2_diu(1)); % fallback: first state if V0 missing

        soc_fill_from_ueff = @(V,q,Uhi,Ulo) deal( ...
            ([V0; V(1:end-1)] + q - Uhi - Vmin)/span, ... % lower bound (max outflow)
            ([V0; V(1:end-1)] + q - Ulo - Vmin)/span);    % upper bound (min outflow)

        [soc_lo_diu_ueff, soc_hi_diu_ueff] = soc_fill_from_ueff(V2_diu, q2_diu, Ueff_diu(:,3), Ueff_diu(:,4));
        [soc_lo_ddu_ueff, soc_hi_ddu_ueff] = soc_fill_from_ueff(V2_ddu, q2_ddu, Ueff_ddu(:,3), Ueff_ddu(:,4));

        % Shaded corridors from Ueff (DIU green / DDU red)
        fill_patch(t, soc_lo_diu_ueff, soc_hi_diu_ueff, [0.70 1.00 0.70], 0.40, 'DIU U_{eff} corridor');
        fill_patch(t, soc_lo_ddu_ueff, soc_hi_ddu_ueff, [1.00 0.60 0.40], 0.40, 'DDU U_{eff} corridor');
    end

    % ---- Corridor B: reconstructed "clamped" bands from storage+ramp logic ----
    % This simulates your commented block post-hoc (does NOT change solved X).
    [flow_lo_diu, flow_hi_diu] = reconstructClampedFlows_U2(X_diu, s(2));
    [flow_lo_ddu, flow_hi_ddu] = reconstructClampedFlows_U2(X_ddu, s(2));

    % Build SoC corridors from reconstructed clamp bands
    V0 = getfield_def(s(2), 'V0', V2_diu(1));
    [soc_lo_diu_clamp, soc_hi_diu_clamp] = soc_fill_from_flows(V2_diu, q2_diu, flow_hi_diu, flow_lo_diu, V0, Vmin, span);
    [soc_lo_ddu_clamp, soc_hi_ddu_clamp] = soc_fill_from_flows(V2_ddu, q2_ddu, flow_hi_ddu, flow_lo_ddu, V0, Vmin, span);

    % Draw clamped corridors as **outlined** patches (to distinguish from Ueff)
    outline_patch(t, soc_lo_diu_clamp, soc_hi_diu_clamp, [0.0 0.5 0.0], 'DIU clamp corridor');
    outline_patch(t, soc_lo_ddu_clamp, soc_hi_ddu_clamp, [0.6 0.0 0.0], 'DDU clamp corridor');

    % ---- Strategy lines ----
    plot(t, soc_diu, 'Color',[0.1 0.5 0.1], 'LineWidth',2.0, 'DisplayName','DIU SoC');
    plot(t, soc_ddu, 'Color',[0.6 0.0 0.0], 'LineWidth',2.0, 'DisplayName','DDU SoC');

    % Event marker on SoC plot (if sigma available)
    if ~isnan(t_evt)
        xline(t_evt, 'k--', 'risk event', 'LineWidth',1, ...
            'LabelOrientation','horizontal','LabelVerticalAlignment','bottom');
    end

    % Legend/title
    leg = {};
    leg{end+1} = 'SoC=0'; leg{end+1} = 'SoC=1';
    if haveUeff
        leg{end+1} = 'DIU U_{eff} corridor';
        leg{end+1} = 'DDU U_{eff} corridor';
    end
    leg{end+1} = 'DIU clamp corridor';
    leg{end+1} = 'DDU clamp corridor';
    leg{end+1} = 'DIU SoC';
    leg{end+1} = 'DDU SoC';

    legend(leg, 'Location','best');
    xlabel('Hour'); ylabel('SoC = (V - V_{min})/(V_{max}-V_{min})');
    title(sprintf('Unit 2 — SoC: DIU (green) vs DDU (red) %s', opts.titleTag));

    % Save (optional)
    if ~isempty(savePath)
        saveas(fh1, fullfile(savePath, 'u2_sigma_diu_vs_ddu.png'));
        saveas(fh2, fullfile(savePath, 'u2_soc_diu_vs_ddu.png'));
    end
end

% ================= helpers =================

function [flow_lo, flow_hi] = reconstructClampedFlows_U2(X, sp2)
% Post-process reconstruction of feasible outflow band for Unit 2
% Implements your commented logic (storage-implied band intersected with ramp limits).
%
% Inputs:
%   X   : [T x 10] decision matrix
%   sp2 : sysparams for unit 2 (fields: min_V, max_V, min_ut, max_ut, RR_dn, RR_up)
%
% Outputs:
%   flow_lo, flow_hi : [T x 1] reconstructed band for u2 at each t

    T = size(X,1);
    flow_lo = zeros(T,1);
    flow_hi = zeros(T,1);

    V2  = X(:,6);
    u2  = X(:,8);
    q2  = X(:,10);

    for t = 1:T
        % Previous state/action (t-1)
        if t == 1
            Vprev2 = V2(1);
            u2_prev = sp2.min_ut; % permissive start; you can switch to u2(1) if you prefer
            u_lo = sp2.min_ut;
            u_hi = sp2.max_ut;
        else
            Vprev2 = V2(t-1);
            u2_prev = u2(t-1);
            % signed convention: RR_dn <= 0, RR_up >= 0
            u_lo = max(sp2.min_ut, u2_prev + sp2.RR_dn);
            u_hi = min(sp2.max_ut, u2_prev + sp2.RR_up);
        end

        % --- Storage-implied band (chance-constraint style) ---
        % ensures V_t >= min_V and V_t <= max_V:
        %   u_t <= V_{t-1} + q_t - min_V
        %   u_t >= V_{t-1} + q_t - max_V
        fmax = Vprev2 + q2(t) - sp2.min_V; % upper bound on u_t
        fmin = Vprev2 + q2(t) - sp2.max_V; % lower bound on u_t

        % --- Intersect with ramps ---
        fmax = min(max(fmax, u_lo), u_hi);
        fmin = min(max(fmin, u_lo), u_hi);

        % --- Keep ordering (clip collapse safe) ---
        if fmin > fmax, fmin = fmax; end

        flow_lo(t) = fmin;
        flow_hi(t) = fmax;
    end
end

function [soc_lo, soc_hi] = soc_fill_from_flows(V, q, u_hi, u_lo, V0, Vmin, span)
% Convert flow bands into SoC corridor using mass balance:
%   V_t in [ V_{t-1} + q_t - u_hi ,  V_{t-1} + q_t - u_lo ]
    Vprev = [V0; V(1:end-1)];
    soc_lo = (Vprev + q - u_hi - Vmin) ./ span; % lower SoC bound (max outflow)
    soc_hi = (Vprev + q - u_lo - Vmin) ./ span; % upper SoC bound (min outflow)
end

function fill_patch(t, ylo, yhi, rgb, alpha, name)
% Shaded corridor (filled)
    patch([t; flipud(t)], [ylo; flipud(yhi)], rgb, ...
          'EdgeColor','none', 'FaceAlpha',alpha, 'DisplayName',name);
end

function outline_patch(t, ylo, yhi, rgb, name)
% Outlined corridor (to distinguish from filled Ueff)
    plot(t, ylo, '-', 'Color', rgb, 'LineWidth', 1.25, 'DisplayName', name);
    plot(t, yhi, '-', 'Color', rgb, 'LineWidth', 1.25, 'HandleVisibility','off');
end

function val = getfield_def(structVal, fieldName, defaultVal)
% Safe getfield with default
    if isfield(structVal, fieldName)
        val = structVal.(fieldName);
    else
        val = defaultVal;
    end
end
