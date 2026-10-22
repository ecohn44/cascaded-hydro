%% ========================================================================
% RUN SETUP
% ========================================================================

function D = dispatchSetup(C, policy, cfg)
% Turn one case, policy and config into the setup a whole run needs
%
% Inputs:
%   C       - the case struct from dataload
%   policy  - 'DIU' or 'DDU'; selects whether sigma responds to decisions
%   cfg     - the run's options; every default is filled in here
%
% Outputs:
%   D       - a flat struct carrying the setup under the names the caller
%             uses, so realtimeSimulation restores its locals in one block
%
% Second of the three tiers, and the only one that reads cfg:
%
%   dataload(root,year,...)        once per sweep   -> C      the water year
%   dispatchSetup(C,policy,cfg)    once per job     -> D      this run's setup
%   realtimeDispatch(t,S,ctx)      once per hour    -> model  the hour's LP
%
% Everything here is fixed for the whole run: the horizon and unit scaling,
% the flat variable index map, the constraint rows and bounds whose
% coefficients do not depend on the real-time hour, the forecast-error model
% and its correlation, and the safety multiplier z. Nothing here looks at
% reservoir state, so a run that fails validation fails before hour 1.
%
% It is separate from dataload for a reason worth keeping: dataload is called
% ONCE per sweep and the result is sliced across a parfor of hundreds of jobs,
% while this depends on cfg and so must be rebuilt for every one of them.
% Folding the two together would reload the water year hundreds of times.
%
% See realtime.md.

%% 1. Dimensions, power model and solver switches

% n comes from the case, not a literal. It is 4 for every published run,
% so this is inert for them, but tileCase builds 8- and 12-unit cascades
% for the scalability study and a hardcoded 4 silently runs only the first
% four units -- then C.qLow (T x 8) meets V (T x 4) in the scoring below.
n = numel(C.units); H = cfg.Horizon; NV = 10; NR = 10;
% PowerModel: 'fixed' (default; head fixed at the implemented state for the
% whole look-ahead) or 'mccormick' (head linear in planned storage; the
% bilinear storage*release product is replaced by its McCormick envelope,
% exactly as in the paper, eq. (2)). The chance constraint, Sigma, SSH and
% all diagnostics are identical for both power models.
if ~isfield(cfg,'PowerModel'), cfg.PowerModel='fixed'; end
mcc=strcmpi(cfg.PowerModel,'mccormick');
if mcc, NV=11; NR=14; end
T = min([cfg.MaxHours, size(C.ref,1), size(C.zseq,1)]);
KAF = C.KAF; CPOW = C.CPOW;
if ~isfield(cfg,'SolutionMethod'), cfg.SolutionMethod='BON'; end
solutionMethod=upper(char(cfg.SolutionMethod));
if ~ismember(solutionMethod,{'BON','SSH'})
    error('SolutionMethod must be ''BON'' or ''SSH'' (received %s).',solutionMethod);
end
if ~isfield(cfg,'SSHMaxIterations'), cfg.SSHMaxIterations=20; end
if ~isfield(cfg,'SSHPhiTolerance'), cfg.SSHPhiTolerance=5e-4; end
if ~isfield(cfg,'SSHBoundaryTolerance'), cfg.SSHBoundaryTolerance=1e-5; end
if ~isfield(cfg,'SSHLookaheadHours'), cfg.SSHLookaheadHours=H; end
cfg.SSHLookaheadHours=max(1,min(cfg.SSHLookaheadHours,H));
if ~isfield(cfg,'SSHInteriorPoint'), cfg.SSHInteriorPoint='bonferroni'; end
% TrackBand: tracking deadband as a fraction of pool size. Deviations
% |V-ref| <= TrackBand*R are free; only the excess is penalized by Theta.
% TrackBand=0 (default) reproduces plain L1 tracking.
if ~isfield(cfg,'TrackBand'), cfg.TrackBand=0; end
% TrackSide: 'both' (default, symmetric L1) or 'below' (penalize only V<ref-band,
% so full pools are never pulled down toward the historical mean).
if ~isfield(cfg,'TrackSide'), cfg.TrackSide='both'; end
trackBelow=strcmpi(cfg.TrackSide,'below');
if ~isfield(cfg,'ProgressHours'), cfg.ProgressHours=168; end
if ~isfield(cfg,'Verbose'), cfg.Verbose=true; end
% ChanceSides: 'lower' (default) enforces the chance constraint only on the
% P20 lower storage bound. The upper bound is the physical maximum, and any
% water above it is spilled by physical recourse, so an upper chance
% constraint only wastes risk budget. 'both' reproduces the original
% two-sided formulation (Bonferroni splits epsilon over 2n tails).
if ~isfield(cfg,'ChanceSides'), cfg.ChanceSides='lower'; end
lowerOnly=strcmpi(cfg.ChanceSides,'lower');
runTimer=tic;
scenarioID=NaN; scenarioSeed=NaN;
if isfield(cfg,'ScenarioID'), scenarioID=cfg.ScenarioID; end
if isfield(cfg,'ScenarioSeed'), scenarioSeed=cfg.ScenarioSeed; end

%% 2. Flat variable index map
% ix(unit,hour,var) addresses the single stacked decision vector; the
% variable and row codes below name its blocks.
VOL=1; REL=2; SPL=3; POW=4; TRK=5; RMP=6; SLO=7; SHI=8; MINSLK=9; WDEF=10;
RUP=1; RDN=2; RPOS=3; RNEG=4; RLO=5; RHI=6; RTP=7; RTN=8; RPW=9; RMIN=10;
PROD=11; MC1=11; MC2=12; MC3=13; MC4=14;   % McCormick only (NV=11, NR=14)
ix = @(i,k,v) ((i-1)*H + (k-1))*NV + v;
nvar = n*H*NV;

%% 3. Forecast-error model
% Fit sigma(|du|) = intercept + slope*|du|.
sig0fit = zeros(1,n); gammaData = zeros(1,n); meanDu = zeros(1,n);
for i = 1:n
    fit = polyfit(C.sigMid(i,:), C.sigSd(i,:), 1);
    gammaData(i) = max(fit(1),0);
    sig0fit(i) = max(fit(2),0.05*mean(C.sigSd(i,:)));
    meanDu(i) = mean(C.sigMid(i,:));
end
sigDIU = sig0fit + gammaData.*meanDu;

% ---------------------------------------------------------------------------
% ARCH-X conditional variance (cfg.SigmaModel = 'binned' default | 'archx').
% 'binned' reproduces every earlier result bit-for-bit.
%   sigma^2_{i,m} = omega_i(kappa) + alpha_i e^2_{i,m-1}
%                   + kappa*( gu_i*ut_{j} + gr_i*|dut_{j}| )
% with ut, |dut| the upstream release level and ramp, min-max normalised on
% the estimation range. Inside the look-ahead e^2 is unknown for m>=2, so the
% standard multi-step recursion replaces it with sigma^2_{m-1}.
if ~isfield(cfg,'SigmaModel'), cfg.SigmaModel='binned'; end
useArchX = strcmpi(cfg.SigmaModel,'archx');
AX = struct();
if useArchX
    % This file sits at the repository root, so one fileparts reaches it.
    axFile = fullfile(fileparts(mfilename('fullpath')), ...
        'data','archx_params.csv');
    Tax = readtable(axFile);
    Tax.unit = string(Tax.unit);
    AX.omega=zeros(1,n); AX.alpha=zeros(1,n); AX.gu=zeros(1,n); AX.gr=zeros(1,n);
    AX.umin=zeros(1,n);  AX.urng=ones(1,n);   AX.rmin=zeros(1,n); AX.rrng=ones(1,n);
    AX.sbar=zeros(1,n);  AX.ubar=zeros(1,n);  AX.dbar=zeros(1,n);
    for i=1:n
        % archx_params.csv carries only the four real projects, while
        % tileCase names replicas MCN_2, JDA_2 and so on. Strip a trailing
        % _<digits> so a replica inherits its base unit's fitted parameters.
        % Inert for the four real units, which carry no suffix. Without this
        % the lookup returns an empty row and the assignment below fails with
        % "left and right sides have a different number of elements".
        baseUnit = regexprep(C.units{i}, '_\d+$', '');
        r = Tax(Tax.unit==string(baseUnit),:);
        assert(height(r) == 1, 'dispatchSetup:noArchXRow', ...
            ['data/archx_params.csv has %d rows for "%s" (unit "%s"); ' ...
             'expected exactly 1.'], height(r), baseUnit, C.units{i});
        AX.omega(i)=r.omega;  AX.alpha(i)=r.alpha;
        AX.gu(i)=r.gamma_release; AX.gr(i)=r.gamma_ramp;
        AX.umin(i)=r.release_min; AX.urng(i)=max(r.release_range,eps);
        AX.rmin(i)=r.ramp_min;    AX.rrng(i)=max(r.ramp_range,eps);
        AX.sbar(i)=r.sigma_DIU;
        % mean normalised operating point, for the kappa-invariance intercept
        if i==1, up=C.upReleaseExo(1:T); else, up=C.outObs(1:T,i-1); end
        AX.ubar(i)=mean(max((up-AX.umin(i))/AX.urng(i),0));
        AX.dbar(i)=mean(max((abs([0;diff(up)])-AX.rmin(i))/AX.rrng(i),0));
    end
    % DIU carries the UNCONDITIONAL residual s.d.
    sigDIU = AX.sbar;
end

%% 4. Belief knobs: Kappa, PolicyKappa, and the sigma evaluator
% Frozen copy of the fitted level. The realized world is ALWAYS built from
% these data values, never from the policy-side quantities below, so the
% belief knobs (Kappa, PolicyKappa) cannot move the testing environment.
% Mirrors sig_diu_data in Ning's dispatch_lookahead.py.
sigDIUdata = sigDIU;

% Kappa changes only the ramp-sensitive component and preserves the
% fitted conditional standard deviation at the mean ramp when possible.
gammaModel = cfg.Kappa*gammaData;
if strcmpi(policy,'DDU')
    gamma = gammaModel;
    sig0 = max(sigDIU-gamma.*meanDu,0.2*sigDIU);
else
    gamma = zeros(1,n);
    sig0 = sigDIU;
end

% PolicyKappa: scales the WHOLE sigma the optimizer assumes (level and ramp
% term together), the belief-side twin of WorldKappa. PolicyKappa = WorldKappa
% is an operator who knows the forecast is k times more uncertain than fitted;
% PolicyKappa = 1 is one who does not. Matches k_policy in Ning's code.
if ~isfield(cfg,'PolicyKappa'), cfg.PolicyKappa=1; end
if cfg.PolicyKappa~=1
    sigDIU = cfg.PolicyKappa*sigDIU;
    sig0   = cfg.PolicyKappa*sig0;
    gamma  = cfg.PolicyKappa*gamma;
end

% Conditional s.d. used everywhere below.
%   sigEval(i,uLev,duRamp,s2prev) -> [sigma, sigma^2]
% 'binned' ignores uLev and s2prev and returns the affine ramp model, so the
% call sites are identical in both modes.
sigCeil = 3.0*sigDIU;                     % ceiling, mirroring the 0.2*sigDIU floor
if useArchX
    kap = cfg.Kappa; if strcmpi(policy,'DIU'), kap = 0; end
    omgK = AX.omega + (1-kap)*(AX.gu.*AX.ubar + AX.gr.*AX.dbar);
    omgK = max(omgK, 1e-6);
    sigEval = @(i,uLev,duR,s2p) forecastError(i,uLev,duR,s2p,AX,omgK,kap, ...
        cfg.PolicyKappa,0.2*sigDIU(i),sigCeil(i));
else
    sigEval = @(i,uLev,duR,s2p) deal(sig0(i)+gamma(i)*duR, ...
        (sig0(i)+gamma(i)*duR)^2);
end

%% 5. Safety multiplier and the inflow-error correlation
% Bonferroni quantile for upper/lower limits at all four projects.
if lowerOnly
    nTails=n;      % one lower tail per project
else
    nTails=2*n;    % lower and upper tails per project
end
z = -sqrt(2)*erfcinv(2*(1-cfg.Epsilon/nTails));
% ChanceMode: 'on' (default) enforces the chance constraint; 'off' is the
% deterministic MPC benchmark (no safety margin, z=0). Epsilon is still used
% to score realised reliability, so metrics stay comparable.
if ~isfield(cfg,'ChanceMode'), cfg.ChanceMode='on'; end
% WorldKappa: severity of the REALIZED world (default 1 = the fitted DDU model).
% The realized error is e = WorldKappa * sigmaWorld(d) * z, i.e. z ~ K*N(0,R) in
% the flowchart notation. It is independent of cfg.Kappa, which sets only what
% the optimizer believes. WorldKappa>1 is a misspecification / stress test in
% which BOTH DIU and DDU are under-protective.
if ~isfield(cfg,'WorldKappa'), cfg.WorldKappa=1; end
if strcmpi(cfg.ChanceMode,'off'), z=0; end
if isfield(C,'Rcorr') && isequal(size(C.Rcorr),[n n])
    Rcorr=nearest_correlation_local(C.Rcorr);
else
    Rcorr=eye(n);
end
%% 6. Benchmark margin rules
% Selects how the predicted-storage row is protected. Every variant shares
% horizon, sampling period, power model, recourse, world and shock path and
% differs ONLY in this scalar, which is what makes the benchmark table one
% code path rather than four controllers.
%
%   'chance' (default) : z = Phi^-1(1 - eps/nTails), the paper's controller
%   'none'             : z = 0, deterministic MPC / certainty equivalence
%   'robust'           : z = cfg.Gamma, box or ellipsoidal robust MPC
%   'scenario'         : z = order statistic over cfg.NScen sampled paths,
%                        discarding cfg.NDiscard of them (Calafiore-Campi,
%                        Campi-Garatti). See scenarioMargin.m.
%
% PLACED HERE, not beside the ChanceMode line, because 'scenario' needs
% Rcorr, which is built just above. cfg.Epsilon is deliberately untouched in
% every case so realised reliability is still scored on the same axis.
%
% 'robust' and 'scenario' DECOUPLE the multiplier from epsilon. That is the
% entire point of both benchmarks: neither carries a calibrated probability,
% so whatever reliability they deliver is whatever the hydrology hands back.
% Do NOT derive Gamma or the sample size from epsilon to make them "fair" --
% that collapses them onto Bonferroni and destroys the comparison.
%
% Ordering note: ChanceMode='off' has already zeroed z above and
% MarginRule='chance' leaves it zeroed, so pre-existing configs are
% unaffected.
if ~isfield(cfg,'MarginRule'), cfg.MarginRule='chance'; end
marginRule=lower(char(cfg.MarginRule));
switch marginRule
    case 'chance'
        % leave z exactly as computed above
    case 'none'
        z=0;
    case 'robust'
        assert(isfield(cfg,'Gamma') && isscalar(cfg.Gamma) && ...
               isfinite(cfg.Gamma) && cfg.Gamma>=0, ...
            ['MarginRule ''robust'' needs a finite scalar cfg.Gamma >= 0. ' ...
             'Gamma is the size of the uncertainty set, NOT a probability.']);
        z=cfg.Gamma;
    case 'scenario'
        assert(isfield(cfg,'NScen') && isscalar(cfg.NScen) && cfg.NScen>=1, ...
            'MarginRule ''scenario'' needs cfg.NScen >= 1 (samples)');
        nDisc=0;
        if isfield(cfg,'NDiscard') && ~isempty(cfg.NDiscard), nDisc=cfg.NDiscard; end
        scenSeed=20260928;
        if isfield(cfg,'ScenSeed') && ~isempty(cfg.ScenSeed), scenSeed=cfg.ScenSeed; end
        % Deliberately NOT seeded from ScenarioID: the multiplier is a
        % property of the sampling design, not of the inflow scenario, so
        % every scenario in an arm must use the identical z or the arm is
        % not internally paired.
        % NOTE the horizon is NOT passed. the chance constraint is
        % per look-ahead step (nTails = n, one lower tail per project), so
        % the scenario statistic is the max over UNITS only. Maxing over
        % units and steps would enforce a joint-over-the-horizon constraint
        % nobody asked for and inflate z from 2.23 to 2.94 at eps = 0.05.
        [z,scenInfo]=scenarioMargin(n,Rcorr,cfg.NScen,nDisc,scenSeed); %#ok<ASGLU>
    otherwise
        error('realtimeSimulation:badMarginRule', ...
            'unknown cfg.MarginRule ''%s'' (chance | none | robust | scenario)', ...
            marginRule);
end

% Built once per run: exact closed form when Rcorr is the identity,
% deterministic fixed-lattice evaluator otherwise.
%% 7. Joint-probability evaluator
probCache=sshProbabilityCache(Rcorr);

%% 8. Rows and bounds that do not change with the hour
Aeq = sparse(n*H,nvar); Aub0 = sparse(NR*n*H,nvar);
lb = -inf(nvar,1); ub = inf(nvar,1);
varNames = cell(nvar,1);
rowNames = cell(n*H+NR*n*H,1);
varLabels = {'V','release','spill','power','tracking','ramp', ...
    'lowSlack','highSlack','minReleaseSlack','forecastWaterDeficit'};
rowLabels = {'rampUp','rampDown','rampAbsPos','rampAbsNeg', ...
    'storageLow','storageHigh','trackPos','trackNeg','powerLimit', ...
    'minimumRelease'};
if mcc
    varLabels{end+1}='product';
    rowLabels=[rowLabels {'mcLower1','mcLower2','mcUpper1','mcUpper2'}];
end
for i = 1:n
    p = C.sys(i);
    for k = 1:H
        r=(i-1)*H+k; b=(r-1)*NR;
        rowNames{r}=sprintf('balance_%s_h%02d',C.units{i},k);
        for v=1:NV
            varNames{ix(i,k,v)}=sprintf('%s_%s_h%02d', ...
                varLabels{v},C.units{i},k);
        end
        for rr=1:NR
            rowNames{n*H+b+rr}=sprintf('%s_%s_h%02d', ...
                rowLabels{rr},C.units{i},k);
        end
        Aeq(r,ix(i,k,VOL))=1;
        Aeq(r,ix(i,k,REL))=1;
        Aeq(r,ix(i,k,SPL))=1;
        % Emergency forecast-water correction enters as additional supply.
        Aeq(r,ix(i,k,WDEF))=-1;
        if k>1, Aeq(r,ix(i,k-1,VOL))=-1; end
        j=C.upstream(i); src=(k-1)-C.lag(i);
        if j>0 && src>=0
            Aeq(r,ix(j,src+1,REL))=-1;
            Aeq(r,ix(j,src+1,SPL))=-1;
        end

        Aub0(b+RPW,ix(i,k,POW))=1;
        Aub0(b+RUP,ix(i,k,REL))=1;
        Aub0(b+RDN,ix(i,k,REL))=-1;
        Aub0(b+RPOS,ix(i,k,REL))=1;
        Aub0(b+RPOS,ix(i,k,RMP))=-1;
        Aub0(b+RNEG,ix(i,k,REL))=-1;
        Aub0(b+RNEG,ix(i,k,RMP))=-1;
        if k>1
            Aub0(b+RUP,ix(i,k-1,REL))=-1;
            Aub0(b+RDN,ix(i,k-1,REL))=1;
            Aub0(b+RPOS,ix(i,k-1,REL))=-1;
            Aub0(b+RNEG,ix(i,k-1,REL))=1;
        end
        Aub0(b+RLO,ix(i,k,VOL))=-1;
        Aub0(b+RLO,ix(i,k,SLO))=-1;
        Aub0(b+RHI,ix(i,k,VOL))=1;
        Aub0(b+RHI,ix(i,k,SHI))=-1;
        Aub0(b+RTP,ix(i,k,VOL))=1;
        Aub0(b+RTP,ix(i,k,TRK))=-1;
        Aub0(b+RTN,ix(i,k,VOL))=-1;
        Aub0(b+RTN,ix(i,k,TRK))=-1;
        % release + emergency slack >= minimum release
        Aub0(b+RMIN,ix(i,k,REL))=-1;
        Aub0(b+RMIN,ix(i,k,MINSLK))=-1;
        if mcc
            % McCormick envelope for w = V*release, V in [0,R], release in
            % [yL,yU]. Storage coefficients (yL, yU) and right-hand sides
            % are set every rolling hour from the ramp-reachable bounds.
            Aub0(b+MC1,ix(i,k,PROD))=-1;
            Aub0(b+MC2,ix(i,k,REL))= p.R; Aub0(b+MC2,ix(i,k,PROD))=-1;
            Aub0(b+MC3,ix(i,k,REL))=-p.R; Aub0(b+MC3,ix(i,k,PROD))= 1;
            Aub0(b+MC4,ix(i,k,PROD))= 1;
        end

        lb(ix(i,k,VOL))=0; ub(ix(i,k,VOL))=p.R;
        lb(ix(i,k,REL))=0; ub(ix(i,k,REL))=p.umax*KAF;
        lb(ix(i,k,SPL))=0;
        lb(ix(i,k,POW))=0; ub(ix(i,k,POW))=p.pmax;
        lb(ix(i,k,TRK))=0; lb(ix(i,k,RMP))=0;
        lb(ix(i,k,SLO))=0; lb(ix(i,k,SHI))=0;
        % Zero upper bound keeps the original minimum release hard.
        % It is unlocked only after the original LP is proven infeasible.
        lb(ix(i,k,MINSLK))=0; ub(ix(i,k,MINSLK))=0;
        lb(ix(i,k,WDEF))=0; ub(ix(i,k,WDEF))=0;
        if mcc
            lb(ix(i,k,PROD))=0; ub(ix(i,k,PROD))=p.R*p.umax*KAF;
        end
    end
end

%% 9. Pack
% Returned under the names the caller already uses, so realtimeSimulation
% restores its locals in one block and the simulation below it is
% untouched by this split.

% dimensions and unit scaling
D.n=n; D.H=H; D.T=T; D.NR=NR;
D.nvar=nvar; D.KAF=KAF; D.CPOW=CPOW; D.mcc=mcc;

% flat index map: ix(unit,hour,var), the variable codes, and the row codes
D.ix=ix; D.VOL=VOL; D.REL=REL; D.SPL=SPL;
D.POW=POW; D.TRK=TRK; D.RMP=RMP; D.SLO=SLO;
D.SHI=SHI; D.MINSLK=MINSLK; D.WDEF=WDEF; D.PROD=PROD;
D.RUP=RUP; D.RDN=RDN; D.RPOS=RPOS; D.RNEG=RNEG;
D.RLO=RLO; D.RHI=RHI; D.RTP=RTP; D.RTN=RTN;
D.RPW=RPW; D.RMIN=RMIN; D.MC1=MC1; D.MC2=MC2;
D.MC3=MC3; D.MC4=MC4;

% rows and bounds whose coefficients never change with the hour
D.Aeq=Aeq; D.Aub0=Aub0; D.lb=lb; D.ub=ub;
D.varNames=varNames; D.rowNames=rowNames;

% forecast-error model, correlation and the safety multiplier
D.sigEval=sigEval; D.sigDIU=sigDIU; D.sigDIUdata=sigDIUdata; D.gammaData=gammaData;
D.meanDu=meanDu; D.useArchX=useArchX; D.AX=AX; D.Rcorr=Rcorr;
D.probCache=probCache; D.z=z;

% resolved switches; cfg comes back because every default is filled in above
D.cfg=cfg; D.solutionMethod=solutionMethod; D.marginRule=marginRule; D.lowerOnly=lowerOnly;
D.trackBelow=trackBelow;

% run bookkeeping
D.runTimer=runTimer; D.scenarioID=scenarioID; D.scenarioSeed=scenarioSeed;

end


function R=nearest_correlation_local(R)
R=(R+R')/2; R(~isfinite(R))=0;
R(1:size(R,1)+1:end)=1;
[Q,D]=eig(R); D=diag(max(diag(D),1e-10)); R=Q*D*Q';
d=sqrt(max(diag(R),1e-12)); R=R./(d*d');
R=(R+R')/2; R(1:size(R,1)+1:end)=1;
end
