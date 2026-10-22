%% ========================================================================
% ROLLING-HORIZON SIMULATION
% ========================================================================

function out = realtimeSimulation(C, policy, cfg)
% Rolling-horizon DIU/DDU dispatch over one water year
%
% Inputs:
%   C       - the case struct from dataload
%   policy  - 'DIU' or 'DDU'
%   cfg     - the run's options; defaults are filled in by dispatchSetup
%
% Outputs:
%   out     - .summary, one struct of scalar metrics per run, and
%             .trajectory, the hourly series behind them
%
% Owns the simulation: the hour loop, reservoir state, the solve and the SSH
% handoff, the emergency recourse ladder, and the summary. Everything fixed
% for the whole run is built by dispatchSetup.m; the per-hour LP is built by
% realtimeDispatch.m, which is pure and holds no state.
%
%   realtimeSimulation -> dispatchSetup      (once, the run's setup)
%                      -> realtimeDispatch   (build the hour's model)
%                      -> applySSHNative     (solve it under a joint CC)
%
% See realtime.md.
%
% cfg.Epsilon is the joint chance-constraint risk tolerance allocated
% equally across the two storage limits and four projects at each step.
% cfg.Kappa scales the fitted ramp-dependent term in the DDU chance
% constraints. It does not change the realized inflow path. Kappa=0
% removes decision dependence and Kappa=1 uses the fitted relationship.
%
% Covariance within a rolling step (paper Sec. IV, eq. 6c and Theorem 1).
% Sigma_t is FIXED within each rolling solve and is identical for BON and
% SSH. Past upstream ramps enter at their implemented values. Upstream
% ramps that fall inside the look-ahead are predicted from this policy's
% previous rolling plan (shifted one hour), which is data at time t, not a
% decision variable. Decision dependence propagates from step to step
% through the implemented releases. This keeps the joint chance set
% convex (log-concave Gaussian CDF, fixed Sigma), so the SSH cuts are
% valid and the Bonferroni set is contained in the SSH set from the same
% state.
%
% Step feasibility (reported, never used to change decisions):
%   BON step-feasible : the Bonferroni plan needs no storage slack and no
%                       dropped margin anywhere in the look-ahead.
%   SSH step-feasible : SSH reaches 1-eps at every look-ahead hour.
% In SSH runs a same-state probe also solves BON from the identical state
% and Sigma, so containment (BON feasible => SSH feasible) and objective
% dominance (SSH objective <= BON objective) are checked every hour.

%% 1. Setup, restored into locals
% Everything that is fixed for the whole run is built in dispatchSetup.m.
D = dispatchSetup(C, policy, cfg);
% Restore the setup's locals so everything below is unchanged by the split.
n=D.n; H=D.H; T=D.T; NR=D.NR;
nvar=D.nvar; KAF=D.KAF; CPOW=D.CPOW; mcc=D.mcc;
ix=D.ix; VOL=D.VOL; REL=D.REL; SPL=D.SPL;
POW=D.POW; TRK=D.TRK; RMP=D.RMP; SLO=D.SLO;
SHI=D.SHI; MINSLK=D.MINSLK; WDEF=D.WDEF; PROD=D.PROD;
RUP=D.RUP; RDN=D.RDN; RPOS=D.RPOS; RNEG=D.RNEG;
RLO=D.RLO; RHI=D.RHI; RTP=D.RTP; RTN=D.RTN;
RPW=D.RPW; RMIN=D.RMIN; MC1=D.MC1; MC2=D.MC2;
MC3=D.MC3; MC4=D.MC4;
Aeq=D.Aeq; Aub0=D.Aub0; lb=D.lb; ub=D.ub;
varNames=D.varNames; rowNames=D.rowNames;
sigEval=D.sigEval; sigDIU=D.sigDIU; sigDIUdata=D.sigDIUdata; gammaData=D.gammaData;
meanDu=D.meanDu; useArchX=D.useArchX; AX=D.AX; Rcorr=D.Rcorr;
probCache=D.probCache; z=D.z;
cfg=D.cfg; solutionMethod=D.solutionMethod; marginRule=D.marginRule; lowerOnly=D.lowerOnly;
trackBelow=D.trackBelow;
runTimer=D.runTimer; scenarioID=D.scenarioID; scenarioSeed=D.scenarioSeed;
clear D

%% 2. Simulation state and diagnostic series
V=zeros(T,n); U=zeros(T,n); Spill=zeros(T,n); Power=zeros(T,n);
Head=zeros(T,n); Qhist=nan(T,n); LocalReal=zeros(T,n); McGap=zeros(T,n);
MinReleaseShortfall=zeros(T,n); EmergencyUsed=false(T,1);
EmergencyForecastWater=zeros(T,1); EmergencyForecastRelease=zeros(T,1);
SolveSeconds=zeros(T,1);
JointReliability=nan(T,1); RiskAllocation=nan(T,n); SSHStrain=nan(T,n);
SSHIterations=zeros(T,1); SSHCuts=zeros(T,1); SSHUsedSlater=false(T,1);
SSHLPSolves=zeros(T,1); SSHTargetUnattainable=false(T,1);
SSHTargetUsed=(1-cfg.Epsilon)*ones(T,1); unattainableCount=0;
SSHMode=repmat({''},T,1);
PlannedLowSlack=zeros(T,n); PlannedHighSlack=zeros(T,n);
Vprev=min(max(C.obs(1,:),0),C.R);
uprev=zeros(1,n);
for i=1:n
    p=C.sys(i);
    uprev(i)=min(max(C.outObs(1,i),p.umin),p.umax);
end
Qprev=C.outObs(1,:); localPrev=C.local(1,:); cumDev=zeros(1,n);
% Previous rolling plan of outflow (release+spill, flow units, the same
% quantity as the implemented du used for known ramps), used only to predict
% upstream ramps inside the look-ahead. No plan exists before hour 1.
PrevPlanRel=repmat(Qprev(:),1,H);
BonStepFeasible=true(T,1); SSHStepFeasible=true(T,1); CrossingAny=false(T,1);
% First-hour (implemented-hour) versions of the same feasibility flags.
CrossingFirst=false(T,1); BonFirstFeasible=true(T,1); SSHFirstFeasible=true(T,1);
ProbeBonFeasible=true(T,1); ProbeBonObj=nan(T,1); ProbeSshObj=nan(T,1);
SSHGuardUsed=false(T,1);
curtailMWh=0; slackHours=0; projectedInflowKaf=0; projectedInflowHours=0;
emergencySolveHours=0; emergencyCurtailMWh=0;
lastResid2=zeros(1,n);   % last observed squared forecast error, ARCH seed
%% 3. Price the planner optimises against
% PriceMode: 'flat' (default; annual mean price, energy maximisation),
% 'hourly' (the hourly price series C.lmp), or 'bpa' (see below).
if ~isfield(cfg,'PriceMode'), cfg.PriceMode='flat'; end
% 'bpa' uses the published BPA PF-24 Load Shaping rate (Schedule PF-24,
% section 2.1.3.1, FY2024-2025, FERC-approved), expanded to an hourly series
% by the GRSP HLH/LLH calendar: hours ending 0700-2200 Pacific, Monday
% through Saturday, excluding NERC holidays. Rates are used as published --
% no rescaling to match any other price profile. Built by
% paper_figures/make_bpa_price.py into data/bpa_pf24_lmp_wy2024.csv.
lmpBpa=[];
if strcmpi(cfg.PriceMode,'bpa')
    bpaFile=fullfile(fileparts(mfilename('fullpath')), ...
        'data','bpa_pf24_lmp_wy2024.csv');
    if ~isfile(bpaFile)
        error('realtimeSimulation:bpaPriceMissing', ...
            'PriceMode ''bpa'' needs %s; run paper_figures/make_bpa_price.py',bpaFile);
    end
    Bt=readtable(bpaFile);
    if ismember('price_usd_per_mwh',Bt.Properties.VariableNames)
        lmpBpa=Bt.price_usd_per_mwh;
    else
        lmpBpa=Bt{:,end};
    end
    if numel(lmpBpa)<T
        error('realtimeSimulation:bpaPriceShort','%s holds %d hours, need %d', ...
            bpaFile,numel(lmpBpa),T);
    end
end
switch lower(char(cfg.PriceMode))
    case 'hourly', lmp=C.lmp(1:T);
    case 'bpa',    lmp=lmpBpa(1:T);
    otherwise,     lmp=mean(C.lmp(1:T))*ones(T,1);
end

%% 4. Settlement price, for scoring only
% Real-time price discrepancy.
%
% The planner optimises against the price it KNOWS -- the flat annual mean or
% the scheduled hourly block curve, held in `lmp` above. Settlement happens at
% a real-time price that departs from that schedule. `lmpReal` is that realised
% price; it enters the reported revenue and NOTHING else, so the dispatch,
% storage, spill and violation metrics are bit-for-bit identical to a run with
% PriceNoiseSigma = 0. Turning this on can never change a reliability result.
%
%   lmpReal_t = C.lmp_t * exp(eta_t - sigma^2/2),   eta AR(1), sd sigma
%
% The -sigma^2/2 term makes E[exp(.)] = 1, so the noise reprices hour to hour
% without inflating or deflating annual revenue on average. The RNG is seeded
% from the scenario seed, so DIU and DDU in the same scenario are settled
% against the SAME realised price path and the comparison stays paired.
%
%   PriceNoiseSigma : log-scale sd of the real-time departure (0 = off)
%   PriceNoiseRho   : hour-to-hour persistence of that departure, in [0,1)
%   PriceNoiseSeed  : override the scenario-derived seed (optional)
if ~isfield(cfg,'PriceNoiseSigma'), cfg.PriceNoiseSigma=0; end
if ~isfield(cfg,'PriceNoiseRho'),   cfg.PriceNoiseRho=0.80; end
if ~isfield(cfg,'PriceNoiseSeed'),  cfg.PriceNoiseSeed=[]; end

% PriceFlatten compresses the price toward its own mean:
%   lambda' = lbar + (1 - f)*(lambda - lbar)
% f = 0 leaves the stylized TOU profile untouched, f = 1 is a flat price.
% This damps the peak/off-peak spread without changing the annual mean, so
% revenue levels stay comparable across settings. It exists because the
% unmodified 6-hour $65 block lets a run-of-river project cycle its whole
% usable storage daily, which real ramp and downstream-flow restrictions
% prevent.
if ~isfield(cfg,'PriceFlatten'), cfg.PriceFlatten=0; end
if cfg.PriceFlatten>0
    f=min(max(cfg.PriceFlatten,0),1);
    lbar=mean(C.lmp(1:T));
    lmp=lbar+(1-f)*(lmp-lbar);
end

lmpSched = C.lmp(1:T);           % the schedule the planner is paid against today
if strcmpi(cfg.PriceMode,'bpa'), lmpSched = lmpBpa(1:T); end
if cfg.PriceFlatten>0
    lmpSched=mean(lmpSched)+(1-min(max(cfg.PriceFlatten,0),1))*(lmpSched-mean(lmpSched));
end
lmpReal  = lmpSched;             % realised settlement price
if cfg.PriceNoiseSigma > 0
    rho = min(max(cfg.PriceNoiseRho,0),0.999);
    if isempty(cfg.PriceNoiseSeed)
        if ~isscalar(scenarioSeed) || ~isfinite(scenarioSeed)
            baseSeed = 20260924;
        else
            baseSeed = scenarioSeed;
        end
    else
        baseSeed = cfg.PriceNoiseSeed;
    end
    % Same stream for DIU and DDU: the seed depends on the scenario only.
    rs  = RandStream('mt19937ar','Seed',mod(baseSeed + 77771,2^32-1));
    eta = zeros(T,1);
    innov = sqrt(1-rho^2) * cfg.PriceNoiseSigma * randn(rs,T,1);
    eta(1) = cfg.PriceNoiseSigma * randn(rs,1);
    for tt = 2:T
        eta(tt) = rho*eta(tt-1) + innov(tt);
    end
    lmpReal = lmpSched .* exp(eta - 0.5*cfg.PriceNoiseSigma^2);
    lmpReal = max(lmpReal, 0);   % no negative settlement in this study
end

%% 5. Solver parameters
gpar.OutputFlag=0; gpar.Threads=1; gpar.Method=1;
% Hang-breaker only (2026-09-26). Measured solve times on this model are
% mean 0.032 s / p95 0.045 s (BON) and 0.180 s / 0.276 s (SSH) per hour,
% so 120 s is >2600x the BON p95 and cannot bind on any healthy solve --
% it changes no result. Its only job is to convert a non-returning solve
% into a non-OPTIMAL status, which raises the error at the fixed-head
% status check and is caught by sweepRunOne, so one job is lost and
% logged instead of the whole parfor freezing with no diagnostic.
gpar.TimeLimit=120;
if cfg.Verbose
    fprintf('%s solution method: %s.\n',upper(policy),solutionMethod);
end
Lh=cfg.SSHLookaheadHours;
nnzRows=NR*n*H;
% Triplet buffers: sparse element-by-element assignment inside the hourly
% loop is O(nnz) per insertion in MATLAB. Collect (row,col,value) instead.
maxRamp=n*H*(H+1)/2;
%% 6. Context handed to realtimeDispatch
% Built once here rather than per hour, because none of it changes with t.
ctx = struct();
ctx.Aeq = Aeq;
ctx.Aub0 = Aub0;
ctx.H = H;
ctx.KAF = KAF;
ctx.RHI = RHI;
ctx.RLO = RLO;
ctx.T = T;
ctx.cfg = cfg;
ctx.lb = lb;
ctx.lowerOnly = lowerOnly;
ctx.mcc = mcc;
ctx.n = n;
ctx.nnzRows = nnzRows;
ctx.nvar = nvar;
ctx.sigDIU = sigDIU;
ctx.sigEval = sigEval;
ctx.ub = ub;
ctx.z = z;
ctx.C = C;
ctx.CPOW = CPOW;
ctx.MC1 = MC1;
ctx.MC2 = MC2;
ctx.MC3 = MC3;
ctx.MC4 = MC4;
ctx.NR = NR;
ctx.POW = POW;
ctx.PROD = PROD;
ctx.RDN = RDN;
ctx.REL = REL;
ctx.RMIN = RMIN;
ctx.RNEG = RNEG;
ctx.RPOS = RPOS;
ctx.RPW = RPW;
ctx.RTN = RTN;
ctx.RTP = RTP;
ctx.RUP = RUP;
ctx.SHI = SHI;
ctx.SLO = SLO;
ctx.SPL = SPL;
ctx.TRK = TRK;
ctx.VOL = VOL;
ctx.ix = ix;
ctx.lmp = lmp;
ctx.maxRamp = maxRamp;
ctx.policy = policy;
ctx.trackBelow = trackBelow;
% ---------------------------------------------------------------------------

%% 7. Rolling-horizon loop
% One iteration is one real-time hour: build, solve, apply recourse,
% implement the first step, advance the state.
for t=1:T
    % t0 is the zero-based hour, used by the recourse and record phases
    % below. realtimeDispatch derives it independently for its own use.
    t0 = t-1;

    % ---- per-hour LP build, in realtimeDispatch.m
    S = struct();
    S.Vprev = Vprev;
    S.uprev = uprev;
    S.Qprev = Qprev;
    S.localPrev = localPrev;
    S.lastResid2 = lastResid2;
    S.CrossingFirst = CrossingFirst;
    S.CrossingAny = CrossingAny;
    S.Qhist = Qhist;
    S.PrevPlanRel = PrevPlanRel;
    [bonModel, B] = realtimeDispatch(t, S, ctx);
    % Restore the build's locals so the solve below is unchanged.
    AubPow = B.AubPow;
    beq = B.beq;
    bub = B.bub;
    headIntercept = B.headIntercept;
    headNow = B.headNow;
    lbNow = B.lbNow;
    mcEmergencyLower = B.mcEmergencyLower;
    nrp = B.nrp;
    obj = B.obj;
    qHighH = B.qHighH;
    qLowH = B.qLowH;
    rpI = B.rpI;
    rpJ = B.rpJ;
    rpK = B.rpK;
    rpV = B.rpV;
    sdCoeff = B.sdCoeff;
    sdConst = B.sdConst;
    ubNow = B.ubNow;
    vIndex=zeros(n,H);
    releaseIndex=zeros(n,H); spillIndex=zeros(n,H); rampIndex=zeros(n,H);
    chanceSlackIndex=zeros(2*n*H,1); slaterRelease=zeros(n*H,1);
    slaterRamp=zeros(n*H,1);
    ss=0;
    for ii=1:n
        p=C.sys(ii); rrDown=abs(p.rr_dn);
        previousSlaterRelease=uprev(ii)*KAF;
        for kk=1:H
            vIndex(ii,kk)=ix(ii,kk,VOL);
            releaseIndex(ii,kk)=ix(ii,kk,REL);
            spillIndex(ii,kk)=ix(ii,kk,SPL);
            rampIndex(ii,kk)=ix(ii,kk,RMP);
            ss=ss+1; chanceSlackIndex(2*ss-1)=ix(ii,kk,SLO);
            chanceSlackIndex(2*ss)=ix(ii,kk,SHI);
            slaterRelease((ii-1)*H+kk)= ...
                max(p.umin,uprev(ii)-kk*rrDown)*KAF;
            slaterRamp((ii-1)*H+kk)=abs( ...
                slaterRelease((ii-1)*H+kk)-previousSlaterRelease);
            previousSlaterRelease=slaterRelease((ii-1)*H+kk);
        end
    end
    sshCtx=struct('n',n,'H',cfg.SSHLookaheadHours, ...
        'Hfull',H,'nvar',nvar,'v_index',vIndex, ...
        'q_low',qLowH,'q_high',chance_upper(qHighH,lowerOnly), ...
        'lower_only',lowerOnly,'sd_const',sdConst, ...
        'sd_coeff',sdCoeff,'sd_floor',1e-10,'Rcorr',Rcorr, ...
        'epsilon',cfg.Epsilon,'phi_tol',cfg.SSHPhiTolerance, ...
        'boundary_tol',cfg.SSHBoundaryTolerance, ...
        'max_iter',cfg.SSHMaxIterations,'max_bisection',60, ...
        'release_index',reshape(releaseIndex',[],1), ...
        'spill_index',reshape(spillIndex',[],1), ...
        'ramp_index',reshape(rampIndex',[],1), ...
        'chance_slack_index',chanceSlackIndex, ...
        'slater_release',slaterRelease,'slater_ramp',slaterRamp, ...
        'prob',probCache,'interior_point',cfg.SSHInteriorPoint, ...
        'interior_margin',1e-6,'NR',NR);

    model=bonModel;
    if strcmp(solutionMethod,'SSH')
        % Remove only the Bonferroni margins. The nominal storage rows,
        % physical bounds, tracking penalties, and emergency variables are
        % unchanged; SSH replaces the margins with joint-probability cuts.
        % Hours beyond SSHLookaheadHours keep their Bonferroni margins,
        % so a shortened SSH look-ahead never leaves later hours with no
        % uncertainty protection at all.
        keepRamp=rpK(1:nrp)>Lh;
        AubSSH=AubPow+sparse(rpI(keepRamp),rpJ(keepRamp),rpV(keepRamp), ...
            nnzRows,nvar);
        bubSSH=bub;
        for ii=1:n
            for kk=1:Lh
                rr=(ii-1)*H+kk; bb=(rr-1)*NR;
                bubSSH(bb+RLO)=-qLowH(kk,ii);
                bubSSH(bb+RHI)=qHighH(kk,ii);
            end
        end
        model=makeModel(Aeq,AubSSH,beq,bubSSH,obj,lbNow,ubNow);
    end
    solveTimer=tic;
    if strcmp(solutionMethod,'SSH')
        [sol,model,sshDiag]=applySSHNative(model,bonModel,gpar,sshCtx);
        if ~strcmp(sol.status,'OPTIMAL')
            % Physically infeasible hour: use the same emergency recourse
            % as BON, on the Bonferroni model, so emergency handling is
            % identical across solution methods.
            model=bonModel;
        end
    else
        sol=gurobi(model,gpar);
        sshDiag=struct('iterations',0,'cuts',0,'used_slater',false, ...
            'lp_solves',1);
    end

    % Emergency fallback: preserve the original hard model whenever it is
    % feasible. If infeasible, first minimize artificial forecast water,
    % then minimize release shortfall at that water correction, and finally
    % optimize the original objective at both minimum emergency quantities.
    if strcmp(sol.status,'INFEASIBLE')
        minSlackIdx=zeros(n*H,1); waterDeficitIdx=zeros(n*H,1); q=0;
        for ii=1:n
            for kk=1:H
                q=q+1;
                minSlackIdx(q)=ix(ii,kk,MINSLK);
                waterDeficitIdx(q)=ix(ii,kk,WDEF);
            end
        end

        emergencyModel=model;
        emergencyModel.ub(minSlackIdx)=inf;
        emergencyModel.ub(waterDeficitIdx)=inf;
        if mcc
            % Releases may now fall below umin: rebuild the McCormick
            % lower-bound rows with the ramp-reachable lower bound.
            eqRows=size(Aeq,1);
            for ii=1:n
                pp=C.sys(ii);
                for kk=1:H
                    rr=(ii-1)*H+kk; bb=(rr-1)*NR; yLe=mcEmergencyLower(ii,kk);
                    emergencyModel.A(eqRows+bb+MC1,ix(ii,kk,VOL))= yLe;
                    emergencyModel.A(eqRows+bb+MC3,ix(ii,kk,VOL))=-yLe;
                    emergencyModel.rhs(eqRows+bb+MC1)=0;
                    emergencyModel.rhs(eqRows+bb+MC3)=-pp.R*yLe;
                end
            end
        end

        % Phase 1: curtail releases as needed before adding forecast water.
        emergencyModel.obj=zeros(nvar,1);
        emergencyModel.obj(waterDeficitIdx)=1;
        phaseWater=gurobi(emergencyModel,gpar);

        if strcmp(phaseWater.status,'OPTIMAL')
            minimumWater=max(sum(phaseWater.x(waterDeficitIdx)),0);
            waterRow=sparse(1,nvar); waterRow(waterDeficitIdx)=1;
            waterTolerance=1e-7*max(1,minimumWater);
            emergencyModel.A=[emergencyModel.A;waterRow];
            emergencyModel.rhs=[emergencyModel.rhs;minimumWater+waterTolerance];
            emergencyModel.sense=[emergencyModel.sense;'<'];

            % Phase 2: minimize minimum-release violations.
            emergencyModel.obj=zeros(nvar,1);
            emergencyModel.obj(minSlackIdx)=1;
            phaseRelease=gurobi(emergencyModel,gpar);

            if strcmp(phaseRelease.status,'OPTIMAL')
                minimumRelease=max(sum(phaseRelease.x(minSlackIdx)),0);
                releaseRow=sparse(1,nvar); releaseRow(minSlackIdx)=1;
                releaseTolerance=1e-7*max(1,minimumRelease);
                emergencyModel.A=[emergencyModel.A;releaseRow];
                emergencyModel.rhs=[emergencyModel.rhs; ...
                    minimumRelease+releaseTolerance];
                emergencyModel.sense=[emergencyModel.sense;'<'];

                % Phase 3: recover the original dispatch objective.
                emergencyModel.obj=obj;
                sol=gurobi(emergencyModel,gpar);

                if strcmp(sol.status,'OPTIMAL')
                    emergencySolveHours=emergencySolveHours+1;
                    EmergencyUsed(t)=true;
                    EmergencyForecastWater(t)=minimumWater;
                    EmergencyForecastRelease(t)=minimumRelease;
                    if emergencySolveHours<=5 || mod(emergencySolveHours,100)==0
                        fprintf(['  emergency relaxation at hour %d: ' ...
                            '24-h forecast water %.4f kaf, release shortfall %.4f kaf\n'], ...
                            t,minimumWater,minimumRelease);
                    end
                end
            end
        end
    end

    if ~strcmp(sol.status,'OPTIMAL')
        % gurobi_iis is a DIAGNOSTIC ONLY and was the one Gurobi call in this
        % file that received no params: no Threads=1, no TimeLimit. Inside a
        % parfor it therefore grabbed every core and could run unbounded,
        % before reaching the error() below -- so the job never reported and
        % the whole pool bogged down with no error file. Gated on Verbose so
        % sweeps (Verbose=false) go straight to the error, and given gpar
        % plus a short TimeLimit when it does run. (2026-09-26)
        if strcmp(sol.status,'INFEASIBLE') && cfg.Verbose
            fprintf('\nIIS for WY%d %s hour %d:\n',C.year,policy,t);
            try
                iisPar=gpar; iisPar.TimeLimit=10;
                iis=gurobi_iis(model,iisPar);
                if isfield(iis,'Arows')
                    badRows=find(iis.Arows);
                    for q=1:numel(badRows)
                        if badRows(q)<=numel(rowNames)
                            fprintf('  constraint: %s\n',rowNames{badRows(q)});
                        else
                            fprintf('  constraint: added row %d\n',badRows(q));
                        end
                    end
                end
                if isfield(iis,'lb')
                    badLb=find(iis.lb);
                    for q=1:numel(badLb)
                        fprintf('  lower bound: %s = %.6g\n', ...
                            varNames{badLb(q)},lb(badLb(q)));
                    end
                end
                if isfield(iis,'ub')
                    badUb=find(iis.ub);
                    for q=1:numel(badUb)
                        fprintf('  upper bound: %s = %.6g\n', ...
                            varNames{badUb(q)},ub(badUb(q)));
                    end
                end
            catch iisError
                fprintf('  IIS diagnostic failed: %s\n',iisError.message);
            end
        end
        error('WY%d %s fixed-head hour %d: %s',C.year,policy,t,sol.status);
    end
    SolveSeconds(t)=toc(solveTimer);
    x=sol.x;
    % ---- step feasibility and same-state probe (diagnostics only)
    lookSlack=0;
    for ii=1:n
        for kk=1:Lh
            lookSlack=lookSlack+x(ix(ii,kk,SLO));
        end
    end
    if strcmp(solutionMethod,'SSH')
        SSHStepFeasible(t)=isfield(sshDiag,'feasible') && sshDiag.feasible;
        SSHGuardUsed(t)=isfield(sshDiag,'guard') && sshDiag.guard;
        probe=gurobi(bonModel,gpar);
        if strcmp(probe.status,'OPTIMAL')
            pSlack=0;
            for ii=1:n
                for kk=1:Lh
                    pSlack=pSlack+probe.x(ix(ii,kk,SLO));
                end
            end
            ProbeBonFeasible(t)=pSlack<=1e-6 && ~CrossingAny(t);
            p1=0; for ii=1:n, p1=p1+probe.x(ix(ii,1,SLO)); end
            BonFirstFeasible(t)=p1<=1e-6 && ~CrossingFirst(t);
            ProbeBonObj(t)=obj'*probe.x;
        else
            ProbeBonFeasible(t)=false; BonFirstFeasible(t)=false;
        end
        ProbeSshObj(t)=obj'*x;
        BonStepFeasible(t)=ProbeBonFeasible(t);
    else
        BonStepFeasible(t)=lookSlack<=1e-6 && ~CrossingAny(t);
        p1=0; for ii=1:n, p1=p1+x(ix(ii,1,SLO)); end
        BonFirstFeasible(t)=p1<=1e-6 && ~CrossingFirst(t);
        ProbeBonFeasible(t)=BonStepFeasible(t);
        ProbeBonObj(t)=obj'*x;
    end
    % Store this plan's releases to predict next hour's look-ahead ramps.
    for ii=1:n
        for kk=1:H
            PrevPlanRel(ii,kk)=(x(ix(ii,kk,REL))+x(ix(ii,kk,SPL)))/KAF;
        end
    end
    % First-look-ahead diagnostics correspond to the decision implemented
    % at this rolling hour; only that hour is recorded, so only it is
    % evaluated (exactly, with analytic sensitivities).
    [phi1,gV1,~,strain1]=sshProbability(x,sshCtx,1);
    weights=abs(gV1(:));
    if sum(weights)>1e-14
        alpha1=weights'/sum(weights);
    else
        alpha1=ones(1,n)/n;
    end
    JointReliability(t)=phi1;
    if strcmp(solutionMethod,'SSH')
        SSHFirstFeasible(t)=phi1>=1-cfg.Epsilon-cfg.SSHPhiTolerance;
    end
    RiskAllocation(t,:)=alpha1;
    SSHStrain(t,:)=strain1(:)';
    SSHIterations(t)=sshDiag.iterations;
    SSHCuts(t)=sshDiag.cuts;
    SSHUsedSlater(t)=sshDiag.used_slater;
    SSHLPSolves(t)=sshDiag.lp_solves;
    % SSH step-infeasible: 1-eps unreachable at some look-ahead hour. SSH
    % then uses the Bonferroni recourse solution for this hour (see
    % applySSHNative); BON is necessarily infeasible there too.
    SSHTargetUnattainable(t)=strcmp(solutionMethod,'SSH') && ~SSHStepFeasible(t);
    if SSHTargetUnattainable(t)
        unattainableCount=unattainableCount+1;
    end
    if isfield(sshDiag,'slater_mode'), SSHMode{t}=sshDiag.slater_mode; end

    u0=zeros(1,n); s0=zeros(1,n); p0=zeros(1,n);
    for i=1:n
        u0(i)=x(ix(i,1,REL));
        s0(i)=x(ix(i,1,SPL));
        p0(i)=x(ix(i,1,POW));
        if mcc
            pm=C.sys(i);
            exactPlan=min(CPOW*max(headIntercept(i)+pm.fb_slope*x(ix(i,1,VOL)),1) ...
                *u0(i)/KAF,pm.pmax);
            McGap(t,i)=p0(i)-exactPlan;
        end
    end
    firstSlack=0;
    for i=1:n
        PlannedLowSlack(t,i)=x(ix(i,1,SLO));
        PlannedHighSlack(t,i)=x(ix(i,1,SHI));
        firstSlack=firstSlack+PlannedLowSlack(t,i)+PlannedHighSlack(t,i);
    end
    slackHours=slackHours+(firstSlack>1e-6);

    Qnow=(u0+s0)/KAF; uImpl=u0; sImpl=s0; pImpl=p0;
    for i=1:n
        j=C.upstream(i);
        if i==1
            routed=C.routedExo(t);
        elseif C.lag(i)==0
            routed=Qnow(j);
        else
            tt=t0-C.lag(i);
            if tt>=0 && isfinite(Qhist(tt+1,j))
                routed=Qhist(tt+1,j);
            else
                routed=Qprev(j);
            end
        end
        if i==1
            a=max(t0-C.lag(i),0)+1;
            b=max(t0-C.lag(i)-1,0)+1;
            du=abs(C.upReleaseExo(a)-C.upReleaseExo(b));
            uLevW=C.upReleaseExo(a);
        elseif C.lag(i)==0
            du=abs(Qnow(j)-known_out(j,t0-1));
            uLevW=Qnow(j);
        else
            tt=t0-C.lag(i);
            du=abs(known_out(j,tt)-known_out(j,tt-1));
            uLevW=known_out(j,tt);
        end

        % Keep the testing environment fixed across kappa values. Kappa
        % changes the uncertainty represented by the DDU optimizer, not
        % the common realized shock path used to compare policies.
        if useArchX
            v=AX.omega(i)+AX.alpha(i)*lastResid2(i) ...
                +AX.gu(i)*max((uLevW-AX.umin(i))/AX.urng(i),0) ...
                +AX.gr(i)*max((du-AX.rmin(i))/AX.rrng(i),0);
            sigmaWorld=min(max(sqrt(max(v,1e-12)),0.2*AX.sbar(i)),3.0*AX.sbar(i));
        else
            sigmaWorld=max(sigDIUdata(i)+gammaData(i)*(du-meanDu(i)), ...
                0.2*sigDIUdata(i));
        end
        e=cfg.WorldKappa*sigmaWorld*C.zseq(t,i);
        rawLocal=C.ar(i,1)+C.ar(i,2)*localPrev(i)+e-cumDev(i)/168;
        lastResid2(i)=e^2;   % seed for the ARCH term at the next solve
        minimumFeasibleLocal=-routed-Vprev(i)/KAF;
        localReal=max(rawLocal,minimumFeasibleLocal);
        projectionKaf=(localReal-rawLocal)*KAF;
        projectedInflowKaf=projectedInflowKaf+projectionKaf;
        projectedInflowHours=projectedInflowHours+(projectionKaf>1e-8);
        cumDev(i)=cumDev(i)+localReal-C.local(t,i);
        LocalReal(t,i)=localReal;

        Vi=Vprev(i)+(localReal+routed)*KAF-uImpl(i)-sImpl(i);
        p=C.sys(i);
        uFloor=max(p.umin,uprev(i)-abs(p.rr_dn))*KAF;   % rr_dn is stored NEGATIVE;
        % without abs() this floor lands ABOVE the current release and the
        % P20 curtailment below can never fire. Every other use of rr_dn in
        % this file already takes abs() (lines ~474, ~507, ~591).
        uCeil=min(p.umax,uprev(i)+p.rr_up)*KAF;

        need=max(C.qLow(t,i)-Vi,0);
        cut=min(need,max(uImpl(i)-uFloor,0));
        uImpl(i)=uImpl(i)-cut; Vi=Vi+cut;

        excess=max(Vi-C.qHigh(t,i),0);
        add=min(excess,max(uCeil-uImpl(i),0));
        uImpl(i)=uImpl(i)+add; Vi=Vi-add; excess=excess-add;
        sImpl(i)=sImpl(i)+excess; Vi=Vi-excess;

        high=max(Vi-C.R(i),0);
        sImpl(i)=sImpl(i)+high; Vi=Vi-high;

        deficit=max(-Vi,0);
        spillCut=min(deficit,sImpl(i));
        sImpl(i)=sImpl(i)-spillCut; Vi=Vi+spillCut; deficit=deficit-spillCut;
        releaseCut=min(deficit,uImpl(i));
        uImpl(i)=uImpl(i)-releaseCut; Vi=Vi+releaseCut;

        remainingDeficit=max(-Vi,0);
        if remainingDeficit>1e-8
            error('Unresolved water deficit at hour %d, unit %s: %.4f kaf', ...
                t,C.units{i},remainingDeficit);
        end

        V(t,i)=min(max(Vi,0),C.R(i));
        if mcc
            hUse=max(headIntercept(i)+p.fb_slope*V(t,i),1);
        else
            hUse=headNow(i);
        end
        Head(t,i)=hUse;
        pImpl(i)=min(CPOW*hUse*max(uImpl(i),0)/KAF,p.pmax);
        curtailMWh=curtailMWh+CPOW*hUse*(cut+releaseCut)/KAF;

        % Record any minimum-release curtailment that was actually
        % implemented, including additional physical recourse after solve.
        MinReleaseShortfall(t,i)=max(p.umin*KAF-uImpl(i),0);
        if MinReleaseShortfall(t,i)>1e-8
            powerAtMinimum=min(CPOW*hUse*p.umin,p.pmax);
            emergencyCurtailMWh=emergencyCurtailMWh+ ...
                max(powerAtMinimum-pImpl(i),0);
        end
        Qnow(i)=(uImpl(i)+sImpl(i))/KAF;
    end

    U(t,:)=uImpl/KAF; Spill(t,:)=sImpl/KAF;
    Power(t,:)=pImpl; Qhist(t,:)=Qnow;
    Vprev=V(t,:); Qprev=Qhist(t,:); localPrev=LocalReal(t,:);
    % The next rolling solve must enforce ramping from the release that was
    % physically implemented after recourse, not the pre-recourse plan.
    uprev=uImpl/KAF;
    if cfg.Verbose && (mod(t,cfg.ProgressHours)==0 || t==T)
        fprintf('  %5d/%d h, energy %.3f TWh, %.1f s elapsed\n', ...
            t,T,sum(Power(1:t,:),'all')/1e6,toc(runTimer));
    end
end

%% 8. Run-level reporting
fprintf('%s inflow projection: %.3f kaf across %d project-hours.\n', ...
    upper(policy),projectedInflowKaf,projectedInflowHours);
minReleaseShortfallKaf=sum(MinReleaseShortfall,'all');
minReleaseShortfallHours=nnz(MinReleaseShortfall>1e-8);
fprintf(['%s emergency minimum-release curtailment: %d rolling solves; ' ...
    'max 24-h forecast-water correction %.3f kaf; ' ...
    '%.3f kaf across %d implemented project-hours (%.3f MWh).\n'], ...
    upper(policy),emergencySolveHours,max(EmergencyForecastWater), ...
    minReleaseShortfallKaf,minReleaseShortfallHours,emergencyCurtailMWh);

%% 9. Reliability scoring and the summary row
lowerDeficit=max(C.qLow(1:T,:)-V,0);
upperExcess=max(V-C.qHigh(1:T,:),0);
% Ignore solver-scale numerical dust when counting boundary exits.
exitToleranceKaf=1e-8;
lowerExitMask=lowerDeficit>exitToleranceKaf;
upperExitMask=upperExcess>exitToleranceKaf;
exitMask=lowerExitMask | upperExitMask;
solveSorted=sort(SolveSeconds);
solveP95=solveSorted(max(1,ceil(0.95*numel(solveSorted))));
s=struct('year',C.year,'scenario_id',scenarioID, ...
    'scenario_seed',scenarioSeed,'policy',upper(policy), ...
    'solution_method',solutionMethod, ...
    'chance_sides',lower(char(cfg.ChanceSides)), ...
    'epsilon',cfg.Epsilon, ...
    'kappa',cfg.Kappa,'hours',T,'energy_MWh',sum(Power,'all'), ...
    'revenue_usd',sum(sum(Power,2).*lmpReal), ...
    'revenue_sched_usd',sum(sum(Power,2).*lmpSched), ...
    'price_noise_sigma',cfg.PriceNoiseSigma, ...
    'price_noise_rho',cfg.PriceNoiseRho, ...
    'price_mode',lower(char(cfg.PriceMode)), ...
    'price_mean_usd_per_mwh',mean(lmp), ...
    'interval_exit_pct',100*mean(exitMask,'all'), ...
    'lower_exit_pct',100*mean(lowerExitMask,'all'), ...
    'upper_exit_pct',100*mean(upperExitMask,'all'), ...
    'lower_ivi_kaf',sum(lowerDeficit,'all'), ...
    'upper_ivi_kaf',sum(upperExcess,'all'), ...
    'energy_gain_pct',NaN,'exit_change_points',NaN, ...
    'margin_rule',marginRule,'margin_z',z, ...
    'curtail_MWh',curtailMWh,'spill_kaf',sum(Spill,'all')*KAF, ...
    'mean_track_error',mean(abs(V-C.ref(1:T,:))./C.R,'all'), ...
    'planned_slack_hours_pct',100*slackHours/T, ...
    'planned_storage_slack_kaf',sum(PlannedLowSlack+PlannedHighSlack,'all'), ...
    'emergency_solve_hours',emergencySolveHours, ...
    'emergency_forecast_water_max_kaf',max(EmergencyForecastWater), ...
    'emergency_forecast_release_max_kaf',max(EmergencyForecastRelease), ...
    'min_release_shortfall_kaf',minReleaseShortfallKaf, ...
    'min_release_shortfall_project_hours',minReleaseShortfallHours, ...
    'emergency_curtail_MWh',emergencyCurtailMWh, ...
    'inflow_projection_kaf',projectedInflowKaf, ...
    'inflow_projection_project_hours',projectedInflowHours, ...
    'solve_time_mean_s',mean(SolveSeconds), ...
    'solve_time_median_s',median(SolveSeconds), ...
    'solve_time_p95_s',solveP95, ...
    'solve_time_max_s',max(SolveSeconds), ...
    'joint_reliability_mean',mean(JointReliability,'omitnan'), ...
    'joint_reliability_min',min(JointReliability,[],'omitnan'), ...
    'joint_target_violation_hours',nnz(JointReliability < ...
        SSHTargetUsed-cfg.SSHPhiTolerance), ...
    'ssh_target_unattainable_hours',nnz(SSHTargetUnattainable), ...
    'ssh_iterations_mean',mean(SSHIterations), ...
    'ssh_iterations_max',max(SSHIterations), ...
    'ssh_cuts_total',sum(SSHCuts), ...
    'ssh_slater_fallback_hours',nnz(SSHUsedSlater), ...
    'ssh_lp_solves_mean',mean(SSHLPSolves), ...
    'ssh_lookahead_hours',Lh, ...
    'chance_infeasible_steps',nnz(~BonStepFeasible & ~strcmp(solutionMethod,'SSH')) ...
        +nnz(~SSHStepFeasible & strcmp(solutionMethod,'SSH')), ...
    'bon_crossing_steps',nnz(CrossingAny), ...
    'first_bon_infeasible_hours',nnz(~BonFirstFeasible), ...
    'first_ssh_infeasible_hours',nnz(~SSHFirstFeasible), ...
    'first_ssh_only_feasible_hours',nnz(~BonFirstFeasible & SSHFirstFeasible), ...
    'probe_bon_infeasible_steps',nnz(~ProbeBonFeasible), ...
    'probe_ssh_infeasible_steps',nnz(~SSHStepFeasible), ...
    'probe_ssh_only_feasible_steps',nnz(~ProbeBonFeasible & SSHStepFeasible), ...
    'probe_containment_violations',nnz(ProbeBonFeasible & ~SSHStepFeasible), ...
    'probe_obj_worse_steps',nnz(ProbeBonFeasible & ...
        ProbeSshObj>ProbeBonObj+1e-6*max(1,abs(ProbeBonObj))), ...
    'probe_obj_gain_pct_mean',100*mean((ProbeBonObj(ProbeBonFeasible)- ...
        ProbeSshObj(ProbeBonFeasible))./max(abs(ProbeBonObj(ProbeBonFeasible)),1),'omitnan'), ...
    'ssh_guard_steps',nnz(SSHGuardUsed), ...
    'wall_seconds',toc(runTimer), ...
    'local_inflow_kaf',sum(LocalReal,'all')*KAF, ...
    'historical_local_inflow_kaf',sum(C.local(1:T,:),'all')*KAF, ...
    'power_model',lower(char(cfg.PowerModel)), ...
    'mc_gap_mean_MW',mean(max(McGap,0),'all'),'mc_gap_max_MW',max(McGap,[],'all'));
%% 10. Return
out.summary=s;
out.trajectory=struct('V',V,'U',U,'spill',Spill,'power',Power, ...
    'head',Head,'qLow',C.qLow(1:T,:),'qHigh',C.qHigh(1:T,:), ...
    'ref',C.ref(1:T,:),'min_release_shortfall_kaf', ...
    MinReleaseShortfall,'emergency_used',EmergencyUsed, ...
    'emergency_forecast_water_kaf',EmergencyForecastWater, ...
    'emergency_forecast_release_kaf',EmergencyForecastRelease, ...
    'local_real',LocalReal,'solve_seconds',SolveSeconds, ...
    'joint_reliability',JointReliability, ...
    'risk_allocation',RiskAllocation, ...
    'risk_budget',cfg.Epsilon*RiskAllocation, ...
    'ssh_strain_sigma',SSHStrain, ...
    'ssh_iterations',SSHIterations,'ssh_cuts',SSHCuts, ...
    'ssh_used_slater',SSHUsedSlater,'ssh_lp_solves',SSHLPSolves, ...
    'ssh_target_unattainable',SSHTargetUnattainable,'ssh_mode',{SSHMode}, ...
    'planned_low_slack_kaf',PlannedLowSlack, ...
    'bon_step_feasible',BonStepFeasible,'ssh_step_feasible',SSHStepFeasible, ...
    'probe_bon_obj',ProbeBonObj,'probe_ssh_obj',ProbeSshObj, ...
    'bon_crossing',CrossingAny,'ssh_guard_used',SSHGuardUsed, ...
    'bon_first_feasible',BonFirstFeasible,'ssh_first_feasible',SSHFirstFeasible, ...
    'planned_high_slack_kaf',PlannedHighSlack, ...
    'solution_method',solutionMethod,'epsilon',cfg.Epsilon, ...
    'units',{C.units},'mc_gap_MW',McGap);


    function value=known_out(j,tt0)
        tt0=min(tt0,T-1);
        if tt0>=0 && isfinite(Qhist(tt0+1,j))
            value=Qhist(tt0+1,j);
        else
            value=Qprev(j);
        end
    end
end

function q=chance_upper(qHigh,lowerOnly)
% Upper limit used in the joint reliability: +inf when only the lower bound
% is a chance constraint (the physical maximum is enforced deterministically).
if lowerOnly
    q=inf(size(qHigh));
else
    q=qHigh;
end
end
