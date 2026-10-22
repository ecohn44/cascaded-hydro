%% ========================================================================
% SEQUENTIAL SUPPORTING HYPERPLANES
% ========================================================================

function [sol,model,diagOut] = applySSHNative(model,bonModel,gpar,ctx)
% Solve one rolling LP under a joint chance constraint by supporting cuts
%
% Inputs:
%   model     - rolling LP with the Bonferroni margins removed for the SSH
%               look-ahead hours (outer approximation of the joint set)
%   bonModel  - the same rolling LP with the Bonferroni margins
%   gpar      - Gurobi parameters
%   ctx       - SSH settings: epsilon, H, max_iter, phi_tol, the index
%               arrays, and the interior-point rule
%
% Outputs:
%   sol       - the returned solution, Gurobi result form
%   model     - the model with the accepted cuts appended
%   diagOut   - per-hour diagnostics: iterations, cuts, reliability,
%               feasible, guard, used_slater, slater_mode
%
% Interior point. SSH needs a point strictly inside the joint chance region.
% Candidates, in order: the Bonferroni optimum, the max-centered Bonferroni
% point (findSlater), and, if neither reaches the target, the
% max-centered point of the SSH model itself (storage margin allowed to be
% negative, so it always exists).
%
% Infeasible step. If no candidate reaches 1-epsilon at every look-ahead
% hour, the joint chance set is empty from this state (the Bonferroni set,
% which it contains, is then empty too). SSH does not relax its target: it
% returns the Bonferroni recourse solution (penalized storage slack), the
% same recourse BON uses, and reports feasible=false.
%
% Dominance guard. Sigma is fixed within the step, so the Bonferroni
% optimum satisfies every valid cut and SSH's objective cannot exceed it.
% If numerics ever produce a worse SSH objective, the Bonferroni optimum is
% returned and guard=true is reported (expected count: 0).
%
% Speed. Reliability for all look-ahead hours is one vectorized call,
% gradients are evaluated only at boundary points, cuts are appended in one
% block per iteration, and re-solves are warm started (dual simplex).
%
% Fallback. If SSH does not converge within ctx.max_iter, the returned
% solution is the better of the Bonferroni optimum (if it meets the
% targets) and the point on the segment toward the interior point where
% every hour meets its target. The hour is flagged (used_slater=true).

%% 1. Settings and the first deterministic LP
target=1-ctx.epsilon;
if ~isfield(ctx,'interior_point'), ctx.interior_point='bonferroni'; end
if ~isfield(ctx,'interior_margin'), ctx.interior_margin=1e-6; end
H=ctx.H;
lpSolves=0;
bon=[]; phiBon=[];
targetVec=target*ones(H,1);

sol=gurobi(model,gpar); lpSolves=lpSolves+1;
if ~strcmp(sol.status,'OPTIMAL')
    % Let realtimeSimulation apply its emergency recourse (as it does for BON).
    diagOut=make_diag(nan(H,1),targetVec,0,0,false,false,zeros(0,1), ...
        nan(H,1),sprintf('base SSH LP %s',sol.status),lpSolves,false);
    diagOut.feasible=false;
    return
end
phiNow=phi_all(sol.x,ctx);
if min(phiNow)>=target-ctx.phi_tol
    diagOut=make_diag(phiNow,targetVec,0,0,false,false,zeros(0,1), ...
        nan(H,1),'not required',lpSolves,true);
    diagOut.feasible=true;
    return
end

%% 2. Interior point
cand={}; candPhi={}; candMode={};
if strcmpi(ctx.interior_point,'bonferroni')
    [bon,phiBon,lpSolves]=solve_bonferroni(bonModel,gpar,ctx,lpSolves);
    if ~isempty(bon)
        cand{end+1}=bon.x; candPhi{end+1}=phiBon;
        candMode{end+1}='Bonferroni optimum';
        if min(phiBon)>=target+ctx.interior_margin
            cand=cand(end); candPhi=candPhi(end); candMode=candMode(end);
        end
    end
end
if isempty(cand) || min(candPhi{1})<target+ctx.interior_margin
    try
        [slater,slaterMode]=findSlater(bonModel,gpar,ctx);
        lpSolves=lpSolves+1;
        cand{end+1}=slater.x; candPhi{end+1}=phi_all(slater.x,ctx);
        candMode{end+1}=slaterMode;
    catch
        % Bonferroni-centered point does not exist; handled below.
    end
    if max(cellfun(@min,[candPhi {-inf}]))<target-ctx.phi_tol
        [xc,etaC,ok]=center_point(model,gpar,ctx);
        lpSolves=lpSolves+1;
        if ok
            cand{end+1}=xc; candPhi{end+1}=phi_all(xc,ctx);
            candMode{end+1}=sprintf('max-centered SSH point, eta %.4g kaf',etaC);
        end
    end
end
if isempty(cand)
    [sol,diagOut,lpSolves]=fallback_bonferroni(bon,phiBon,bonModel,gpar, ...
        ctx,lpSolves,targetVec,'infeasible: no interior candidate');
    return
end
% Interior = candidate with the highest worst-hour reliability.
[~,best]=max(cellfun(@min,candPhi));
interior=cand{best}; interiorPhi=candPhi{best}; interiorMode=candMode{best};

% Joint chance set empty from this state: no candidate reaches 1-eps at
% every look-ahead hour. Use the Bonferroni recourse solution.
if min(interiorPhi)<target+ctx.interior_margin
    [sol,diagOut,lpSolves]=fallback_bonferroni(bon,phiBon,bonModel,gpar, ...
        ctx,lpSolves,targetVec,sprintf('infeasible: best interior %.4f (%s)', ...
        min(interiorPhi),interiorMode));
    return
end
isRelaxed=false;

%% 3. Cut loop
cuts=0; converged=false; failed=false;
dxHistory=zeros(ctx.max_iter,1); iterations=0;
for iter=1:ctx.max_iter
    violated=find(phiNow<targetVec-ctx.phi_tol);
    if isempty(violated)
        converged=true;
        break
    end
    xCurrent=sol.x;
    newRows=zeros(numel(violated),ctx.nvar);
    newRhs=zeros(numel(violated),1);
    nNew=0;
    for kk=violated(:)'
        xBoundary=boundary_point(xCurrent,interior,ctx,kk,targetVec(kk));
        [~,gV,gSd]=sshProbability(xBoundary,ctx,kk);
        gradient=zeros(ctx.nvar,1);
        gradient(ctx.v_index(:,kk))=gV;
        rows=(kk-1)*ctx.n+(1:ctx.n);
        gradient=gradient+ctx.sd_coeff(rows,:)'*gSd;
        scale=max(abs(gradient));
        if ~(scale>1e-12) || any(~isfinite(gradient))
            continue
        end
        gradient=gradient/scale;
        % Chance-feasible superlevel set lies in g'(x-xB)>=0; Gurobi uses <=.
        nNew=nNew+1;
        newRows(nNew,:)=-gradient';
        newRhs(nNew)=-gradient'*xBoundary;
    end
    if nNew==0
        failed=true;
        break
    end
    model.A=[model.A; sparse(newRows(1:nNew,:))];
    model.rhs=[model.rhs; newRhs(1:nNew)];
    model.sense=[model.sense; repmat('<',nNew,1)];
    if isfield(model,'constrnames')
        model.constrnames=[model.constrnames; ...
            arrayfun(@(q)sprintf('ssh_cut%04d',q),(cuts+1:cuts+nNew)', ...
            'UniformOutput',false)];
    end
    cuts=cuts+nNew;

    % Warm start the dual simplex from the previous optimal basis. The new
    % cut rows start basic (their slacks), which is dual feasible.
    if isfield(sol,'vbasis') && isfield(sol,'cbasis')
        model.vbasis=sol.vbasis;
        model.cbasis=[sol.cbasis; zeros(nNew,1)];
    end
    next=gurobi(model,gpar); lpSolves=lpSolves+1;
    iterations=iter;
    if ~strcmp(next.status,'OPTIMAL')
        failed=true;
        break
    end
    dxHistory(iter)=norm(next.x-xCurrent,inf);
    sol=next;
    phiNow=phi_all(sol.x,ctx);
end
if ~converged && ~failed && all(phiNow>=targetVec-ctx.phi_tol)
    converged=true;
end
model=strip_basis(model);
dxHistory=dxHistory(1:iterations);

if converged
    diagOut=make_diag(phiNow,targetVec,iterations,cuts,false,isRelaxed, ...
        dxHistory,interiorPhi,interiorMode,lpSolves,true);
    diagOut.feasible=true;
    [sol,diagOut]=dominance_guard(sol,diagOut,bon,phiBon,bonModel,target,ctx);
    return
end

%% 4. Feasible fallback
lambda=0;
for kk=find(phi_all(sol.x,ctx)<targetVec)'
    [~,lamK]=boundary_point(sol.x,interior,ctx,kk,targetVec(kk),true);
    lambda=max(lambda,lamK);
end
xLine=(1-lambda)*sol.x+lambda*interior;
phiLine=phi_all(xLine,ctx);
candidates={}; objectives=[]; phis={}; labels={};
if all(phiLine>=targetVec-ctx.phi_tol)
    candidates{end+1}=xLine; objectives(end+1)=model.obj'*xLine;
    phis{end+1}=phiLine; labels{end+1}='segment point';
end
if isempty(bon)
    [bon,phiBon,lpSolves]=solve_bonferroni(bonModel,gpar,ctx,lpSolves);
end
if ~isempty(bon) && all(phiBon>=targetVec-ctx.phi_tol)
    candidates{end+1}=bon.x; objectives(end+1)=model.obj'*bon.x;
    phis{end+1}=phiBon; labels{end+1}='Bonferroni optimum';
end
if isempty(candidates)
    candidates{1}=interior; objectives=model.obj'*interior;
    phis{1}=interiorPhi; labels{1}='interior point';
end
[~,best]=min(objectives);
sol=struct('status','OPTIMAL','x',candidates{best}, ...
    'objval',objectives(best));
diagOut=make_diag(phis{best},targetVec,iterations,cuts,true,isRelaxed, ...
    dxHistory,interiorPhi,sprintf('SSH not converged; used %s (%s)', ...
    labels{best},interiorMode),lpSolves,false);
diagOut.feasible=all(phis{best}>=targetVec-ctx.phi_tol);
[sol,diagOut]=dominance_guard(sol,diagOut,bon,phiBon,bonModel,target,ctx);
end

function [sol,d]=dominance_guard(sol,d,bon,phiBon,bonModel,target,ctx)
% Return the Bonferroni optimum if (numerically) it beats SSH. Only a
% Bonferroni plan that meets 1-eps everywhere is a valid replacement.
d.guard=false;
if isempty(bon) || isempty(phiBon) || any(phiBon<target-ctx.phi_tol)
    return
end
objSSH=bonModel.obj'*sol.x; objBON=bonModel.obj'*bon.x;
if objSSH>objBON+1e-6*max(1,abs(objBON))
    sol=bon; d.guard=true; d.phi=phiBon;
end
end

% =====================================================================
function phi=phi_all(x,ctx)
phi=sshProbability(x,ctx,1:ctx.H);
end

function [x,eta,ok]=center_point(model,gpar,ctx)
% Maximize the smallest storage margin (kaf) in the SSH model
% The margin may be negative, so this point exists whenever the SSH LP is
% feasible. Planned chance slack is removed from these auxiliary rows.
NR=10; if isfield(ctx,'NR'), NR=ctx.NR; end; RLO=5; RHI=6; n=ctx.n; Hf=ctx.Hfull;
rowsIdx=zeros(2*n*Hf,1); q=0;
for i=1:n
    for k=1:Hf
        rr=(i-1)*Hf+k;
        q=q+1; rowsIdx(q)=n*Hf+(rr-1)*NR+RLO;
        q=q+1; rowsIdx(q)=n*Hf+(rr-1)*NR+RHI;
    end
end
if isfield(ctx,'lower_only') && ctx.lower_only
    rowsIdx=rowsIdx(1:2:end);   % margin to the P20 floor only
end
c=strip_basis(model);
nv=numel(c.obj);
newRows=c.A(rowsIdx,:);
newRows(:,ctx.chance_slack_index)=0;
c.A=[c.A sparse(size(c.A,1),1); newRows sparse(ones(numel(rowsIdx),1))];
c.rhs=[c.rhs; c.rhs(rowsIdx)];
c.sense=[c.sense; repmat('<',numel(rowsIdx),1)];
c.obj=[zeros(nv,1); -1];
c.lb=[c.lb; -inf]; c.ub=[c.ub; inf];
if isfield(c,'varnames'), c=rmfield(c,'varnames'); end
if isfield(c,'constrnames'), c=rmfield(c,'constrnames'); end
r=gurobi(c,gpar);
ok=strcmp(r.status,'OPTIMAL');
if ok
    x=r.x(1:nv); eta=r.x(end);
else
    x=[]; eta=NaN;
end
end

function [xB,lambda]=boundary_point(x,interior,ctx,kk,target,feasibleSide)
% Bisection on the segment x -> interior for hour kk
% By default returns the point just on the INFEASIBLE side (phi<=target),
% which keeps the supporting hyperplane a valid outer cut.
if nargin<6, feasibleSide=false; end
lo=0; hi=1;
for bis=1:ctx.max_bisection
    lambda=0.5*(lo+hi);
    phiB=sshProbability((1-lambda)*x+lambda*interior,ctx,kk);
    if abs(phiB-target)<=ctx.boundary_tol
        lo=lambda; hi=lambda;
        break
    end
    if phiB<target, lo=lambda; else, hi=lambda; end
end
if feasibleSide, lambda=hi; else, lambda=lo; end
xB=(1-lambda)*x+lambda*interior;
end

function [bon,phiBon,lpSolves]=solve_bonferroni(bonModel,gpar,ctx,lpSolves)
bon=gurobi(bonModel,gpar); lpSolves=lpSolves+1;
if strcmp(bon.status,'OPTIMAL')
    phiBon=phi_all(bon.x,ctx);
else
    bon=[]; phiBon=[];
end
end

function [sol,diagOut,lpSolves]=fallback_bonferroni(bon,phiBon,bonModel, ...
    gpar,ctx,lpSolves,targetVec,mode)
if isempty(bon)
    [bon,phiBon,lpSolves]=solve_bonferroni(bonModel,gpar,ctx,lpSolves);
end
if isempty(bon)
    % Physically infeasible hour: realtimeSimulation applies its emergency
    % recourse on the Bonferroni model, exactly as for BON.
    sol=struct('status','INFEASIBLE','x',[]);
    diagOut=make_diag(nan(ctx.H,1),targetVec,0,0,false,false,zeros(0,1), ...
        nan(ctx.H,1),[mode '; Bonferroni LP infeasible'],lpSolves,false);
    return
end
sol=bon;
% Infeasible step (not a convergence failure): used_slater stays false;
% the step is reported through feasible=false.
diagOut=make_diag(phiBon,targetVec,0,0,false,false,zeros(0,1),phiBon, ...
    mode,lpSolves,false);
diagOut.feasible=false;
end

function model=strip_basis(model)
if isfield(model,'vbasis'), model=rmfield(model,'vbasis'); end
if isfield(model,'cbasis'), model=rmfield(model,'cbasis'); end
end

function d=make_diag(phi,targetVec,iterations,cuts,usedFallback,relaxed, ...
    dx,slaterPhi,mode,lpSolves,converged)
d=struct('phi',phi,'target',targetVec,'relaxed',relaxed, ...
    'alpha',[],'strain',[],'iterations',iterations, ...
    'cuts',cuts,'used_slater',usedFallback,'dx_history',dx, ...
    'slater_phi',slaterPhi,'slater_mode',mode,'lp_solves',lpSolves, ...
    'converged',converged,'feasible',false,'guard',false);
end
